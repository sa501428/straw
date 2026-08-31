/*
  The MIT License (MIT)

  Copyright (c) 2017-2021 Aiden Lab, Rice University, Baylor College of Medicine

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <cmath>
#include <set>
#include <utility>
#include <vector>
#include <streambuf>
#include <curl/curl.h>
#include <iterator>
#include <algorithm>
#include "zlib.h"
#include "zstd.h"
#include "straw.h"
#include "straw_v10.h"
#include <thread>
#include <mutex>
#include <future>
#include <queue>
#include <condition_variable>
#include "hic_slice.h"

using namespace std;

StrawException::StrawException(StrawErrorCode code, const string &message)
    : runtime_error(message), errorCode(code) {}

StrawErrorCode StrawException::code() const noexcept { return errorCode; }

/*
  Straw: fast C++ implementation of dump. Not as fully featured as the
  Java version. Reads the .hic file, finds the appropriate matrix and slice
  of data, and outputs as text in sparse upper triangular format.

  Currently only supporting matrices.

  Usage: straw [observed/oe/expected] <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>
 */

// callback for libcurl. data written to this buffer
static size_t WriteMemoryCallback(void *contents, size_t size, size_t nmemb, void *userp) {
    size_t realsize = size * nmemb;
    struct MemoryStruct *mem;
    mem = (struct MemoryStruct *) userp;

    char *resized = static_cast<char *>(realloc(mem->memory, mem->size + realsize + 1));
    if (resized == nullptr) {
        /* out of memory! */
        printf("not enough memory (realloc returned NULL)\n");
        return 0;
    }
    mem->memory = resized;

    std::memcpy(&(mem->memory[mem->size]), contents, realsize);
    mem->size += realsize;
    mem->memory[mem->size] = 0;

    return realsize;
}

// Get a buffer that can be used as an input stream from the URL.
//
// `received` reports how many bytes actually arrived, which is not always the
// number requested: several callers deliberately ask for a fixed window that may
// extend past EOF. Callers must bound their memstream by `received`, not by the
// requested size, or they read past the end of the allocation.
char *getData(CURL *curl, int64_t position, int64_t chunksize, int64_t &received) {
    received = 0;
    if (chunksize <= 0) {
        throw runtime_error("Invalid HTTP range length requested");
    }
    std::ostringstream oss;
    struct MemoryStruct chunk{};
    chunk.memory = static_cast<char *>(malloc(1));
    if (chunk.memory == nullptr) {
        throw bad_alloc();
    }
    chunk.memory[0] = 0;
    chunk.size = 0;    /* no data at this point */
    oss << position << "-" << position + chunksize - 1;
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void *) &chunk);
    curl_easy_setopt(curl, CURLOPT_RANGE, oss.str().c_str());
    CURLcode res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
        const string message = string("curl_easy_perform() failed: ") + curl_easy_strerror(res);
        free(chunk.memory);
        throw runtime_error(message);
    }

    long status = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &status);
    // A server that ignores Range answers 200 with the whole file. Accepting that
    // would silently mis-parse every offset and re-download the file per request.
    if (status != 206 && static_cast<int64_t>(chunk.size) > chunksize) {
        free(chunk.memory);
        throw runtime_error("Server ignored the HTTP range request (status " +
                            to_string(status) + "); .hic over HTTP requires range support");
    }
    if (static_cast<int64_t>(chunk.size) > chunksize) {
        free(chunk.memory);
        throw runtime_error("Server returned more data than the requested HTTP range");
    }
    if (chunk.size == 0) {
        free(chunk.memory);
        throw runtime_error("Server returned an empty HTTP range response");
    }

    received = static_cast<int64_t>(chunk.size);
    return chunk.memory;
}

bool readMagicString(istream &fin) {
    string str;
    getline(fin, str, '\0');
    return str.size() >= 3 && str[0] == 'H' && str[1] == 'I' && str[2] == 'C';
}

// A short read leaves the destination indeterminate, and no caller inspects the
// stream state. Zero-initializing keeps a truncated or corrupt file from feeding
// uninitialized record counts and bin offsets into the block parser.
char readCharFromFile(istream &fin) {
    char tempChar = 0;
    fin.read(&tempChar, sizeof(char));
    if (!fin) return 0;
    return tempChar;
}

int16_t readInt16FromFile(istream &fin) {
    int16_t tempInt16 = 0;
    fin.read((char *) &tempInt16, sizeof(int16_t));
    if (!fin) return 0;
    return tempInt16;
}

int32_t readInt32FromFile(istream &fin) {
    int32_t tempInt32 = 0;
    fin.read((char *) &tempInt32, sizeof(int32_t));
    if (!fin) return 0;
    return tempInt32;
}

int64_t readInt64FromFile(istream &fin) {
    int64_t tempInt64 = 0;
    fin.read((char *) &tempInt64, sizeof(int64_t));
    if (!fin) return 0;
    return tempInt64;
}

float readFloatFromFile(istream &fin) {
    float tempFloat = 0.0f;
    fin.read((char *) &tempFloat, sizeof(float));
    if (!fin) return 0.0f;
    return tempFloat;
}

double readDoubleFromFile(istream &fin) {
    double tempDouble = 0.0;
    fin.read((char *) &tempDouble, sizeof(double));
    if (!fin) return 0.0;
    return tempDouble;
}

void convertGenomeToBinPos(const int64_t origRegionIndices[4], int64_t regionIndices[4], int32_t resolution) {
    for (uint16_t q = 0; q < 4; q++) {
        // used to find the blocks we need to access
        regionIndices[q] = origRegionIndices[q] / resolution;
    }
}

static CURL *initCURL(const char *url) {
    static once_flag curlInitFlag;
    call_once(curlInitFlag, [] {
        const CURLcode result = curl_global_init(CURL_GLOBAL_DEFAULT);
        if (result != CURLE_OK) {
            throw runtime_error(string("Unable to initialize libcurl: ") + curl_easy_strerror(result));
        }
    });

    CURL *curl = curl_easy_init();
    if (curl) {
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteMemoryCallback);
        curl_easy_setopt(curl, CURLOPT_URL, url);
        curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
        curl_easy_setopt(curl, CURLOPT_USERAGENT, "straw");
    } else {
        throw runtime_error("Unable to initialize curl");
    }
    return curl;
}

class HiCFileStream {
public:
    string prefix = "http"; // HTTP code
    ifstream fin;
    CURL *curl = nullptr;
    bool isHttp = false;

    explicit HiCFileStream(const string &fileName) {
        if (std::strncmp(fileName.c_str(), prefix.c_str(), prefix.size()) == 0) {
            isHttp = true;
            curl = initCURL(fileName.c_str());
            if (!curl) {
                throw runtime_error("URL " + fileName + " cannot be opened for reading");
            }
        } else {
            fin.open(fileName, fstream::in | fstream::binary);
            if (!fin) {
                throw runtime_error("File " + fileName + " cannot be opened for reading");
            }
        }
    }

    ~HiCFileStream() {
        close();
    }

    void close() {
        if (isHttp) {
            if (curl != nullptr) {
                curl_easy_cleanup(curl);
                curl = nullptr;
            }
        } else if (fin.is_open()) {
            fin.close();
        }
    }

    char *readCompressedBytes(indexEntry idx) {
        if (isHttp) {
            int64_t received = 0;
            char *buffer = getData(curl, idx.position, idx.size, received);
            if (received != idx.size) {
                free(buffer);
                throw runtime_error("Short HTTP read for .hic block");
            }
            return buffer;
        } else {
            char *buffer = static_cast<char *>(malloc(static_cast<size_t>(idx.size)));
            if (buffer == nullptr) {
                throw bad_alloc();
            }
            fin.seekg(idx.position, ios::beg);
            fin.read(buffer, idx.size);
            if (!fin) {
                free(buffer);
                throw runtime_error("Unable to read compressed .hic block");
            }
            return buffer;
        }
    }
};

char *readCompressedBytesFromFile(const string &fileName, indexEntry idx) {
    HiCFileStream stream(fileName);
    return stream.readCompressedBytes(idx);
}

// reads the header, storing the positions of the normalization vectors and returning the masterIndexPosition pointer
map<string, chromosome> readHeader(istream &fin, int64_t &masterIndexPosition, string &genomeID,
                                   int32_t &numChromosomes, int32_t &version, int64_t &nviPosition,
                                   int64_t &nviLength,
                                   vector<pair<string, string>> *attributes) {
    map<string, chromosome> chromosomeMap;
    if (!readMagicString(fin)) {
        cerr << "Hi-C magic string is missing, does not appear to be a hic file" << endl;
        masterIndexPosition = -1;
        return chromosomeMap;
    }

    version = readInt32FromFile(fin);
    if (version == 10) {
        map<string, chromosome> chromosomes;
        straw_v10::readHeader(fin, masterIndexPosition, genomeID, numChromosomes, nviPosition, nviLength, chromosomes);
        return chromosomes;
    }
    if (version < 6) {
        cerr << "Version " << version << " no longer supported" << endl;
        masterIndexPosition = -1;
        return chromosomeMap;
    }
    masterIndexPosition = readInt64FromFile(fin);
    getline(fin, genomeID, '\0');

    if (version > 8) {
        nviPosition = readInt64FromFile(fin);
        nviLength = readInt64FromFile(fin);
    }

    int32_t nattributes = readInt32FromFile(fin);

    // reading and ignoring attribute-value dictionary
    for (int i = 0; i < nattributes; i++) {
        string key, value;
        getline(fin, key, '\0');
        getline(fin, value, '\0');
        if (attributes) attributes->emplace_back(key, value);
    }

    numChromosomes = readInt32FromFile(fin);
    // chromosome map for finding matrixType
    for (int i = 0; i < numChromosomes; i++) {
        string name;
        int64_t length;
        getline(fin, name, '\0');
        if (version > 8) {
            length = readInt64FromFile(fin);
        } else {
            length = (int64_t) readInt32FromFile(fin);
        }

        chromosome chr;
        chr.index = i;
        chr.name = name;
        chr.length = length;
        chromosomeMap[name] = chr;
    }
    return chromosomeMap;
}

void readResolutionsFromHeader(istream &fin, vector<int32_t> &bp, vector<int32_t> &frag) {
    int numBpResolutions = readInt32FromFile(fin);
    for (int i = 0; i < numBpResolutions; i++) {
        bp.push_back(readInt32FromFile(fin));
    }
    int numFragResolutions = readInt32FromFile(fin);
    for (int i = 0; i < numFragResolutions; i++) {
        frag.push_back(readInt32FromFile(fin));
    }
}

vector<int32_t> readResolutionsFromHeader(istream &fin) {
    vector<int32_t> bp, frag;
    readResolutionsFromHeader(fin, bp, frag);
    return bp;
}

//https://www.techiedelight.com/get-slice-sub-vector-from-vector-cpp/
vector<double> sliceVector(vector<double> &v, int64_t m, int64_t n) {
    vector<double> vec;
    copy(v.begin() + m, v.begin() + n + 1, back_inserter(vec));
    return vec;
}

void populateVectorWithFloats(istream &fin, vector<double> &vector, int64_t nValues) {
    for (int j = 0; j < nValues; j++) {
        double v = readFloatFromFile(fin);
        vector.push_back(v);
    }
}

void populateVectorWithDoubles(istream &fin, vector<double> &vector, int64_t nValues) {
    for (int j = 0; j < nValues; j++) {
        double v = readDoubleFromFile(fin);
        vector.push_back(v);
    }
}

int64_t readThroughExpectedVectorURL(CURL *curl, int64_t currentPointer, int32_t version, vector<double> &expectedValues, int64_t nValues,
                               bool store, int32_t /*resolution*/) {
    if (store) {
        int64_t bufferSize = static_cast<int64_t>(nValues) * sizeof(double) + 10000;
        if (version > 8) {
            bufferSize = static_cast<int64_t>(nValues) * sizeof(float) + 10000;
        }
        int64_t received = 0;
        char *buffer = getData(curl, currentPointer, bufferSize, received);
        memstream fin(buffer, received);

        if (version > 8) {
            populateVectorWithFloats(fin, expectedValues, nValues);
        } else {
            populateVectorWithDoubles(fin, expectedValues, nValues);
        }
        free(buffer);
    }

    if (version > 8) {
        return nValues * sizeof(float);
    } else {
        return nValues * sizeof(double);
    }
}

void readThroughExpectedVector(int32_t version, istream &fin, vector<double> &expectedValues, int64_t nValues,
                               bool store, int32_t /*resolution*/) {
    if (store) {
        if (version > 8) {
            populateVectorWithFloats(fin, expectedValues, nValues);
        } else {
            populateVectorWithDoubles(fin, expectedValues, nValues);
        }
    } else if (nValues > 0) {
        if (version > 8) {
            fin.seekg(nValues * sizeof(float), ios_base::cur);
        } else {
            fin.seekg(nValues * sizeof(double), ios_base::cur);
        }
    }
}

int64_t readThroughNormalizationFactorsURL(CURL *curl, int64_t currentPointer, int32_t version, bool store, vector<double> &expectedValues,
                                     int32_t c1, int32_t nNormalizationFactors) {

    if (store) {
        int64_t bufferSize = static_cast<int64_t>(nNormalizationFactors) *
                                 (sizeof(int32_t) + sizeof(double)) + 10000;
        if (version > 8) {
            bufferSize = static_cast<int64_t>(nNormalizationFactors) *
                             (sizeof(int32_t) + sizeof(float)) + 10000;
        }
        int64_t received = 0;
        char *buffer = getData(curl, currentPointer, bufferSize, received);
        memstream fin(buffer, received);

        for (int j = 0; j < nNormalizationFactors; j++) {
            int32_t chrIdx = readInt32FromFile(fin);
            double v;
            if (version > 8) {
                v = readFloatFromFile(fin);
            } else {
                v = readDoubleFromFile(fin);
            }
            if (chrIdx == c1) {
                for (double &expectedValue : expectedValues) {
                    expectedValue = expectedValue / v;
                }
            }
        }
        free(buffer);
    }

    if (version > 8) {
        return nNormalizationFactors * (sizeof(int32_t) + sizeof(float));
    } else {
        return nNormalizationFactors * (sizeof(int32_t) + sizeof(double));
    }
}

void readThroughNormalizationFactors(istream &fin, int32_t version, bool store, vector<double> &expectedValues,
                                     int32_t c1) {
    int32_t nNormalizationFactors = readInt32FromFile(fin);
    if (store) {
        for (int j = 0; j < nNormalizationFactors; j++) {
            int32_t chrIdx = readInt32FromFile(fin);
            double v;
            if (version > 8) {
                v = readFloatFromFile(fin);
            } else {
                v = readDoubleFromFile(fin);
            }
            if (chrIdx == c1) {
                for (double &expectedValue : expectedValues) {
                    expectedValue = expectedValue / v;
                }
            }
        }
    } else if (nNormalizationFactors > 0) {
        if (version > 8) {
            fin.seekg(nNormalizationFactors * (sizeof(int32_t) + sizeof(float)), ios_base::cur);
        } else {
            fin.seekg(nNormalizationFactors * (sizeof(int32_t) + sizeof(double)), ios_base::cur);
        }
    }
}

int64_t readStringFromURL(istream &fin, string &basicString) {
    getline(fin, basicString, '\0');
    return (basicString.length() + 1);
}

// reads the footer from the master pointer location. takes in the chromosomes,
// norm, unit (BP or FRAG) and resolution or binsize, and sets the file
// position of the matrix and the normalization vectors for those chromosomes
// at the given normalization and resolution
bool readFooterURL(CURL *curl, int64_t master, int32_t version, int32_t c1, int32_t c2, const string &matrixType,
                const string &norm, const string &unit, int32_t resolution, int64_t &myFilePos,
                indexEntry &c1NormEntry, indexEntry &c2NormEntry, vector<double> &expectedValues) {

    int64_t currentPointer = master;

    int64_t received = 0;
    char *buffer = getData(curl, currentPointer, 100, received);
    memstream newFin(buffer, received);

    if (version > 8) {
        readInt64FromFile(newFin);
        currentPointer += 8;
    } else {
        readInt32FromFile(newFin);
        currentPointer += 4;
    }

    stringstream ss;
    ss << c1 << "_" << c2;
    string key = ss.str();


    int32_t nEntries = readInt32FromFile(newFin);

    currentPointer += 4;
    free(buffer);

    // 50 bytes/entry is a heuristic for a variable-length key plus 12 fixed
    // bytes. Compute in int64 so a large nEntries cannot overflow into a
    // negative (and therefore malformed) range length.
    int64_t bufferSize0 = max<int64_t>(1024, static_cast<int64_t>(nEntries) * 50);
    buffer = getData(curl, currentPointer, bufferSize0, received);

    memstream newFin2(buffer, received);


    bool found = false;
    for (int i = 0; i < nEntries; i++) {
        string keyStr;

        currentPointer += readStringFromURL(newFin2, keyStr);

        int64_t fpos = readInt64FromFile(newFin2);
        int32_t sizeinbytes = readInt32FromFile(newFin2);

        currentPointer += 12;
        if (keyStr == key) {
            myFilePos = fpos;
            found = true;
        }
    }
    free(buffer);
    if (!found) {
        cerr << "Remote file doesn't have the given chr_chr map " << key << endl;
        return false;
    }

    if ((matrixType == "observed" && norm == "NONE") ||
        ((matrixType == "oe" || matrixType == "expected") && norm == "NONE" && c1 != c2))
        return true; // no need to read norm vector index

    // read in and ignore expected value maps; don't store; reading these to
    // get to norm vector index
    buffer = getData(curl, currentPointer, 100, received);

    memstream newFin3(buffer, received);

    int32_t nExpectedValues = readInt32FromFile(newFin3);

    currentPointer += 4;
    free(buffer);
    for (int i = 0; i < nExpectedValues; i++) {

        buffer = getData(curl, currentPointer, 1000, received);

        memstream newFin4(buffer, received);

        string unit0;
        currentPointer += readStringFromURL(newFin4, unit0);

        int32_t binSize = readInt32FromFile(newFin4);

        currentPointer += 4;

        int64_t nValues;
        if (version > 8) {

            nValues = readInt64FromFile(newFin4);
            currentPointer += 8;
        } else {
            nValues = (int64_t) readInt32FromFile(newFin4);

            currentPointer += 4;
        }

        free(buffer);

        bool store = c1 == c2 && (matrixType == "oe" || matrixType == "expected") && norm == "NONE" && unit0 == unit &&
                     binSize == resolution;

        currentPointer += readThroughExpectedVectorURL(curl, currentPointer, version, expectedValues, nValues, store, resolution);

        buffer = getData(curl, currentPointer, 100, received);

        memstream newFin5(buffer, received);
        int32_t nNormalizationFactors = readInt32FromFile(newFin5);

        currentPointer += 4;
        free(buffer);

        currentPointer += readThroughNormalizationFactorsURL(curl, currentPointer, version, store, expectedValues, c1, nNormalizationFactors);
    }

    if (c1 == c2 && (matrixType == "oe" || matrixType == "expected") && norm == "NONE") {
        if (expectedValues.empty()) {
            cerr << "Remote file did not contain expected values vectors at " << resolution << " " << unit << endl;
            return false;
        }
        return true;
    }

    buffer = getData(curl, currentPointer, 100, received);

    memstream newFin6(buffer, received);
    nExpectedValues = readInt32FromFile(newFin6);

    currentPointer += 4;
    free(buffer);
    for (int i = 0; i < nExpectedValues; i++) {
        buffer = getData(curl, currentPointer, 1000, received);

        memstream newFin7(buffer, received);

        string nType, unit0;
        currentPointer += readStringFromURL(newFin7, nType);
        currentPointer += readStringFromURL(newFin7, unit0);

        int32_t binSize = readInt32FromFile(newFin7);

        currentPointer += 4;

        int64_t nValues;
        if (version > 8) {

            nValues = readInt64FromFile(newFin7);
            currentPointer += 8;
        } else {
            nValues = (int64_t) readInt32FromFile(newFin7);

            currentPointer += 4;
        }
        bool store = c1 == c2 && (matrixType == "oe" || matrixType == "expected") && nType == norm && unit0 == unit &&
                     binSize == resolution;

        free(buffer);

        currentPointer += readThroughExpectedVectorURL(curl, currentPointer, version, expectedValues, nValues, store, resolution);

        buffer = getData(curl, currentPointer, 100, received);

        memstream newFin8(buffer, received);
        int32_t nNormalizationFactors = readInt32FromFile(newFin8);

        currentPointer += 4;
        free(buffer);

        currentPointer += readThroughNormalizationFactorsURL(curl, currentPointer, version, store, expectedValues, c1, nNormalizationFactors);
    }

    if (c1 == c2 && (matrixType == "oe" || matrixType == "expected") && norm != "NONE") {
        if (expectedValues.empty()) {
            cerr << "Remote file did not contain normalized expected values vectors at " << resolution << " " << unit << endl;
            return false;
        }
    }

    buffer = getData(curl, currentPointer, 100, received);

    memstream newFin9(buffer, received);
    nEntries = readInt32FromFile(newFin9);

    currentPointer += 4;
    free(buffer);

    bool found1 = false;
    bool found2 = false;
    int64_t bufferSize2 = max<int64_t>(1024, static_cast<int64_t>(nEntries) * 60);
    buffer = getData(curl, currentPointer, bufferSize2, received);

    memstream newFin10(buffer, received);

    for (int i = 0; i < nEntries; i++) {
        string normtype;
        currentPointer += readStringFromURL(newFin10, normtype);

        int32_t chrIdx = readInt32FromFile(newFin10);
        currentPointer += 4;
        string unit1;
        currentPointer += readStringFromURL(newFin10, unit1);

        int32_t resolution1 = readInt32FromFile(newFin10);
        int64_t filePosition = readInt64FromFile(newFin10);

        currentPointer += 12;

        int64_t sizeInBytes;
        if (version > 8) {

            sizeInBytes = readInt64FromFile(newFin10);
            currentPointer += 8;
        } else {
            sizeInBytes = (int64_t) readInt32FromFile(newFin10);

            currentPointer += 4;
        }

        if (chrIdx == c1 && normtype == norm && unit1 == unit && resolution1 == resolution) {
            c1NormEntry.position = filePosition;
            c1NormEntry.size = sizeInBytes;
            found1 = true;
        }
        if (chrIdx == c2 && normtype == norm && unit1 == unit && resolution1 == resolution) {
            c2NormEntry.position = filePosition;
            c2NormEntry.size = sizeInBytes;
            found2 = true;
        }
    }
    free(buffer);
    if (!found1 || !found2) {
        cerr << "Remote file did not contain " << norm << " normalization vectors for one or both chromosomes at "
             << resolution << " " << unit << endl;
        return false;
    }
    return true;
}

bool readFooter(istream &fin, int64_t master, int32_t version, int32_t c1, int32_t c2, const string &matrixType,
                const string &norm, const string &unit, int32_t resolution, int64_t &myFilePos,
                indexEntry &c1NormEntry, indexEntry &c2NormEntry, vector<double> &expectedValues) {

    if (version > 8) {
        int64_t nBytes = readInt64FromFile(fin);
    } else {
        int32_t nBytes = readInt32FromFile(fin);
    }

    stringstream ss;
    ss << c1 << "_" << c2;
    string key = ss.str();

    int32_t nEntries = readInt32FromFile(fin);
    bool found = false;
    for (int i = 0; i < nEntries; i++) {
        string keyStr;
        getline(fin, keyStr, '\0');
        int64_t fpos = readInt64FromFile(fin);
        int32_t sizeinbytes = readInt32FromFile(fin);
        if (keyStr == key) {
            myFilePos = fpos;
            found = true;
        }
    }
    if (!found) {
        cerr << "File doesn't have the given chr_chr map " << key << endl;
        return false;
    }

    if ((matrixType == "observed" && norm == "NONE") ||
        ((matrixType == "oe" || matrixType == "expected") && norm == "NONE" && c1 != c2))
        return true; // no need to read norm vector index

    // read in and ignore expected value maps; don't store; reading these to
    // get to norm vector index
    int32_t nExpectedValues = readInt32FromFile(fin);
    for (int i = 0; i < nExpectedValues; i++) {
        string unit0;
        getline(fin, unit0, '\0'); //unit
        int32_t binSize = readInt32FromFile(fin);

        int64_t nValues;
        if (version > 8) {
            nValues = readInt64FromFile(fin);
        } else {
            nValues = (int64_t) readInt32FromFile(fin);
        }

        bool store = c1 == c2 && (matrixType == "oe" || matrixType == "expected") && norm == "NONE" && unit0 == unit &&
                     binSize == resolution;
        readThroughExpectedVector(version, fin, expectedValues, nValues, store, resolution);
        readThroughNormalizationFactors(fin, version, store, expectedValues, c1);
    }

    if (c1 == c2 && (matrixType == "oe" || matrixType == "expected") && norm == "NONE") {
        if (expectedValues.empty()) {
            cerr << "File did not contain expected values vectors at " << resolution << " " << unit << endl;
            return false;
        }
        return true;
    }

    nExpectedValues = readInt32FromFile(fin);
    for (int i = 0; i < nExpectedValues; i++) {
        string nType, unit0;
        getline(fin, nType, '\0'); //typeString
        getline(fin, unit0, '\0'); //unit
        int32_t binSize = readInt32FromFile(fin);

        int64_t nValues;
        if (version > 8) {
            nValues = readInt64FromFile(fin);
        } else {
            nValues = (int64_t) readInt32FromFile(fin);
        }
        bool store = c1 == c2 && (matrixType == "oe" || matrixType == "expected") && nType == norm && unit0 == unit &&
                     binSize == resolution;
        readThroughExpectedVector(version, fin, expectedValues, nValues, store, resolution);
        readThroughNormalizationFactors(fin, version, store, expectedValues, c1);
    }

    if (c1 == c2 && (matrixType == "oe" || matrixType == "expected") && norm != "NONE") {
        if (expectedValues.empty()) {
            cerr << "File did not contain normalized expected values vectors at " << resolution << " " << unit << endl;
            return false;
        }
    }

    // Index of normalization vectors
    nEntries = readInt32FromFile(fin);
    bool found1 = false;
    bool found2 = false;
    for (int i = 0; i < nEntries; i++) {
        string normtype;
        getline(fin, normtype, '\0'); //normalization type
        int32_t chrIdx = readInt32FromFile(fin);
        string unit1;
        getline(fin, unit1, '\0'); //unit
        int32_t resolution1 = readInt32FromFile(fin);
        int64_t filePosition = readInt64FromFile(fin);
        int64_t sizeInBytes;
        if (version > 8) {
            sizeInBytes = readInt64FromFile(fin);
        } else {
            sizeInBytes = (int64_t) readInt32FromFile(fin);
        }

        if (chrIdx == c1 && normtype == norm && unit1 == unit && resolution1 == resolution) {
            c1NormEntry.position = filePosition;
            c1NormEntry.size = sizeInBytes;
            found1 = true;
        }
        if (chrIdx == c2 && normtype == norm && unit1 == unit && resolution1 == resolution) {
            c2NormEntry.position = filePosition;
            c2NormEntry.size = sizeInBytes;
            found2 = true;
        }
    }
    if (!found1 || !found2) {
        cerr << "File did not contain " << norm << " normalization vectors for one or both chromosomes at "
             << resolution << " " << unit << endl;
        return false;
    }
    return true;
}

indexEntry readIndexEntry(istream &fin) {
    int64_t filePosition = readInt64FromFile(fin);
    int32_t blockSizeInBytes = readInt32FromFile(fin);
    indexEntry entry = indexEntry();
    entry.size = (int64_t) blockSizeInBytes;
    entry.position = filePosition;
    return entry;
}

void setValuesForMZD(istream &fin, const string &myunit, float &mySumCounts, int32_t &mybinsize,
                     int32_t &myBlockBinCount, int32_t &myBlockColumnCount, bool &found) {
    string unit;
    getline(fin, unit, '\0'); // unit
    readInt32FromFile(fin); // Old "zoom" index -- not used
    float sumCounts = readFloatFromFile(fin); // sumCounts
    readFloatFromFile(fin); // occupiedCellCount
    readFloatFromFile(fin); // stdDev
    readFloatFromFile(fin); // percent95
    int32_t binSize = readInt32FromFile(fin);
    int32_t blockBinCount = readInt32FromFile(fin);
    int32_t blockColumnCount = readInt32FromFile(fin);
    found = false;
    if (myunit == unit && mybinsize == binSize) {
        mySumCounts = sumCounts;
        myBlockBinCount = blockBinCount;
        myBlockColumnCount = blockColumnCount;
        found = true;
    }
}

void populateBlockMap(istream &fin, int32_t nBlocks, map<int32_t, indexEntry> &blockMap) {
    for (int b = 0; b < nBlocks; b++) {
        int32_t blockNumber = readInt32FromFile(fin);
        blockMap[blockNumber] = readIndexEntry(fin);
    }
}

// reads the raw binned contact matrix at specified resolution, setting the block bin count and block column count
map<int32_t, indexEntry> readMatrixZoomData(istream &fin, const string &myunit, int32_t mybinsize, float &mySumCounts,
                                            int32_t &myBlockBinCount, int32_t &myBlockColumnCount, bool &found) {

    map<int32_t, indexEntry> blockMap;
    setValuesForMZD(fin, myunit, mySumCounts, mybinsize, myBlockBinCount, myBlockColumnCount, found);

    int32_t nBlocks = readInt32FromFile(fin);
    if (found) {
        populateBlockMap(fin, nBlocks, blockMap);
    } else {
        fin.seekg(nBlocks * (sizeof(int32_t) + sizeof(int64_t) + sizeof(int32_t)), ios_base::cur);
    }
    return blockMap;
}

// reads the raw binned contact matrix at specified resolution, setting the block bin count and block column count
map<int32_t, indexEntry> readMatrixZoomDataHttp(CURL *curl, int64_t &myFilePosition, const string &myunit,
                                                int32_t mybinsize, float &mySumCounts, int32_t &myBlockBinCount,
                                                int32_t &myBlockColumnCount, bool &found) {
    map<int32_t, indexEntry> blockMap;
    int32_t header_size = 5 * sizeof(int32_t) + 4 * sizeof(float);
    int64_t received = 0;
    char *first = getData(curl, myFilePosition, 1, received);
    if (first[0] == 'B') {
        header_size += 3;
    } else if (first[0] == 'F') {
        header_size += 5;
    } else {
        cerr << "Unit not understood" << endl;
        return blockMap;
    }
    free(first);
    char *buffer = getData(curl, myFilePosition, header_size, received);
    memstream fin(buffer, received);
    setValuesForMZD(fin, myunit, mySumCounts, mybinsize, myBlockBinCount, myBlockColumnCount, found);
    int32_t nBlocks = readInt32FromFile(fin);
    free(buffer);

    if (found) {
        int64_t chunkSize = static_cast<int64_t>(nBlocks) *
                                (sizeof(int32_t) + sizeof(int64_t) + sizeof(int32_t));
        buffer = getData(curl, myFilePosition + header_size, chunkSize, received);
        memstream fin2(buffer, received);
        populateBlockMap(fin2, nBlocks, blockMap);
        free(buffer);
    } else {
        myFilePosition = myFilePosition + header_size
                         + (nBlocks * (sizeof(int32_t) + sizeof(int64_t) + sizeof(int32_t)));
    }
    return blockMap;
}

// goes to the specified file pointer in http and finds the raw contact matrixType at specified resolution, calling readMatrixZoomData.
// sets blockbincount and blockcolumncount
map<int32_t, indexEntry> readMatrixHttp(CURL *curl, int64_t myFilePosition, const string &unit, int32_t resolution,
                                        float &mySumCounts, int32_t &myBlockBinCount, int32_t &myBlockColumnCount,
                                        bool &foundResolution) {
    int32_t size = sizeof(int32_t) * 3;
    int64_t received = 0;
    char *buffer = getData(curl, myFilePosition, size, received);
    memstream bufin(buffer, received);

    int32_t c1 = readInt32FromFile(bufin);
    int32_t c2 = readInt32FromFile(bufin);
    int32_t nRes = readInt32FromFile(bufin);
    int32_t i = 0;
    bool found = false;
    myFilePosition = myFilePosition + size;
    free(buffer);
    map<int32_t, indexEntry> blockMap;

    while (i < nRes && !found) {
        // myFilePosition gets updated within call
        blockMap = readMatrixZoomDataHttp(curl, myFilePosition, unit, resolution, mySumCounts, myBlockBinCount,
                                          myBlockColumnCount, found);
        i++;
    }
    if (!found) {
        cerr << "Error finding block data" << endl;
    }
    foundResolution = found;
    return blockMap;
}

// goes to the specified file pointer and finds the raw contact matrixType at specified resolution, calling readMatrixZoomData.
// sets blockbincount and blockcolumncount
map<int32_t, indexEntry> readMatrix(istream &fin, int64_t myFilePosition, const string &unit, int32_t resolution,
                                    float &mySumCounts, int32_t &myBlockBinCount, int32_t &myBlockColumnCount,
                                    bool &foundResolution) {
    map<int32_t, indexEntry> blockMap;

    fin.seekg(myFilePosition, ios::beg);
    int32_t c1 = readInt32FromFile(fin);
    int32_t c2 = readInt32FromFile(fin);
    int32_t nRes = readInt32FromFile(fin);
    int32_t i = 0;
    bool found = false;
    while (i < nRes && !found) {
        blockMap = readMatrixZoomData(fin, unit, resolution, mySumCounts, myBlockBinCount, myBlockColumnCount, found);
        i++;
    }
    if (!found) {
        cerr << "Error finding block data" << endl;
    }
    foundResolution = found;
    return blockMap;
}

// gets the blocks that need to be read for this slice of the data.  needs blockbincount, blockcolumncount, and whether
// or not this is intrachromosomal.
set<int32_t> getBlockNumbersForRegionFromBinPosition(const int64_t *regionIndices, int32_t blockBinCount,
                                                     int32_t blockColumnCount, bool intra) {
    int32_t col1, col2, row1, row2;
    col1 = static_cast<int32_t>(regionIndices[0] / blockBinCount);
    col2 = static_cast<int32_t>((regionIndices[1] + 1) / blockBinCount);
    row1 = static_cast<int32_t>(regionIndices[2] / blockBinCount);
    row2 = static_cast<int32_t>((regionIndices[3] + 1) / blockBinCount);

    set<int32_t> blocksSet;
    // first check the upper triangular matrixType
    for (int r = row1; r <= row2; r++) {
        for (int c = col1; c <= col2; c++) {
            int32_t blockNumber = r * blockColumnCount + c;
            blocksSet.insert(blockNumber);
        }
    }
    // check region part that overlaps with lower left triangle but only if intrachromosomal
    if (intra) {
        for (int r = col1; r <= col2; r++) {
            for (int c = row1; c <= row2; c++) {
                int32_t blockNumber = r * blockColumnCount + c;
                blocksSet.insert(blockNumber);
            }
        }
    }
    return blocksSet;
}

set<int32_t> getBlockNumbersForRegionFromBinPositionV9Intra(int64_t *regionIndices, int32_t blockBinCount,
                                                            int32_t blockColumnCount) {
    // regionIndices is binX1 binX2 binY1 binY2
    set<int32_t> blocksSet;
    int32_t translatedLowerPAD, translatedHigherPAD, translatedNearerDepth, translatedFurtherDepth;
    translatedLowerPAD = static_cast<int32_t>((regionIndices[0] + regionIndices[2]) / 2 / blockBinCount);
    translatedHigherPAD = static_cast<int32_t>((regionIndices[1] + regionIndices[3]) / 2 / blockBinCount + 1);
    translatedNearerDepth = static_cast<int32_t>(log2(
            1 + abs(regionIndices[0] - regionIndices[3]) / sqrt(2) / blockBinCount));
    translatedFurtherDepth = static_cast<int32_t>(log2(
            1 + abs(regionIndices[1] - regionIndices[2]) / sqrt(2) / blockBinCount));

    // because code above assume above diagonal; but we could be below diagonal
    int32_t nearerDepth = min(translatedNearerDepth, translatedFurtherDepth);
    if ((regionIndices[0] > regionIndices[3] && regionIndices[1] < regionIndices[2]) ||
        (regionIndices[1] > regionIndices[2] && regionIndices[0] < regionIndices[3])) {
        nearerDepth = 0;
    }
    int32_t furtherDepth = max(translatedNearerDepth, translatedFurtherDepth) + 1; // +1; integer divide rounds down

    for (int depth = nearerDepth; depth <= furtherDepth; depth++) {
        for (int pad = translatedLowerPAD; pad <= translatedHigherPAD; pad++) {
            int32_t blockNumber = depth * blockColumnCount + pad;
            blocksSet.insert(blockNumber);
        }
    }

    return blocksSet;
}

// The v7+ row/column loops advance `index` independently of the record count in
// the block header, so a corrupt or truncated block can drive it past the end of
// the vector. Reject that instead of writing out of bounds.
void appendRecord(vector<contactRecord> &vector, int32_t index, int32_t binX, int32_t binY, float counts) {
    if (index < 0 || static_cast<size_t>(index) >= vector.size()) {
        throw runtime_error("Corrupt .hic block: more records than the block header declares");
    }
    contactRecord record = contactRecord();
    record.binX = binX;
    record.binY = binY;
    record.counts = counts;
    vector[index] = record;
}

// Upper bound on how far a stored block expands. Both the output buffer and the
// decompressor's capacity are derived from this single constant.
static const int64_t DECOMPRESSION_FACTOR = 10; // biggest ratio seen so far is 3

// Auto-detect zstd vs zlib from magic bytes: zstd frames start with 0xFD 0x2F 0xB5 0x28
static bool isZstdCompressed(const char *data, int32_t size) {
    return size >= 4 &&
           (unsigned char)data[0] == 0xFD &&
           (unsigned char)data[1] == 0x2F &&
           (unsigned char)data[2] == 0xB5 &&
           (unsigned char)data[3] == 0x28;
}

// Returns the decompressed size, or -1 when the block could not be fully
// decompressed. A partial result must never be reported as success: the record
// parser trusts the length and would read past the end of the payload.
int32_t decompressBlock(indexEntry idx, char *compressedBytes, char *uncompressedBytes) {
    const size_t dstCapacity = static_cast<size_t>(idx.size) * DECOMPRESSION_FACTOR;
    if (isZstdCompressed(compressedBytes, idx.size)) {
        size_t const result = ZSTD_decompress(uncompressedBytes, dstCapacity,
                                               compressedBytes, static_cast<size_t>(idx.size));
        if (ZSTD_isError(result)) {
            cerr << "zstd decompression error: " << ZSTD_getErrorName(result) << endl;
            return -1;
        }
        return static_cast<int32_t>(result);
    }
    z_stream infstream;
    infstream.zalloc = Z_NULL;
    infstream.zfree = Z_NULL;
    infstream.opaque = Z_NULL;
    infstream.avail_in = static_cast<uInt>(idx.size); // size of input
    infstream.next_in = (Bytef *) compressedBytes; // input char array
    infstream.avail_out = static_cast<uInt>(dstCapacity); // size of output
    infstream.next_out = (Bytef *) uncompressedBytes; // output char array
    // the actual decompression work.
    if (inflateInit(&infstream) != Z_OK) {
        cerr << "zlib initialization error" << endl;
        return -1;
    }
    const int status = inflate(&infstream, Z_FINISH);
    const uLong produced = infstream.total_out;
    inflateEnd(&infstream);
    // Z_FINISH reports Z_STREAM_END only when the whole block fit in the output
    // buffer. Anything else means the block needed more than DECOMPRESSION_FACTOR
    // times its stored size, which historically was silently truncated here.
    if (status != Z_STREAM_END) {
        cerr << "zlib decompression error: block did not fit in the output buffer" << endl;
        return -1;
    }
    return static_cast<int32_t>(produced);
}

long getNumRecordsInBlock(HiCFileStream &stream, indexEntry idx, int32_t version){
    if (idx.size <= 0) {
        return 0;
    }
    char *compressedBytes = stream.readCompressedBytes(idx);
    char *uncompressedBytes = new char[idx.size * DECOMPRESSION_FACTOR];
    int32_t uncompressedSize = decompressBlock(idx, compressedBytes, uncompressedBytes);
    if (uncompressedSize < 0) {
        free(compressedBytes);
        delete[] uncompressedBytes;
        throw runtime_error("Unable to decompress .hic block");
    }

    // create stream from buffer for ease of use
    memstream bufferin(uncompressedBytes, uncompressedSize);
    uint64_t nRecords;
    nRecords = static_cast<uint64_t>(readInt32FromFile(bufferin));
    free(compressedBytes);
    delete[] uncompressedBytes; // don't forget to delete your heap arrays in C++!
    return nRecords;
}

// this is the meat of reading the data.  takes in the block number and returns the set of contact records corresponding to
// that block.  the block data is compressed and must be decompressed using the zlib library functions
vector<contactRecord> readBlock(HiCFileStream &stream, indexEntry idx, int32_t version) {
    if (idx.size <= 0) {
        vector<contactRecord> v;
        return v;
    }
    char *compressedBytes = stream.readCompressedBytes(idx);
    char *uncompressedBytes = new char[idx.size * DECOMPRESSION_FACTOR];
    int32_t uncompressedSize = decompressBlock(idx, compressedBytes, uncompressedBytes);
    if (uncompressedSize < 0) {
        free(compressedBytes);
        delete[] uncompressedBytes;
        throw runtime_error("Unable to decompress .hic block");
    }

    // create stream from buffer for ease of use
    memstream bufferin(uncompressedBytes, uncompressedSize);
    uint64_t nRecords;
    nRecords = static_cast<uint64_t>(readInt32FromFile(bufferin));
    // A record cannot be encoded in fewer than 3 bytes, so a declared count that
    // exceeds the payload means the header is corrupt.
    if (nRecords > static_cast<uint64_t>(uncompressedSize)) {
        free(compressedBytes);
        delete[] uncompressedBytes;
        throw runtime_error("Corrupt .hic block: record count exceeds block size");
    }
    vector<contactRecord> v(nRecords);
    int32_t index = 0;
    // different versions have different specific formats
    if (version < 7) {
        for (uInt i = 0; i < nRecords; i++) {
            int32_t binX = readInt32FromFile(bufferin);
            int32_t binY = readInt32FromFile(bufferin);
            float counts = readFloatFromFile(bufferin);
            appendRecord(v, i, binX, binY, counts);
        }
    } else {
        int32_t binXOffset = readInt32FromFile(bufferin);
        int32_t binYOffset = readInt32FromFile(bufferin);
        bool useShort = readCharFromFile(bufferin) == 0; // yes this is opposite of usual

        bool useShortBinX = true;
        bool useShortBinY = true;
        if (version > 8) {
            useShortBinX = readCharFromFile(bufferin) == 0;
            useShortBinY = readCharFromFile(bufferin) == 0;
        }

        char type = readCharFromFile(bufferin);

        bool allCountsOne = false;
        bool useDeltaColumn = false;
        if (version > 9) {
            allCountsOne  = readCharFromFile(bufferin) != 0;
            useDeltaColumn = readCharFromFile(bufferin) != 0;
        }

        if (type == 1) {
            if (useShortBinX && useShortBinY) {
                int16_t rowCount = readInt16FromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    int32_t binY = binYOffset + readInt16FromFile(bufferin);
                    int16_t colCount = readInt16FromFile(bufferin);
                    int32_t prevColVal = 0;
                    for (int j = 0; j < colCount; j++) {
                        int32_t colVal = (int32_t)readInt16FromFile(bufferin);
                        if (useDeltaColumn) prevColVal += colVal; else prevColVal = colVal;
                        int32_t binX = binXOffset + prevColVal;
                        float counts = 1.0f;
                        if (!allCountsOne) {
                            counts = useShort ? (float)readInt16FromFile(bufferin) : readFloatFromFile(bufferin);
                        }
                        appendRecord(v, index++, binX, binY, counts);
                    }
                }
            } else if (useShortBinX && !useShortBinY) {
                int32_t rowCount = readInt32FromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    int32_t binY = binYOffset + readInt32FromFile(bufferin);
                    int16_t colCount = readInt16FromFile(bufferin);
                    int32_t prevColVal = 0;
                    for (int j = 0; j < colCount; j++) {
                        int32_t colVal = (int32_t)readInt16FromFile(bufferin);
                        if (useDeltaColumn) prevColVal += colVal; else prevColVal = colVal;
                        int32_t binX = binXOffset + prevColVal;
                        float counts = 1.0f;
                        if (!allCountsOne) {
                            counts = useShort ? (float)readInt16FromFile(bufferin) : readFloatFromFile(bufferin);
                        }
                        appendRecord(v, index++, binX, binY, counts);
                    }
                }
            } else if (!useShortBinX && useShortBinY) {
                int16_t rowCount = readInt16FromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    int32_t binY = binYOffset + readInt16FromFile(bufferin);
                    int32_t colCount = readInt32FromFile(bufferin);
                    int32_t prevColVal = 0;
                    for (int j = 0; j < colCount; j++) {
                        int32_t colVal = readInt32FromFile(bufferin);
                        if (useDeltaColumn) prevColVal += colVal; else prevColVal = colVal;
                        int32_t binX = binXOffset + prevColVal;
                        float counts = 1.0f;
                        if (!allCountsOne) {
                            counts = useShort ? (float)readInt16FromFile(bufferin) : readFloatFromFile(bufferin);
                        }
                        appendRecord(v, index++, binX, binY, counts);
                    }
                }
            } else {
                int32_t rowCount = readInt32FromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    int32_t binY = binYOffset + readInt32FromFile(bufferin);
                    int32_t colCount = readInt32FromFile(bufferin);
                    int32_t prevColVal = 0;
                    for (int j = 0; j < colCount; j++) {
                        int32_t colVal = readInt32FromFile(bufferin);
                        if (useDeltaColumn) prevColVal += colVal; else prevColVal = colVal;
                        int32_t binX = binXOffset + prevColVal;
                        float counts = 1.0f;
                        if (!allCountsOne) {
                            counts = useShort ? (float)readInt16FromFile(bufferin) : readFloatFromFile(bufferin);
                        }
                        appendRecord(v, index++, binX, binY, counts);
                    }
                }
            }
        } else if (type == 2) {
            int32_t nPts = readInt32FromFile(bufferin);
            int16_t w = readInt16FromFile(bufferin);

            for (int i = 0; i < nPts; i++) {
                //int32_t idx = (p.y - binOffset2) * w + (p.x - binOffset1);
                int32_t row = i / w;
                int32_t col = i - row * w;
                int32_t bin1 = binXOffset + col;
                int32_t bin2 = binYOffset + row;

                float counts;
                if (useShort) {
                    int16_t c = readInt16FromFile(bufferin);
                    if (c != -32768) {
                        appendRecord(v, index++, bin1, bin2, c);
                    }
                } else {
                    counts = readFloatFromFile(bufferin);
                    if (!isnan(counts)) {
                        appendRecord(v, index++, bin1, bin2, counts);
                    }
                }
            }
        }
    }
    free(compressedBytes);
    if (version >= 7) {
        v.resize(static_cast<size_t>(index));
    }
    delete[] uncompressedBytes; // don't forget to delete your heap arrays in C++!
    return v;
}

vector<contactRecord> readBlock(const string &fileName, indexEntry idx, int32_t version) {
    HiCFileStream stream(fileName);
    return readBlock(stream, idx, version);
}

// reads the normalization vector from the file at the specified location
vector<double> readNormalizationVector(istream &bufferin, int32_t version) {
    int64_t nValues;
    if (version > 8) {
        nValues = readInt64FromFile(bufferin);
    } else {
        nValues = (int64_t) readInt32FromFile(bufferin);
    }

    uint64_t numValues;
    numValues = static_cast<uint64_t>(nValues);
    vector<double> values(numValues);

    if (version > 8) {
        for (int i = 0; i < nValues; i++) {
            values[i] = (double) readFloatFromFile(bufferin);
        }
    } else {
        for (int i = 0; i < nValues; i++) {
            values[i] = readDoubleFromFile(bufferin);
        }
    }

    return values;
}

// Modify the BlockResult struct to include block number for sorting
struct BlockResult {
    vector<contactRecord> records;
    // int64_t: this used to be int32_t while holding a file offset, which
    // truncated (and could go negative) for any .hic larger than 2 GB.
    int64_t blockNumber;
};

// Add a comparison function for sorting BlockResults
bool compareBlockResults(const BlockResult &a, const BlockResult &b) {
    return a.blockNumber < b.blockNumber;
}

// Add this helper function that processes a single block
// Takes an already-open stream: constructing one per block meant a fresh
// ifstream, or over HTTP a fresh connection and TLS handshake, for every block.
BlockResult processBlock(HiCFileStream &stream, indexEntry idx, int32_t blockNumber, int32_t version,
                       int64_t *regionIndices, int32_t resolution,
                       const string &norm, vector<double> &c1Norm, vector<double> &c2Norm,
                       bool isIntra, const string &matrixType, vector<double> &expectedValues,
                       double avgCount) {
    BlockResult result;
    vector<contactRecord> records = readBlock(stream, idx, version);
    vector<contactRecord> filteredRecords;
    
    for (const contactRecord &rec : records) {
        // Widen before multiplying: binX and resolution are both int32_t, so the
        // product overflows for genomes past ~2.1 Gb before it ever reaches int64_t.
        int64_t x = static_cast<int64_t>(rec.binX) * resolution;
        int64_t y = static_cast<int64_t>(rec.binY) * resolution;

        if ((x >= regionIndices[0] && x <= regionIndices[1] &&
             y >= regionIndices[2] && y <= regionIndices[3]) ||
            (isIntra && y >= regionIndices[0] && y <= regionIndices[1] && 
             x >= regionIndices[2] && x <= regionIndices[3])) {

            float c = rec.counts;
            if (norm != "NONE") {
                // .at() so a short or corrupt normalization vector raises instead
                // of reading out of bounds, matching the V10 reader.
                c = static_cast<float>(c / (c1Norm.at(rec.binX) * c2Norm.at(rec.binY)));
            }
            if (matrixType == "oe") {
                if (isIntra) {
                    c = static_cast<float>(c / expectedValues[min(expectedValues.size() - 1,
                                                              (size_t) floor(abs(y - x) /
                                                                         resolution))]);
                } else {
                    c = static_cast<float>(c / avgCount);
                }
            } else if (matrixType == "expected") {
                if (isIntra) {
                    c = static_cast<float>(expectedValues[min(expectedValues.size() - 1,
                                                              (size_t) floor(abs(y - x) /
                                                                         resolution))]);
                } else {
                    c = static_cast<float>(avgCount);
                }
            }

            if (!isnan(c) && !isinf(c)) {
                if (x > INT32_MAX || y > INT32_MAX) {
                    throw runtime_error("Genomic coordinate exceeds the legacy contactRecord range");
                }
                contactRecord record = contactRecord();
                record.binX = static_cast<int32_t>(x);
                record.binY = static_cast<int32_t>(y);
                record.counts = c;
                filteredRecords.push_back(record);
            }
        }
    }
    
    result.records = filteredRecords;
    // Sort by the block number, which is what actually orders the output. The
    // file offset used here before was truncated into an int32_t sort key.
    result.blockNumber = blockNumber;
    return result;
}

void processBlockRecords(HiCFileStream &stream, indexEntry idx, int32_t version,
                         int64_t *regionIndices, int32_t resolution,
                         const string &norm, vector<double> &c1Norm, vector<double> &c2Norm,
                         bool isIntra, const string &matrixType, vector<double> &expectedValues,
                         double avgCount, const StrawRecordCallback &callback) {
    vector<contactRecord> records = readBlock(stream, idx, version);

    for (const contactRecord &rec : records) {
        // Widen before multiplying: binX and resolution are both int32_t, so the
        // product overflows for genomes past ~2.1 Gb before it ever reaches int64_t.
        int64_t x = static_cast<int64_t>(rec.binX) * resolution;
        int64_t y = static_cast<int64_t>(rec.binY) * resolution;

        if ((x >= regionIndices[0] && x <= regionIndices[1] &&
             y >= regionIndices[2] && y <= regionIndices[3]) ||
            (isIntra && y >= regionIndices[0] && y <= regionIndices[1] &&
             x >= regionIndices[2] && x <= regionIndices[3])) {

            float c = rec.counts;
            if (norm != "NONE") {
                // .at() so a short or corrupt normalization vector raises instead
                // of reading out of bounds, matching the V10 reader.
                c = static_cast<float>(c / (c1Norm.at(rec.binX) * c2Norm.at(rec.binY)));
            }
            if (matrixType == "oe") {
                if (isIntra) {
                    c = static_cast<float>(c / expectedValues[min(expectedValues.size() - 1,
                                                              (size_t) floor(abs(y - x) /
                                                                         resolution))]);
                } else {
                    c = static_cast<float>(c / avgCount);
                }
            } else if (matrixType == "expected") {
                if (isIntra) {
                    c = static_cast<float>(expectedValues[min(expectedValues.size() - 1,
                                                              (size_t) floor(abs(y - x) /
                                                                         resolution))]);
                } else {
                    c = static_cast<float>(avgCount);
                }
            }

            if (!isnan(c) && !isinf(c)) {
                if (x > INT32_MAX || y > INT32_MAX) {
                    throw runtime_error("Genomic coordinate exceeds the legacy contactRecord range");
                }
                contactRecord record = contactRecord();
                record.binX = static_cast<int32_t>(x);
                record.binY = static_cast<int32_t>(y);
                record.counts = c;
                callback(record);
            }
        }
    }
}

class ThreadPool {
public:
    explicit ThreadPool(size_t numThreads) : stop(false) {
        for(size_t i = 0; i < numThreads; ++i) {
            workers.emplace_back([this] {
                while(true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(queueMutex);
                        condition.wait(lock, [this] { 
                            return stop || !tasks.empty(); 
                        });
                        if(stop && tasks.empty()) {
                            return;
                        }
                        task = std::move(tasks.front());
                        tasks.pop();
                    }
                    task();
                }
            });
        }
    }

    template<class F>
    auto enqueue(F&& f) -> std::future<typename std::result_of<F()>::type> {
        using return_type = typename std::result_of<F()>::type;
        auto task = std::make_shared<std::packaged_task<return_type()>>(std::forward<F>(f));
        std::future<return_type> res = task->get_future();
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            if(stop) {
                throw std::runtime_error("enqueue on stopped ThreadPool");
            }
            tasks.emplace([task]() { (*task)(); });
        }
        condition.notify_one();
        return res;
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            stop = true;
        }
        condition.notify_all();
        for(std::thread &worker: workers) {
            worker.join();
        }
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queueMutex;
    std::condition_variable condition;
    bool stop;
};

class MatrixZoomData {
public:
    bool isIntra;
    string fileName;
    int64_t myFilePos = 0LL;
    vector<double> expectedValues;
    bool foundFooter = false;
    vector<double> c1Norm;
    vector<double> c2Norm;
    int32_t c1 = 0;
    int32_t c2 = 0;
    string matrixType;
    string norm;
    int32_t version = 0;
    int32_t resolution = 0;
    int32_t numBins1 = 0;
    int32_t numBins2 = 0;
    float sumCounts = 0.0f;
    // Only set when the requested resolution is present in the matrix record;
    // getBlockNumbers divides by blockBinCount, so they must never be read raw.
    int32_t blockBinCount = 0, blockColumnCount = 0;
    map<int32_t, indexEntry> blockMap;
    double avgCount = 0.0;

    MatrixZoomData(const chromosome &chrom1, const chromosome &chrom2, const string &matrixType,
                   const string &norm, const string &unit, int32_t resolution,
                   int32_t &version, int64_t &master, int64_t &totalFileSize,
                   const string &fileName) {
        this->version = version;
        this->fileName = fileName;
        int32_t c01 = chrom1.index;
        int32_t c02 = chrom2.index;
        if (c01 <= c02) { // default is ok
            this->c1 = c01;
            this->c2 = c02;
            this->numBins1 = static_cast<int32_t>(chrom1.length / resolution);
            this->numBins2 = static_cast<int32_t>(chrom2.length / resolution);
        } else { // flip
            this->c1 = c02;
            this->c2 = c01;
            this->numBins1 = static_cast<int32_t>(chrom2.length / resolution);
            this->numBins2 = static_cast<int32_t>(chrom1.length / resolution);
        }
        isIntra = c1 == c2;

        this->matrixType = matrixType;
        this->norm = norm;
        this->resolution = resolution;

        HiCFileStream stream(fileName);
        indexEntry c1NormEntry{}, c2NormEntry{};

        if (stream.isHttp) {
            foundFooter = readFooterURL(stream.curl, master, version, c1, c2, matrixType, norm, unit,
                                     resolution,
                                     myFilePos,
                                     c1NormEntry, c2NormEntry, expectedValues);
        } else {
            stream.fin.seekg(master, ios::beg);
            foundFooter = readFooter(stream.fin, master, version, c1, c2, matrixType, norm,
                                     unit,
                                     resolution, myFilePos,
                                     c1NormEntry, c2NormEntry, expectedValues);
        }

        if (!foundFooter) {
            return;
        }

        if (norm != "NONE") {
            if (c1NormEntry.size <= 0 || (!isIntra && c2NormEntry.size <= 0)) {
                cerr << "Normalization " << norm << " is not available for the requested chromosomes, unit, and resolution" << endl;
                foundFooter = false;
                return;
            }
            c1Norm = readNormalizationVectorFromFooter(c1NormEntry, version, fileName);
            if (isIntra) {
                c2Norm = c1Norm;
            } else {
                c2Norm = readNormalizationVectorFromFooter(c2NormEntry, version, fileName);
            }
        }

        // Reuse the footer stream rather than opening a second one; over HTTP
        // that second open cost an extra connection for every matrix.
        bool foundResolution = false;
        if (stream.isHttp) {
            // readMatrix will assign blockBinCount and blockColumnCount
            blockMap = readMatrixHttp(stream.curl, myFilePos, unit, resolution, sumCounts,
                                      blockBinCount,
                                      blockColumnCount, foundResolution);
        } else {
            // readMatrix will assign blockBinCount and blockColumnCount
            blockMap = readMatrix(stream.fin, myFilePos, unit, resolution, sumCounts,
                                  blockBinCount,
                                  blockColumnCount, foundResolution);
        }
        stream.close();

        // Without the matrix record there is no block geometry, so every later
        // query would divide by an unset blockBinCount.
        if (!foundResolution || blockBinCount <= 0 || blockColumnCount <= 0) {
            foundFooter = false;
            return;
        }

        if (!isIntra && numBins1 > 0 && numBins2 > 0) {
            avgCount = (sumCounts / numBins1) / numBins2;   // <= trying to avoid overflows
        }
    }

    static vector<double> readNormalizationVectorFromFooter(indexEntry cNormEntry, int32_t &version,
                                                            const string &fileName) {
        char *buffer = readCompressedBytesFromFile(fileName, cNormEntry);
        memstream bufferin(buffer, cNormEntry.size);
        vector<double> cNorm = readNormalizationVector(bufferin, version);
        free(buffer);
        return cNorm;
    }

    static bool isInRange(int32_t r, int32_t c, int32_t numRows, int32_t numCols) {
        return 0 <= r && r < numRows && 0 <= c && c < numCols;
    }

    set<int32_t> getBlockNumbers(int64_t *regionIndices) const {
        if (version > 8 && isIntra) {
            return getBlockNumbersForRegionFromBinPositionV9Intra(regionIndices, blockBinCount, blockColumnCount);
        } else {
            return getBlockNumbersForRegionFromBinPosition(regionIndices, blockBinCount, blockColumnCount, isIntra);
        }
    }

    vector<double> getNormVector(int32_t index) {
        if (index == c1) {
            return c1Norm;
        } else if (index == c2) {
            return c2Norm;
        }
        cerr << "Invalid index provided: " << index << endl;
        cerr << "Should be either " << c1 << " or " << c2 << endl;
        vector<double> v;
        return v;
    }

    vector<double> getExpectedValues() {
        return expectedValues;
    }

    vector<contactRecord> getRecords(int64_t gx0, int64_t gx1, int64_t gy0, int64_t gy1) {
        if (!foundFooter) {
            vector<contactRecord> v;
            return v;
        }
        int64_t origRegionIndices[] = {gx0, gx1, gy0, gy1};
        int64_t regionIndices[4];
        convertGenomeToBinPos(origRegionIndices, regionIndices, resolution);

        set<int32_t> blockNumbers = getBlockNumbers(regionIndices);

        // Pre-filter to blocks that actually exist in the map.
        // blockMap[key] inserts a default entry when key is missing, which
        // is a data race when multiple threads call it concurrently.
        // Filtering first keeps all subsequent lambda accesses read-only.
        vector<int32_t> existingBlocks;
        for (int32_t blockNumber : blockNumbers) {
            if (blockMap.count(blockNumber) > 0) {
                existingBlocks.push_back(blockNumber);
            }
        }

        vector<BlockResult> allResults;
        if (existingBlocks.empty()) {
            return vector<contactRecord>();
        }

        // Adjust thread count based on block count and available cores
        // hardware_concurrency() is allowed to return 0; subtracting first would
        // underflow to UINT_MAX and spawn one thread per block.
        const unsigned int cores = max(1u, thread::hardware_concurrency());
        unsigned int maxThreads = max(1u, cores - 1);
        unsigned int numThreads = max(1u, min(
            maxThreads,                // Don't use more than available cores minus one
            static_cast<unsigned int>(existingBlocks.size())  // Don't create more threads than blocks
        ));

        // One task per worker over a contiguous slice of the block list, so each
        // worker opens a single stream instead of one per block. A stream cannot
        // be shared across threads, hence one per task rather than one per query.
        vector<future<vector<BlockResult>>> futures;
        {
            ThreadPool pool(numThreads);
            const size_t total = existingBlocks.size();
            for (unsigned int t = 0; t < numThreads; ++t) {
                const size_t begin = total * t / numThreads;
                const size_t end = total * (t + 1) / numThreads;
                if (begin >= end) continue;
                futures.push_back(
                    pool.enqueue([this, begin, end, &existingBlocks, &origRegionIndices]() {
                        vector<BlockResult> partial;
                        partial.reserve(end - begin);
                        HiCFileStream stream(fileName);
                        for (size_t i = begin; i < end; ++i) {
                            const int32_t blockNumber = existingBlocks[i];
                            partial.push_back(processBlock(
                                stream, blockMap.at(blockNumber), blockNumber, version,
                                origRegionIndices, resolution,
                                norm, c1Norm, c2Norm, isIntra,
                                matrixType, expectedValues, avgCount
                            ));
                        }
                        return partial;
                    })
                );
            }

            // Collect all results
            allResults.reserve(existingBlocks.size());
            for (auto& future : futures) {
                vector<BlockResult> partial = future.get();
                for (auto &entry : partial) {
                    allResults.push_back(std::move(entry));
                }
            }
        }

        // Sort results by block number to maintain consistent order
        sort(allResults.begin(), allResults.end(), compareBlockResults);

        // Combine all records in sorted order
        vector<contactRecord> records;
        size_t totalSize = 0;
        for (const auto& result : allResults) {
            totalSize += result.records.size();
        }
        records.reserve(totalSize);
        
        for (const auto& result : allResults) {
            records.insert(records.end(), 
                          result.records.begin(), 
                          result.records.end());
        }

        return records;
    }

    void streamRecords(int64_t gx0, int64_t gx1, int64_t gy0, int64_t gy1,
                       const StrawRecordCallback &callback) {
        if (!foundFooter) {
            return;
        }
        int64_t origRegionIndices[] = {gx0, gx1, gy0, gy1};
        int64_t regionIndices[4];
        convertGenomeToBinPos(origRegionIndices, regionIndices, resolution);

        set<int32_t> blockNumbers = getBlockNumbers(regionIndices);
        // One stream for the whole query rather than one per block.
        HiCFileStream stream(fileName);
        for (int32_t blockNumber : blockNumbers) {
            auto found = blockMap.find(blockNumber);
            if (found == blockMap.end()) {
                continue;
            }
            processBlockRecords(
                stream, found->second, version,
                origRegionIndices, resolution,
                norm, c1Norm, c2Norm, isIntra,
                matrixType, expectedValues, avgCount,
                callback
            );
        }
    }

    vector<vector<float> > getRecordsAsMatrix(int64_t gx0, int64_t gx1, int64_t gy0, int64_t gy1) {
        vector<contactRecord> records = this->getRecords(gx0, gx1, gy0, gy1);
        if (records.empty()) {
            vector<vector<float> > res = vector<vector<float> >(1, vector<float>(1, 0));
            return res;
        }

        int64_t origRegionIndices[] = {gx0, gx1, gy0, gy1};
        int64_t regionIndices[4];
        convertGenomeToBinPos(origRegionIndices, regionIndices, resolution);

        int64_t originR = regionIndices[0];
        int64_t endR = regionIndices[1];
        int64_t originC = regionIndices[2];
        int64_t endC = regionIndices[3];
        int32_t numRows = endR - originR + 1;
        int32_t numCols = endC - originC + 1;
        vector<vector<float> > matrix;
        for (int32_t i = 0; i < numRows; i++) {
            matrix.emplace_back(numCols, 0);
        }

        for (contactRecord cr : records) {
            if (isnan(cr.counts) || isinf(cr.counts)) continue;
            int32_t r = cr.binX / resolution - originR;
            int32_t c = cr.binY / resolution - originC;
            if (isInRange(r, c, numRows, numCols)) {
                matrix[r][c] = cr.counts;
            }
            if (isIntra) {
                r = cr.binY / resolution - originR;
                c = cr.binX / resolution - originC;
                if (isInRange(r, c, numRows, numCols)) {
                    matrix[r][c] = cr.counts;
                }
            }
        }
        return matrix;
    }

    int64_t getNumberOfTotalRecords() {
        if (!foundFooter) {
            return 0;
        }
        int64_t regionIndices[4] = {0, numBins1, 0, numBins2};
        set<int32_t> blockNumbers = getBlockNumbers(regionIndices);
        int64_t total = 0;
        HiCFileStream stream(fileName);
        for (int32_t blockNumber : blockNumbers) {
            // find, not operator[]: the latter inserts an empty entry for every
            // block number the region covers but the file does not contain.
            const auto found = blockMap.find(blockNumber);
            if (found == blockMap.end()) continue;
            total += getNumRecordsInBlock(stream, found->second, version);
        }
        return total;
    }
};

class HiCFile {
public:
    string prefix = "http"; // HTTP code
    int64_t master = 0LL;
    map<string, chromosome> chromosomeMap;
    string genomeID;
    int32_t numChromosomes = 0;
    int32_t version = 0;
    int64_t nviPosition = 0LL;
    int64_t nviLength = 0LL;
    vector<int32_t> resolutions;
    vector<int32_t> fragResolutions;
    vector<pair<string, string>> attributes;
    static thread_local int64_t totalFileSize;
    string fileName;

    static size_t hdf(char *b, size_t size, size_t nitems, void *userdata) {
        (void)userdata;
        size_t numbytes = size * nitems;
        string s(b, numbytes);
        int32_t found;
        found = static_cast<int32_t>(s.find("content-range"));
        if ((size_t) found == string::npos) {
            found = static_cast<int32_t>(s.find("Content-Range"));
        }
        if ((size_t) found != string::npos) {
            int32_t found2;
            found2 = static_cast<int32_t>(s.find('/'));
            //content-range: bytes 0-100000/891471462
            if ((size_t) found2 != string::npos) {
                string total = s.substr(found2 + 1);
                totalFileSize = stol(total);
            }
        }

        return numbytes;
    }

    static CURL *oneTimeInitCURL(const char *url) {
        CURL *curl = initCURL(url);
        curl_easy_setopt(curl, CURLOPT_HEADERFUNCTION, hdf);
        return curl;
    }

    explicit HiCFile(const string &fileName) {
        this->fileName = fileName;

        // read header into buffer; 100K should be sufficient
        if (std::strncmp(fileName.c_str(), prefix.c_str(), prefix.size()) == 0) {
            CURL *curl;
            curl = oneTimeInitCURL(fileName.c_str());
            int64_t received = 0;
            char *buffer = getData(curl, 0, 100000, received);
            memstream bufin(buffer, received);
            chromosomeMap = readHeader(bufin, master, genomeID, numChromosomes,
                                       version, nviPosition, nviLength, &attributes);
            if (master < 0) {
                curl_easy_cleanup(curl);
                free(buffer);
                throw runtime_error("invalid or unsupported .hic header");
            }
            readResolutionsFromHeader(bufin, resolutions, fragResolutions);
            curl_easy_cleanup(curl);
            free(buffer);
        } else {
            ifstream fin;
            fin.open(fileName, fstream::in | fstream::binary);
            if (!fin) {
                throw runtime_error("File " + fileName + " cannot be opened for reading");
            }
            chromosomeMap = readHeader(fin, master, genomeID, numChromosomes,
                                       version, nviPosition, nviLength, &attributes);
            if (master < 0) throw runtime_error("invalid or unsupported .hic header");
            readResolutionsFromHeader(fin, resolutions, fragResolutions);
            fin.close();
        }
    }

    string getGenomeID() const {
        return genomeID;
    }

    vector<int32_t> getResolutions() const {
        return resolutions;
    }

    vector<chromosome> getChromosomes() {
        // chromosomeMap is keyed by name, so duplicate names in the header make it
        // smaller than numChromosomes while indices still run to numChromosomes-1.
        // Size by the index range and drop anything outside it.
        size_t count = chromosomeMap.size();
        for (const auto &entry : chromosomeMap) {
            if (entry.second.index >= 0 && static_cast<size_t>(entry.second.index) + 1 > count) {
                count = static_cast<size_t>(entry.second.index) + 1;
            }
        }
        vector<chromosome> chromosomes(count);
        for (const auto &entry : chromosomeMap) {
            const chromosome &chrom = entry.second;
            if (chrom.index < 0 || static_cast<size_t>(chrom.index) >= chromosomes.size()) {
                cerr << "Ignoring chromosome " << chrom.name << " with out-of-range index "
                     << chrom.index << endl;
                continue;
            }
            chromosomes[chrom.index] = chrom;
        }

        return chromosomes;
    }

    MatrixZoomData * getMatrixZoomData(const string &chr1, const string &chr2, const string &matrixType,
                                       const string &norm, const string &unit, int32_t resolution) {
        // operator[] would insert a default chromosome (index 0, length 0) for an
        // unknown name and silently return an empty result instead of erroring.
        const auto first = chromosomeMap.find(chr1);
        const auto second = chromosomeMap.find(chr2);
        if (first == chromosomeMap.end()) {
            throw runtime_error("chromosome " + chr1 + " not found in the file");
        }
        if (second == chromosomeMap.end()) {
            throw runtime_error("chromosome " + chr2 + " not found in the file");
        }
        return new MatrixZoomData(first->second, second->second, (matrixType), (norm), (unit),
                                  resolution, version, master, totalFileSize, fileName);
    }
};

thread_local int64_t HiCFile::totalFileSize = 0LL;

// The chromosome map is passed by const reference: taking it by value copied the
// whole map on every call.
void parsePositions(const string &chrLoc, string &chrom, int64_t &pos1, int64_t &pos2,
                    const std::map<string, chromosome> &map) {
    string x, y;
    stringstream ss(chrLoc);
    getline(ss, chrom, ':');
    if (map.count(chrom) == 0) {
        throw runtime_error("chromosome " + chrom + " not found in the file");
    }

    if (getline(ss, x, ':') && getline(ss, y, ':')) {
        pos1 = stol(x);
        pos2 = stol(y);
    } else {
        pos1 = 0LL;
        pos2 = map.at(chrom).length;
    }
}

struct StrawPreparedQuery::Impl {
    string fileName;
    string matrixType;
    string normalization;
    string firstChromosome;
    string secondChromosome;
    string unit;
    int32_t resolution;
    bool transpose = false;
    bool intra = false;
    unique_ptr<HiCFile> legacyFile;
    unique_ptr<MatrixZoomData> legacyZoom;
    unique_ptr<straw_v10::File> v10File;

    Impl(const string &path, const string &type, const string &norm,
         const string &first, const string &second, const string &queryUnit,
         int32_t queryResolution)
        : fileName(path), matrixType(type), normalization(norm),
          firstChromosome(first), secondChromosome(second), unit(queryUnit),
          resolution(queryResolution) {
        if (path.empty() || first.empty() || second.empty() || queryResolution <= 0)
            throw StrawException(StrawErrorCode::InvalidArgument, "invalid prepared query arguments");
        if (queryUnit != "BP" && queryUnit != "FRAG")
            throw StrawException(StrawErrorCode::InvalidArgument, "unit must be BP or FRAG");
        if (type != "observed" && type != "oe" && type != "expected")
            throw StrawException(StrawErrorCode::InvalidArgument,
                                 "matrix type must be observed, oe, or expected");

        const bool remote = path.compare(0, 4, "http") == 0;
        if (!remote) {
            ifstream input(path, ios::binary);
            if (!input)
                throw StrawException(StrawErrorCode::Io, "file cannot be opened for reading");
        }

        try {
            if (straw_v10::isV10(path)) {
                v10File.reset(new straw_v10::File(path));
                const vector<chromosome> chromosomes = v10File->chromosomes();
                auto hasChromosome = [&](const string &name) {
                    return any_of(chromosomes.begin(), chromosomes.end(),
                                  [&](const chromosome &value) { return value.name == name; });
                };
                if (!hasChromosome(first) || !hasChromosome(second))
                    throw StrawException(StrawErrorCode::NotFound,
                                         "requested chromosome is not present in the file");
                const vector<int32_t> resolutions = v10File->resolutions(queryUnit);
                if (find(resolutions.begin(), resolutions.end(), queryResolution) == resolutions.end())
                    throw StrawException(StrawErrorCode::Unavailable,
                                         "requested resolution is not available");
                vector<string> norms = v10File->normalizations();
                if (norm != "NONE" && find(norms.begin(), norms.end(), norm) == norms.end())
                    throw StrawException(StrawErrorCode::Unavailable,
                                         "requested normalization is not available");
                return;
            }

            legacyFile.reset(new HiCFile(path));
            if (legacyFile->version < 6 || legacyFile->version > 9)
                throw StrawException(StrawErrorCode::UnsupportedVersion,
                                     "unsupported .hic version");
            auto firstIt = legacyFile->chromosomeMap.find(first);
            auto secondIt = legacyFile->chromosomeMap.find(second);
            if (firstIt == legacyFile->chromosomeMap.end() ||
                secondIt == legacyFile->chromosomeMap.end())
                throw StrawException(StrawErrorCode::NotFound,
                                     "requested chromosome is not present in the file");
            const vector<int32_t> &resolutions = queryUnit == "BP"
                ? legacyFile->resolutions : legacyFile->fragResolutions;
            if (find(resolutions.begin(), resolutions.end(), queryResolution) == resolutions.end())
                throw StrawException(StrawErrorCode::Unavailable,
                                     "requested resolution is not available");
            transpose = firstIt->second.index > secondIt->second.index;
            intra = first == second;
            const string &storedFirst = transpose ? second : first;
            const string &storedSecond = transpose ? first : second;
            legacyZoom.reset(legacyFile->getMatrixZoomData(
                storedFirst, storedSecond, type, norm, queryUnit, queryResolution));
            if (!legacyZoom || !legacyZoom->foundFooter)
                throw StrawException(StrawErrorCode::Unavailable,
                                     "requested matrix or normalization is not available");
        } catch (const StrawException &) {
            throw;
        } catch (const exception &error) {
            throw StrawException(remote ? StrawErrorCode::Io : StrawErrorCode::CorruptFile,
                                 error.what());
        }
    }

    void stream(int64_t xStart, int64_t xEnd, int64_t yStart, int64_t yEnd,
                const StrawRecordCallback &callback) {
        try {
            if (v10File) {
                v10File->stream(matrixType, normalization,
                    firstChromosome + ":" + to_string(xStart) + ":" + to_string(xEnd),
                    secondChromosome + ":" + to_string(yStart) + ":" + to_string(yEnd),
                    unit, resolution, callback);
                return;
            }
            auto orientedCallback = [&](const contactRecord &input) {
                contactRecord record = input;
                if (transpose) {
                    swap(record.binX, record.binY);
                } else if (intra) {
                    const bool direct = record.binX >= xStart && record.binX <= xEnd &&
                                        record.binY >= yStart && record.binY <= yEnd;
                    if (!direct) swap(record.binX, record.binY);
                }
                callback(record);
            };
            if (transpose)
                legacyZoom->streamRecords(yStart, yEnd, xStart, xEnd, orientedCallback);
            else
                legacyZoom->streamRecords(xStart, xEnd, yStart, yEnd, orientedCallback);
        } catch (const StrawException &) {
            throw;
        } catch (const exception &error) {
            throw StrawException(StrawErrorCode::CorruptFile, error.what());
        }
    }
};

StrawPreparedQuery::StrawPreparedQuery(const string &fileName, const string &matrixType,
                                       const string &normalization, const string &firstChromosome,
                                       const string &secondChromosome, const string &unit,
                                       int32_t resolution)
    : impl(new Impl(fileName, matrixType, normalization, firstChromosome,
                    secondChromosome, unit, resolution)) {}

StrawPreparedQuery::~StrawPreparedQuery() = default;

void StrawPreparedQuery::streamWindow(int64_t xStart, int64_t xEnd,
                                      int64_t yStart, int64_t yEnd,
                                      const StrawRecordCallback &callback) {
    if (xStart < 0 || yStart < 0 || xEnd < xStart || yEnd < yStart)
        throw StrawException(StrawErrorCode::InvalidArgument, "invalid query window");
    impl->stream(xStart, xEnd, yStart, yEnd, callback);
}

bool strawStream(const string &matrixType, const string &norm, const string &fileName, const string &chr1loc,
                 const string &chr2loc, const string &unit, int32_t binsize,
                 const StrawRecordCallback &callback) {
    if (straw_v10::isV10(fileName)) {
        straw_v10::File(fileName).stream(matrixType, norm, chr1loc, chr2loc, unit, binsize, callback);
        return true;
    }

    if (!(unit == "BP" || unit == "FRAG")) {
        cerr << "Norm specified incorrectly, must be one of <BP/FRAG>" << endl;
        cerr << "Usage: straw [observed/oe/expected] <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>"
             << endl;
        return false;
    }

    // Owned: parsePositions throws for an unknown chromosome, which previously
    // leaked the HiCFile (and any MatrixZoomData) on the way out.
    unique_ptr<HiCFile> hiCFile(new HiCFile(fileName));
    string chr1, chr2;
    int64_t origRegionIndices[4] = {-100LL, -100LL, -100LL, -100LL};
    parsePositions((chr1loc), chr1, origRegionIndices[0], origRegionIndices[1], hiCFile->chromosomeMap);
    parsePositions((chr2loc), chr2, origRegionIndices[2], origRegionIndices[3], hiCFile->chromosomeMap);

    // parsePositions already verified both names are present.
    bool transpose = hiCFile->chromosomeMap.at(chr1).index > hiCFile->chromosomeMap.at(chr2).index;
    bool intra = chr1 == chr2;
    auto orientedCallback = [&](const contactRecord &input) {
        contactRecord record = input;
        if (transpose) {
            // Matrix storage is ordered by chromosome index, but the public API
            // returns coordinates in the chromosome order requested by the
            // caller, as the V10 reader does.
            std::swap(record.binX, record.binY);
        } else if (intra) {
            bool direct = record.binX >= origRegionIndices[0] &&
                          record.binX <= origRegionIndices[1] &&
                          record.binY >= origRegionIndices[2] &&
                          record.binY <= origRegionIndices[3];
            if (!direct) {
                // Cis contacts are stored above the diagonal. Reflect contacts
                // selected through the symmetric half of the query so a below-
                // diagonal request is returned in the requested orientation.
                std::swap(record.binX, record.binY);
            }
        }
        callback(record);
    };

    if (transpose) {
        unique_ptr<MatrixZoomData> mzd(hiCFile->getMatrixZoomData(chr2, chr1, matrixType, norm, unit, binsize));
        mzd->streamRecords(origRegionIndices[2], origRegionIndices[3], origRegionIndices[0], origRegionIndices[1],
                           orientedCallback);
    } else {
        unique_ptr<MatrixZoomData> mzd(hiCFile->getMatrixZoomData(chr1, chr2, matrixType, norm, unit, binsize));
        mzd->streamRecords(origRegionIndices[0], origRegionIndices[1], origRegionIndices[2], origRegionIndices[3],
                           orientedCallback);
    }
    return true;
}

vector<contactRecord> straw(const string &matrixType, const string &norm, const string &fileName, const string &chr1loc,
                            const string &chr2loc, const string &unit, int32_t binsize) {
    vector<contactRecord> records;
    strawStream(matrixType, norm, fileName, chr1loc, chr2loc, unit, binsize,
                [&records](const contactRecord &record) {
                    records.push_back(record);
                });
    return records;
}

vector<vector<float> > strawAsMatrix(const string &matrixType, const string &norm, const string &fileName, const string &chr1loc,
                   const string &chr2loc, const string &unit, int32_t binsize) {
    if (straw_v10::isV10(fileName))
        return straw_v10::File(fileName).matrix(matrixType, norm, chr1loc, chr2loc, unit, binsize);

    if (!(unit == "BP" || unit == "FRAG")) {
        cerr << "Norm specified incorrectly, must be one of <BP/FRAG>" << endl;
        cerr << "Usage: straw [observed/oe/expected] <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>"
             << endl;
        vector<vector<float> > res = vector<vector<float> >(1, vector<float>(1, 0));
        return res;
    }

    unique_ptr<HiCFile> hiCFile(new HiCFile(fileName));
    string chr1, chr2;
    int64_t origRegionIndices[4] = {-100LL, -100LL, -100LL, -100LL};
    parsePositions((chr1loc), chr1, origRegionIndices[0], origRegionIndices[1], hiCFile->chromosomeMap);
    parsePositions((chr2loc), chr2, origRegionIndices[2], origRegionIndices[3], hiCFile->chromosomeMap);

    vector<vector<float> > result;
    if (hiCFile->chromosomeMap.at(chr1).index > hiCFile->chromosomeMap.at(chr2).index) {
        unique_ptr<MatrixZoomData> mzd(hiCFile->getMatrixZoomData(chr2, chr1, matrixType, norm, unit, binsize));
        vector<vector<float> > stored = mzd->getRecordsAsMatrix(
            origRegionIndices[2], origRegionIndices[3], origRegionIndices[0], origRegionIndices[1]);
        if (stored.size() == 1 && stored[0].size() == 1) {
            result = std::move(stored);
        } else {
            size_t rows = stored.size(), cols = stored[0].size();
            result.assign(cols, vector<float>(rows, 0));
            for (size_t r = 0; r < rows; ++r)
                for (size_t c = 0; c < cols; ++c)
                    result[c][r] = stored[r][c];
        }
    } else {
        unique_ptr<MatrixZoomData> mzd(hiCFile->getMatrixZoomData(chr1, chr2, matrixType, norm, unit, binsize));
        result = mzd->getRecordsAsMatrix(
            origRegionIndices[0], origRegionIndices[1], origRegionIndices[2], origRegionIndices[3]);
    }
    return result;
}

vector<chromosome> getChromosomesForFile(const string &fileName) {
    if (straw_v10::isV10(fileName)) return straw_v10::File(fileName).chromosomes();

    HiCFile hiCFile(fileName);
    return hiCFile.getChromosomes();
}

vector<int32_t> getResolutionsForFile(const string &fileName, const string &unit) {
    if (unit != "BP" && unit != "FRAG") throw invalid_argument("unit must be BP or FRAG");
    if (straw_v10::isV10(fileName)) return straw_v10::File(fileName).resolutions(unit);

    HiCFile hiCFile(fileName);
    return unit == "BP" ? hiCFile.resolutions : hiCFile.fragResolutions;
}

string getGenomeForFile(const string &fileName) {
    if (straw_v10::isV10(fileName)) return straw_v10::File(fileName).genome();
    return HiCFile(fileName).genomeID;
}

int32_t getVersionForFile(const string &fileName) {
    if (straw_v10::isV10(fileName)) return 10;
    return HiCFile(fileName).version;
}

vector<string> getNormalizationsForFile(const string &fileName) {
    if (straw_v10::isV10(fileName)) {
        vector<string> values = straw_v10::File(fileName).normalizations();
        if (find(values.begin(), values.end(), "NONE") == values.end()) values.push_back("NONE");
        return values;
    }
    // Legacy normalization indexes are resolution-specific and are queried by
    // name. NONE is the only normalization guaranteed to exist in every file.
    return {"NONE"};
}

vector<pair<string, string>> getAttributesForFile(const string &fileName) {
    if (straw_v10::isV10(fileName)) return straw_v10::File(fileName).attributes();
    return HiCFile(fileName).attributes;
}

StrawFileInfo getFileInfo(const string &fileName) {
    StrawFileInfo info;
    if (straw_v10::isV10(fileName)) {
        straw_v10::File file(fileName);
        info.version = 10;
        info.genome = file.genome();
        info.chromosomes = file.chromosomes();
        info.bpResolutions = file.resolutions("BP");
        info.fragResolutions = file.resolutions("FRAG");
        info.normalizations = file.normalizations();
        if (find(info.normalizations.begin(), info.normalizations.end(), "NONE") ==
            info.normalizations.end()) {
            info.normalizations.push_back("NONE");
        }
        info.attributes = file.attributes();
        return info;
    }
    HiCFile file(fileName);
    info.version = file.version;
    info.genome = file.genomeID;
    info.chromosomes = file.getChromosomes();
    info.bpResolutions = file.resolutions;
    info.fragResolutions = file.fragResolutions;
    // Legacy normalization indexes are resolution-specific and are queried by
    // name; NONE is the only one guaranteed to exist in every file.
    info.normalizations = {"NONE"};
    info.attributes = file.attributes;
    return info;
}

void forEachRawObservedBlock(const string &fileName,
                             const string &chr1,
                             const string &chr2,
                             int32_t binsize,
                             const StrawBlockCallback &processor) {
    if (straw_v10::isV10(fileName)) {
        straw_v10::File file(fileName);
        vector<contactRecord> batch;
        file.raw(chr1, chr2, "BP", binsize, 0, UINT32_MAX, 0, UINT32_MAX,
                 [&](const straw_v10::Record& r) {
                     if (r.binX > INT32_MAX || r.binY > INT32_MAX)
                         throw runtime_error("V10: bin coordinate exceeds legacy API; use File::raw");
                     batch.push_back({static_cast<int32_t>(r.binX), static_cast<int32_t>(r.binY),
                                      r.isScore ? r.score : static_cast<float>(r.count)});
                     if (batch.size() == 8192) { processor(batch); batch.clear(); }
                 });
        if (!batch.empty()) processor(batch);
        return;
    }

    HiCFile hiCFile(fileName);
    string first = chr1;
    string second = chr2;
    if (!hiCFile.chromosomeMap.count(first) || !hiCFile.chromosomeMap.count(second)) {
        throw runtime_error("chromosome not found in the file");
    }
    const bool transpose = hiCFile.chromosomeMap.at(first).index > hiCFile.chromosomeMap.at(second).index;
    if (transpose) {
        swap(first, second);
    }

    MatrixZoomData *mzd = hiCFile.getMatrixZoomData(first, second, "observed", "NONE", "BP", binsize);
    if (mzd == nullptr || !mzd->foundFooter) {
        delete mzd;
        return;
    }

    int64_t regionIndices[4] = {0, mzd->numBins1, 0, mzd->numBins2};
    set<int32_t> blockNumbers = mzd->getBlockNumbers(regionIndices);
    HiCFileStream stream(fileName);
    for (int32_t blockNumber : blockNumbers) {
        const auto found = mzd->blockMap.find(blockNumber);
        if (found == mzd->blockMap.end()) {
            continue;
        }

        vector<contactRecord> blockRecords = readBlock(stream, found->second, mzd->version);
        if (!blockRecords.empty()) {
            if (transpose) {
                for (auto &record : blockRecords) {
                    std::swap(record.binX, record.binY);
                }
            }
            processor(blockRecords);
        }
    }

    delete mzd;
}

struct StrawRawReader::Impl {
    HiCFile file;
    HiCFileStream stream;
    explicit Impl(const string &fileName) : file(fileName), stream(fileName) {}
};

StrawRawReader::StrawRawReader(const string &fileName) : impl(new Impl(fileName)) {}
StrawRawReader::~StrawRawReader() = default;

bool StrawRawReader::forEachBlock(const string &chr1, const string &chr2, int32_t binsize,
                                  const StrawBlockCallback &processor) {
    HiCFile &hiCFile = impl->file;
    if (!hiCFile.chromosomeMap.count(chr1) || !hiCFile.chromosomeMap.count(chr2)) {
        return false;
    }
    string first = chr1;
    string second = chr2;
    const bool transpose = hiCFile.chromosomeMap.at(first).index > hiCFile.chromosomeMap.at(second).index;
    if (transpose) {
        swap(first, second);
    }

    unique_ptr<MatrixZoomData> mzd(
        hiCFile.getMatrixZoomData(first, second, "observed", "NONE", "BP", binsize));
    if (!mzd || !mzd->foundFooter) {
        return false;
    }

    int64_t regionIndices[4] = {0, mzd->numBins1, 0, mzd->numBins2};
    const set<int32_t> blockNumbers = mzd->getBlockNumbers(regionIndices);
    for (int32_t blockNumber : blockNumbers) {
        const auto found = mzd->blockMap.find(blockNumber);
        if (found == mzd->blockMap.end()) {
            continue;
        }
        vector<contactRecord> blockRecords = readBlock(impl->stream, found->second, mzd->version);
        if (!blockRecords.empty()) {
            if (transpose) {
                for (auto &record : blockRecords) {
                    std::swap(record.binX, record.binY);
                }
            }
            processor(blockRecords);
        }
    }
    return true;
}

bool forEachRawObservedBlockWithNorm(const string &fileName,
                                     const string &chromosomeName,
                                     int32_t binsize,
                                     const string &norm,
                                     vector<double> &normVector,
                                     const StrawBlockCallback &processor) {
    if (straw_v10::isV10(fileName)) {
        normVector = straw_v10::File(fileName).normalization(chromosomeName, "BP", binsize, norm);
        forEachRawObservedBlock(fileName, chromosomeName, chromosomeName, binsize, processor);
        return true;
    }

    HiCFile hiCFile(fileName);
    const auto chromosomeIt = hiCFile.chromosomeMap.find(chromosomeName);
    if (chromosomeIt == hiCFile.chromosomeMap.end()) {
        cerr << "chromosome " << chromosomeName << " not found in the file." << endl;
        return false;
    }

    MatrixZoomData *mzd = hiCFile.getMatrixZoomData(
        chromosomeName, chromosomeName, "observed", norm, "BP", binsize);
    if (mzd == nullptr || !mzd->foundFooter) {
        delete mzd;
        return false;
    }

    if (norm == "NONE") {
        normVector.assign(static_cast<size_t>(mzd->numBins1) + 1, 1.0);
    } else {
        normVector = mzd->c1Norm;
        if (normVector.empty()) {
            delete mzd;
            return false;
        }
    }

    int64_t regionIndices[4] = {0, mzd->numBins1, 0, mzd->numBins2};
    const set<int32_t> blockNumbers = mzd->getBlockNumbers(regionIndices);
    HiCFileStream stream(fileName);
    for (int32_t blockNumber : blockNumbers) {
        const auto found = mzd->blockMap.find(blockNumber);
        if (found == mzd->blockMap.end()) {
            continue;
        }
        vector<contactRecord> blockRecords = readBlock(stream, found->second, mzd->version);
        if (!blockRecords.empty()) {
            processor(blockRecords);
        }
    }

    delete mzd;
    return true;
}

namespace {
class ScopedCerrSilence {
    ostringstream sink;
    streambuf *previous;
public:
    ScopedCerrSilence() : previous(cerr.rdbuf(sink.rdbuf())) {}
    ~ScopedCerrSilence() { cerr.rdbuf(previous); }
};
}

bool getNormalizationVectorForFile(const string &fileName, const string &chromosomeName,
                                   int32_t binsize, const string &norm,
                                   vector<double> &values, const string &unit) {
    values.clear();
    if (straw_v10::isV10(fileName)) {
        try {
            values = straw_v10::File(fileName).normalization(chromosomeName, unit, binsize, norm);
            return true;
        } catch (const exception &) { return false; }
    }
    HiCFile file(fileName);
    if (!file.chromosomeMap.count(chromosomeName)) return false;
    // Capability probing is expected to fail often (normalizations are
    // resolution-dependent), so do not emit the legacy query path's diagnostics.
    MatrixZoomData *mzd;
    {
        ScopedCerrSilence silence;
        mzd = file.getMatrixZoomData(chromosomeName, chromosomeName,
                                     "observed", norm, unit, binsize);
    }
    if (!mzd || !mzd->foundFooter || (norm != "NONE" && mzd->c1Norm.empty())) {
        delete mzd;
        return false;
    }
    if (norm == "NONE") values.assign(static_cast<size_t>(mzd->numBins1) + 1, 1.0);
    else values = mzd->c1Norm;
    delete mzd;
    return true;
}

bool getExpectedVectorForFile(const string &fileName, const string &chromosomeName,
                              int32_t binsize, const string &norm,
                              vector<double> &values, const string &unit) {
    values.clear();
    if (straw_v10::isV10(fileName)) {
        try {
            values = straw_v10::File(fileName).expected(chromosomeName, unit, binsize, norm);
            return true;
        } catch (const exception &) { return false; }
    }
    HiCFile file(fileName);
    if (!file.chromosomeMap.count(chromosomeName)) return false;
    MatrixZoomData *mzd;
    {
        ScopedCerrSilence silence;
        mzd = file.getMatrixZoomData(chromosomeName, chromosomeName,
                                     "expected", norm, unit, binsize);
    }
    if (!mzd || !mzd->foundFooter || mzd->expectedValues.empty()) {
        delete mzd;
        return false;
    }
    values = mzd->expectedValues;
    delete mzd;
    return true;
}

bool strawStreamRegions(const string &fileName,
                        const string &chromosomeName,
                        int32_t binsize,
                        const string &norm,
                        const vector<StrawRegion> &regions,
                        const StrawRegionRecordCallback &callback) {
    if (straw_v10::isV10(fileName)) {
        straw_v10::File file(fileName);
        for (size_t i = 0; i < regions.size(); ++i) {
            const auto& r = regions[i];
            file.stream("observed", norm,
                        chromosomeName + ":" + to_string(r.xStart) + ":" + to_string(r.xEnd),
                        chromosomeName + ":" + to_string(r.yStart) + ":" + to_string(r.yEnd),
                        "BP", binsize, [&](const contactRecord& record) { callback(i, record); });
        }
        return true;
    }

    HiCFile hiCFile(fileName);
    if (hiCFile.chromosomeMap.count(chromosomeName) == 0) {
        return false;
    }
    MatrixZoomData *mzd = hiCFile.getMatrixZoomData(
        chromosomeName, chromosomeName, "observed", norm, "BP", binsize);
    if (mzd == nullptr || !mzd->foundFooter) {
        delete mzd;
        return false;
    }
    for (size_t regionIndex = 0; regionIndex < regions.size(); ++regionIndex) {
        const StrawRegion &region = regions[regionIndex];
        mzd->streamRecords(region.xStart, region.xEnd, region.yStart, region.yEnd,
                           [&](const contactRecord &record) {
                               contactRecord oriented = record;
                               const bool direct = oriented.binX >= region.xStart &&
                                                   oriented.binX <= region.xEnd &&
                                                   oriented.binY >= region.yStart &&
                                                   oriented.binY <= region.yEnd;
                               if (!direct) {
                                   std::swap(oriented.binX, oriented.binY);
                               }
                               callback(regionIndex, oriented);
                           });
    }
    delete mzd;
    return true;
}

int64_t getNumRecordsForFile(const string &fileName, int32_t binsize, bool interOnly) {
    if (straw_v10::isV10(fileName)) return straw_v10::File(fileName).countRecords(binsize, interOnly);

    HiCFile hiCFile(fileName);
    int64_t totalNumRecords = 0;

    size_t indexOffset = 0;
    if (interOnly){
        indexOffset = 1;
    }

    vector<chromosome> chromosomes = hiCFile.getChromosomes();
    for(size_t i = 0; i < chromosomes.size(); i++){
        if(chromosomes[i].index <= 0) continue;
        for(size_t j = i + indexOffset; j < chromosomes.size(); j++){
            if(chromosomes[j].index <= 0) continue;
            // Previously leaked one MatrixZoomData per chromosome pair.
            unique_ptr<MatrixZoomData> mzd;
            if(chromosomes[i].index > chromosomes[j].index){
                mzd.reset(hiCFile.getMatrixZoomData(chromosomes[j].name, chromosomes[i].name, "observed", "NONE", "BP", binsize));
            } else {
                mzd.reset(hiCFile.getMatrixZoomData(chromosomes[i].name, chromosomes[j].name, "observed", "NONE", "BP", binsize));
            }
            totalNumRecords += mzd->getNumberOfTotalRecords();
        }
    }

    return totalNumRecords;
}

int64_t getNumRecordsForChromosomes(const string &fileName, int32_t binsize, bool interOnly) {
    if (straw_v10::isV10(fileName)) return straw_v10::File(fileName).countRecords(binsize, false, true);

    HiCFile hiCFile(fileName);
    vector<chromosome> chromosomes = hiCFile.getChromosomes();
    for(size_t i = 0; i < chromosomes.size(); i++){
        if(chromosomes[i].index <= 0) continue;
        unique_ptr<MatrixZoomData> mzd(hiCFile.getMatrixZoomData(chromosomes[i].name, chromosomes[i].name, "observed", "NONE", "BP", binsize));
        int64_t totalNumRecords = mzd->getNumberOfTotalRecords();
        cout << chromosomes[i].name << " " << totalNumRecords << " ";
        cout << totalNumRecords*12/1000/1000/1000 << " GB" << endl;
    }
    return 0;
}

vector<pair<string, int64_t>> getRecordCountsByChromosome(const string &fileName, int32_t binsize) {
    if (straw_v10::isV10(fileName)) return straw_v10::File(fileName).countRecordsByChromosome(binsize);

    HiCFile hiCFile(fileName);
    vector<pair<string, int64_t>> result;
    vector<chromosome> chromosomes = hiCFile.getChromosomes();
    for (size_t i = 0; i < chromosomes.size(); i++) {
        if (chromosomes[i].index <= 0) continue;
        unique_ptr<MatrixZoomData> mzd(hiCFile.getMatrixZoomData(chromosomes[i].name, chromosomes[i].name,
                                                                 "observed", "NONE", "BP", binsize));
        result.emplace_back(chromosomes[i].name, mzd->getNumberOfTotalRecords());
    }
    return result;
}

void writeCompressedBuffer(gzFile& file, const char* buffer, size_t size) {
    gzwrite(file, buffer, size);
}

void writeContactRecord(FILE* file, const CompressedContactRecord& record) {
    fwrite(&record, sizeof(CompressedContactRecord), 1, file);
}

void writeCompressedHeader(gzFile& file, const HicSliceHeader& header) {
    gzwrite(file, HICSLICE_MAGIC.c_str(), HICSLICE_MAGIC.length());
    gzwrite(file, (char*)&header.resolution, sizeof(int32_t));
    gzwrite(file, (char*)&header.numChromosomes, sizeof(int32_t));
    for (const auto& chr : header.chromosomeKeys) {
        int32_t nameLength = chr.first.length();
        gzwrite(file, (char*)&nameLength, sizeof(int32_t));
        gzwrite(file, chr.first.c_str(), nameLength);
        gzwrite(file, (char*)&chr.second, sizeof(int16_t));
    }
}

void writeUncompressedHeader(FILE* file, const HicSliceHeader& header) {
    fwrite(HICSLICE_MAGIC.c_str(), 1, HICSLICE_MAGIC.length(), file);
    fwrite(&header.resolution, sizeof(int32_t), 1, file);
    fwrite(&header.numChromosomes, sizeof(int32_t), 1, file);
    for (const auto& chr : header.chromosomeKeys) {
        int32_t nameLength = chr.first.length();
        fwrite(&nameLength, sizeof(int32_t), 1, file);
        fwrite(chr.first.c_str(), 1, nameLength, file);
        fwrite(&chr.second, sizeof(int16_t), 1, file);
    }
}

bool shouldKeepRecord(const contactRecord& rec, const chromosome& chr1, const chromosome& chr2, 
                     int32_t resolution, ContactFilter filter) {
    const int32_t FIVE_MB = 5000000;
    
    if (filter == ContactFilter::ALL) {
        return true;
    }
    
    if (filter == ContactFilter::INTER) {
        return chr1.name != chr2.name;
    }
    
    // For both INTRA cases, first check if same chromosome
    if (chr1.name == chr2.name) {
        if (filter == ContactFilter::INTRA) {
            return true;
        }

        int32_t distance = abs(rec.binX - rec.binY);
        if (filter == ContactFilter::INTRA_SHORT) {
            return distance < FIVE_MB/resolution;
        } else { // INTRA_LONG
            return distance > FIVE_MB/resolution;
        }
    }
    
    return false;  // Different chromosomes for INTRA cases
}

void dumpGenomeWideDataAtResolution(const std::string& matrixType,
                                  const std::string& norm,
                                  const std::string& filePath,
                                  const std::string& unit,
                                  int32_t resolution,
                                  const std::string& outputPath,
                                  bool compressed,
                                  ContactFilter filter) {
    if (straw_v10::isV10(filePath)) {
        dumpV10(matrixType, norm, filePath, unit, resolution, outputPath, compressed, filter);
        return;
    }

    HiCFile hicFile(filePath);

    // Create header
    HicSliceHeader header;
    header.resolution = resolution;

    // Get chromosomes and create mapping
    std::vector<chromosome> chromosomes = hicFile.getChromosomes();
    int16_t chrKey = 0;
    for (const auto& chr : chromosomes) {
        if (chr.index > 0) {
            header.chromosomeKeys[chr.name] = chrKey++;
        }
    }
    header.numChromosomes = static_cast<int32_t>(header.chromosomeKeys.size());

    gzFile gzOut = nullptr;
    FILE *plainOut = nullptr;
    if (compressed) {
        gzOut = gzopen(outputPath.c_str(), "wb");
        if (!gzOut) {
            std::cerr << "Error: Could not open compressed output file " << outputPath << std::endl;
            return;
        }
    } else {
        plainOut = fopen(outputPath.c_str(), "wb");
        if (!plainOut) {
            std::cerr << "Error: Could not open uncompressed output file " << outputPath << std::endl;
            return;
        }
    }

    if (compressed) {
        writeCompressedHeader(gzOut, header);
    } else {
        writeUncompressedHeader(plainOut, header);
    }

    auto writeRecord = [&](const CompressedContactRecord& record) {
        if (compressed) {
            gzwrite(gzOut, (char*)&record, sizeof(CompressedContactRecord));
        } else {
            fwrite(&record, sizeof(CompressedContactRecord), 1, plainOut);
        }
    };

    for (const auto& chr1 : chromosomes) {
        if (chr1.index <= 0) continue;

        for (const auto& chr2 : chromosomes) {
            if (chr2.index <= 0 || chr1.index > chr2.index) continue;

            // Skip chromosome pairs that don't match filter
            if (filter == ContactFilter::INTER && chr1.name == chr2.name) continue;
            if ((filter == ContactFilter::INTRA_SHORT || filter == ContactFilter::INTRA_LONG || filter == ContactFilter::INTRA)
                && chr1.name != chr2.name) continue;

            try {
                MatrixZoomData mzd(chr1, chr2, matrixType, norm, unit, resolution,
                                   hicFile.version, hicFile.master,
                                   HiCFile::totalFileSize, hicFile.fileName);
                if (!mzd.foundFooter) continue;

                // Stream through MatrixZoomData rather than reading blocks
                // directly: this is what applies `norm` and `matrixType`, which
                // the previous block-level loop silently ignored.
                mzd.streamRecords(0, chr1.length, 0, chr2.length,
                                  [&](const contactRecord& rec) {
                    if (!(rec.counts > 0) || isnan(rec.counts) || isinf(rec.counts)) return;
                    // streamRecords yields genomic positions; the slice format
                    // stores bin indices, as the V10 writer also does.
                    const int32_t binX = rec.binX / resolution;
                    const int32_t binY = rec.binY / resolution;
                    contactRecord binned = rec;
                    binned.binX = binX;
                    binned.binY = binY;
                    if (!shouldKeepRecord(binned, chr1, chr2, resolution, filter)) return;

                    // Value-initialized: the struct carries padding, and writing
                    // it uninitialized produced nondeterministic output bytes.
                    CompressedContactRecord compressedRecord{};
                    compressedRecord.chr1Key = header.chromosomeKeys[chr1.name];
                    compressedRecord.binX = binX;
                    compressedRecord.chr2Key = header.chromosomeKeys[chr2.name];
                    compressedRecord.binY = binY;
                    compressedRecord.value = rec.counts;
                    writeRecord(compressedRecord);
                });
            } catch (const std::exception& e) {
                std::cerr << "Skipping chromosome pair " << chr1.name << "-" << chr2.name
                         << ": " << e.what() << std::endl;
            }
        }
    }

    if (compressed) {
        gzclose(gzOut);
    } else {
        fclose(plainOut);
    }
}
