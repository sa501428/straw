/*
  The MIT License (MIT)

  Copyright (c) 2011-2016 Broad Institute, Aiden Lab

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
#ifndef STRAW_H
#define STRAW_H

#include <fstream>
#include <functional>
#include <set>
#include <vector>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <curl/curl.h>

// Forward declarations
class HiCFile;
class MatrixZoomData;

// pointer structure for reading blocks or matrices, holds the size and position
struct indexEntry {
    int64_t size;
    int64_t position;
};

// sparse matrixType entry
struct contactRecord {
    int32_t binX;
    int32_t binY;
    float counts;
};

// chromosome
struct chromosome {
    std::string name;
    int32_t index;
    int64_t length;
};

// this is for creating a stream from a byte array for ease of use
// see https://stackoverflow.com/questions/41141175/how-to-implement-seekg-seekpos-on-an-in-memory-buffer
// Length is size_t: buffer sizes here come from int64_t index entries, and
// narrowing to int32_t could yield a negative length and egptr() < eback().
struct membuf : std::streambuf {
    membuf(char *begin, std::size_t l) {
        setg(begin, begin, begin + l);
    }
};

struct memstream : virtual membuf, std::istream {
    memstream(char *begin, std::size_t l) :
            membuf(begin, l),
            std::istream(static_cast<std::streambuf*>(this)) {
    }

    std::istream::pos_type seekpos(std::istream::pos_type sp, std::ios_base::openmode which) override {
        return seekoff(sp - std::istream::pos_type(std::istream::off_type(0)), std::ios_base::beg, which);
    }

    // Every seek is clamped to [eback(), egptr()]. Callers routinely skip over
    // sections using sizes read from the file, so an out-of-range offset must
    // saturate rather than move gptr outside the buffer.
    std::istream::pos_type seekoff(std::istream::off_type off,
                                    std::ios_base::seekdir dir,
                                    std::ios_base::openmode which = std::ios_base::in) override {
        (void)which;
        const std::istream::off_type size = egptr() - eback();
        std::istream::off_type target;
        if (dir == std::ios_base::cur)
            target = (gptr() - eback()) + off;
        else if (dir == std::ios_base::end)
            target = size + off;
        else
            target = off;
        if (target < 0) target = 0;
        if (target > size) target = size;
        setg(eback(), eback() + target, egptr());
        return target;
    }
};

// for holding data from URL call
struct MemoryStruct {
    char *memory;
    size_t size;
};

// Function declarations
std::vector<contactRecord> straw(const std::string& matrixType, const std::string& norm, const std::string& fname, 
                               const std::string& chr1loc, const std::string& chr2loc, const std::string& unit, 
                               int32_t binsize);

using StrawRecordCallback = std::function<void(const contactRecord&)>;

enum class StrawErrorCode {
    InvalidArgument,
    Io,
    UnsupportedVersion,
    CorruptFile,
    NotFound,
    Unavailable,
    UnsupportedOperation,
    Internal
};

class StrawException : public std::runtime_error {
public:
    StrawException(StrawErrorCode code, const std::string &message);
    StrawErrorCode code() const noexcept;
private:
    StrawErrorCode errorCode;
};

// Retains the parsed file header and prepared matrix/zoom state across window
// calls. This is the native state held by the stable C straw_query_t handle.
class StrawPreparedQuery {
public:
    StrawPreparedQuery(const std::string &fileName,
                       const std::string &matrixType,
                       const std::string &normalization,
                       const std::string &firstChromosome,
                       const std::string &secondChromosome,
                       const std::string &unit,
                       int32_t resolution);
    ~StrawPreparedQuery();
    StrawPreparedQuery(const StrawPreparedQuery &) = delete;
    StrawPreparedQuery &operator=(const StrawPreparedQuery &) = delete;

    void streamWindow(int64_t xStart, int64_t xEnd, int64_t yStart, int64_t yEnd,
                      const StrawRecordCallback &callback);

private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};

bool strawStream(const std::string& matrixType, const std::string& norm, const std::string& fname,
                 const std::string& chr1loc, const std::string& chr2loc, const std::string& unit,
                 int32_t binsize, const StrawRecordCallback& callback);

std::vector<std::vector<float>> strawAsMatrix(const std::string& matrixType, const std::string& norm, 
                                            const std::string& fileName, const std::string& chr1loc, 
                                            const std::string& chr2loc, const std::string& unit, 
                                            int32_t binsize);

using StrawBlockCallback = std::function<void(const std::vector<contactRecord>&)>;

std::vector<chromosome> getChromosomesForFile(const std::string& fileName);

std::vector<int32_t> getResolutionsForFile(const std::string& fileName,
                                           const std::string& unit = "BP");

std::string getGenomeForFile(const std::string& fileName);

int32_t getVersionForFile(const std::string& fileName);

std::vector<std::string> getNormalizationsForFile(const std::string& fileName);

std::vector<std::pair<std::string, std::string>> getAttributesForFile(
    const std::string& fileName);

// Everything the per-file accessors above return, gathered from a single open.
// Calling them individually costs two file opens each (an isV10 probe plus a
// reader), which over HTTP is two round trips per property.
struct StrawFileInfo {
    int32_t version = 0;
    std::string genome;
    std::vector<chromosome> chromosomes;
    std::vector<int32_t> bpResolutions;
    std::vector<int32_t> fragResolutions;
    std::vector<std::string> normalizations;
    std::vector<std::pair<std::string, std::string>> attributes;
};

StrawFileInfo getFileInfo(const std::string& fileName);

void forEachRawObservedBlock(const std::string& fileName,
                             const std::string& chr1,
                             const std::string& chr2,
                             int32_t binsize,
                             const StrawBlockCallback& processor);

// Streams raw observed records while also returning the normalization vector
// loaded from the same MatrixZoomData.  This lets callers retain integer-like
// raw counts for statistical tests and derive normalized values without a
// second pass through the contact blocks.
bool forEachRawObservedBlockWithNorm(const std::string& fileName,
                                     const std::string& chromosomeName,
                                     int32_t binsize,
                                     const std::string& norm,
                                     std::vector<double>& normVector,
                                     const StrawBlockCallback& processor);

// Reads raw observed blocks for many chromosome pairs while holding the parsed
// header and a single open stream. forEachRawObservedBlock re-opens the file and
// re-parses the header for every call, which costs one round trip per pair over
// HTTP; callers sweeping the whole genome should use this instead.
class StrawRawReader {
public:
    explicit StrawRawReader(const std::string& fileName);
    ~StrawRawReader();
    StrawRawReader(const StrawRawReader&) = delete;
    StrawRawReader& operator=(const StrawRawReader&) = delete;

    // Returns false when the pair or resolution is unavailable. Records are
    // returned in the caller's chromosome order, as forEachRawObservedBlock does.
    bool forEachBlock(const std::string& chr1, const std::string& chr2, int32_t binsize,
                      const StrawBlockCallback& processor);

private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};

// Load a complete vector without reading contact blocks. Returns false when
// the requested capability is not present in the file.
bool getNormalizationVectorForFile(const std::string& fileName,
                                   const std::string& chromosomeName,
                                   int32_t binsize,
                                   const std::string& norm,
                                   std::vector<double>& values,
                                   const std::string& unit = "BP");

bool getExpectedVectorForFile(const std::string& fileName,
                              const std::string& chromosomeName,
                              int32_t binsize,
                              const std::string& norm,
                              std::vector<double>& values,
                              const std::string& unit = "BP");

struct StrawRegion {
    int64_t xStart;
    int64_t xEnd;
    int64_t yStart;
    int64_t yEnd;
};

using StrawRegionRecordCallback = std::function<void(size_t, const contactRecord&)>;

// Query many regions from one matrix setup. Records are normalized according
// to `norm`, and coordinates are returned in base pairs, as in strawStream.
bool strawStreamRegions(const std::string& fileName,
                        const std::string& chromosomeName,
                        int32_t binsize,
                        const std::string& norm,
                        const std::vector<StrawRegion>& regions,
                        const StrawRegionRecordCallback& callback);

int64_t getNumRecordsForFile(const std::string& filename, int32_t binsize, bool interOnly);

int64_t getNumRecordsForChromosomes(const std::string& filename, int32_t binsize, bool interOnly);

std::vector<std::pair<std::string, int64_t>> getRecordCountsByChromosome(
    const std::string& filename, int32_t binsize);

// Add readHeader declaration
std::map<std::string, chromosome> readHeader(std::istream &fin, int64_t &masterIndexPosition, 
                                           std::string &genomeID, int32_t &numChromosomes, 
                                           int32_t &version, int64_t &nviPosition, 
                                           int64_t &nviLength,
                                           std::vector<std::pair<std::string, std::string>> *attributes = nullptr);

#endif
