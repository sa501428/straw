#pragma once

#include "straw.h"
#include <memory>

// The V10 wire format is independent of V6-V9. Keep its parser and query path
// out of straw.cpp. See hic-format/HiCFormatV10.md (consolidated specification).
namespace straw_v10 {

// Unlike contactRecord, this interface does not round integer counts to float.
struct Record {
    uint32_t binX, binY;
    uint64_t count = 0;
    float score = 0;
    bool isScore = false;
};
using Callback = std::function<void(const Record &)>;

bool isV10(const std::string &path);

class File {
  public:
    explicit File(const std::string &path);
    ~File();
    File(const File &) = delete;
    File &operator=(const File &) = delete;
    std::vector<chromosome> chromosomes() const;
    std::vector<int32_t> resolutions(const std::string &unit = "BP") const;
    std::string genome() const;
    std::vector<std::string> normalizations() const;
    std::vector<std::pair<std::string, std::string>> attributes() const;
    // Coordinates are bin indices; end positions are exclusive. Raw counts
    // remain uint64_t, including during derived-resolution aggregation.
    void raw(const std::string &chr1, const std::string &chr2, const std::string &unit,
             int32_t resolution, uint64_t x0, uint64_t x1, uint64_t y0, uint64_t y1,
             const Callback &callback);
    // Exact raw query using straw-style chromosome:start:end locations.
    // Returned coordinates remain bin indices, as in raw().
    void streamRaw(const std::string &chr1loc, const std::string &chr2loc, const std::string &unit,
                   int32_t resolution, const Callback &callback);
    std::vector<double> normalization(const std::string &chr, const std::string &unit,
                                      int32_t resolution, const std::string &norm);
    std::vector<double> expected(const std::string &chr, const std::string &unit,
                                 int32_t resolution, const std::string &norm);
    void stream(const std::string &matrixType, const std::string &norm, const std::string &chr1loc,
                const std::string &chr2loc, const std::string &unit, int32_t resolution,
                const StrawRecordCallback &callback);
    std::vector<std::vector<float>> matrix(const std::string &matrixType, const std::string &norm,
                                           const std::string &chr1loc, const std::string &chr2loc,
                                           const std::string &unit, int32_t resolution);
    int64_t countRecords(int32_t resolution, bool interOnly, bool printByChromosome = false);

  private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};

void readHeader(std::istream &input, int64_t &master, std::string &genome, int32_t &nChromosomes,
                int64_t &nvi, int64_t &nviLength, std::map<std::string, chromosome> &chromosomes);
} // namespace straw_v10
