#pragma once
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

// Stored vectors only (no synthetic NONE vectors). Works without cis matrices
// or expected arrays. V9 words and V10 decoded words retain their exact bits.
struct StrawNormalizationVector {
    std::string normalization, chromosome, unit;
    int32_t resolution = 0;
    uint64_t count = 0, binCount = 0;
};
class StrawNormalizationVectors {
  public:
    explicit StrawNormalizationVectors(const std::string &path);
    ~StrawNormalizationVectors();
    const std::vector<StrawNormalizationVector> &entries() const;
    void stream(size_t index, const std::function<void(uint32_t)> &consume);
  private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};
int dumpNormsMain(int argc, char **argv);
