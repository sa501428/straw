#pragma once

#include "straw.h"
#include <zlib.h>

// HBS v1: see HBS_FORMAT.md. Output is staged and renamed only after gzclose.
class HbsWriter {
public:
    HbsWriter(const std::string& input, const std::string& output,
              uint32_t resolution, const std::vector<chromosome>& chromosomes);
    ~HbsWriter();
    HbsWriter(const HbsWriter&) = delete;
    HbsWriter& operator=(const HbsWriter&) = delete;
    // Numeric positions are consumed as supplied in [0, chromosome length];
    // the terminal endpoint aliases the final real matrix bin.
    void record(const chromosome& a, uint64_t x, const chromosome& b, uint64_t y, uint64_t count);
    void finish();
private:
    void bytes(const void* data, unsigned size);
    void number(uint64_t value, unsigned size);
    std::string output_, temporary_;
    gzFile file_ = nullptr;
    uint32_t resolution_;
    std::map<int32_t, uint16_t> ids_;
};

bool isHbsPath(const std::string& path);
