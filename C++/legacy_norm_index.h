#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

namespace straw_v9 {
struct NormalizationEntry {
    std::string normalization, unit;
    int32_t chromosome, resolution;
    int64_t position, size;
};

// The direct v9 locator describes the index, not the vector payloads. Parse
// that exact interval without fetching a header from every indexed vector.
inline std::vector<NormalizationEntry> normalizationIndex(const uint8_t *data, size_t length) {
    size_t at = 0;
    auto integer = [&](size_t bytes) {
        if (bytes > length - at) throw std::runtime_error("truncated v9 normalization index");
        uint64_t value = 0;
        for (size_t i = 0; i < bytes; ++i) value |= uint64_t(data[at++]) << (8 * i);
        return value;
    };
    auto string = [&]() {
        size_t first = at;
        while (at < length && data[at]) ++at;
        if (at == length) throw std::runtime_error("unterminated v9 normalization index string");
        std::string value(reinterpret_cast<const char *>(data + first), at - first);
        ++at;
        return value;
    };
    uint64_t count = integer(4);
    if (count > (length - at) / 27) throw std::runtime_error("invalid v9 normalization index count");
    std::vector<NormalizationEntry> entries;
    for (uint64_t i = 0; i < count; ++i) {
        NormalizationEntry entry;
        entry.normalization = string();
        entry.chromosome = static_cast<int32_t>(integer(4));
        entry.unit = string();
        entry.resolution = static_cast<int32_t>(integer(4));
        entry.position = static_cast<int64_t>(integer(8));
        entry.size = static_cast<int64_t>(integer(8));
        if (entry.normalization.empty() || entry.chromosome < 0 || entry.resolution <= 0 ||
            (entry.unit != "BP" && entry.unit != "FRAG") || entry.position < 0 || entry.size < 8)
            throw std::runtime_error("invalid v9 normalization index entry");
        entries.push_back(std::move(entry));
    }
    if (at != length) throw std::runtime_error("v9 normalization index length mismatch");
    return entries;
}

template <class Entry>
bool findNormalizationEntries(const std::vector<NormalizationEntry> &entries, int32_t c1, int32_t c2,
                              const std::string &norm, const std::string &unit, int32_t resolution,
                              Entry &first, Entry &second) {
    bool found1 = false, found2 = false;
    for (const auto &entry : entries) {
        if (entry.normalization != norm || entry.unit != unit || entry.resolution != resolution) continue;
        // Last entry wins for superseded duplicate keys, like the legacy scan.
        if (entry.chromosome == c1) {
            first.position = entry.position; first.size = entry.size; found1 = true;
        }
        if (entry.chromosome == c2) {
            second.position = entry.position; second.size = entry.size; found2 = true;
        }
    }
    return found1 && found2;
}
} // namespace straw_v9
