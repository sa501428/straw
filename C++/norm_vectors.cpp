#include "norm_vectors.h"
#include "straw_v10.h"
#include "file_source.h"
#include <array>
#include <map>
#include <tuple>

namespace {
// Buffered absolute cursor for the legacy header/index. Site lists and expected
// arrays are skipped by range; no matrices or contact blocks are decoded.
struct NormCursor {
    straw_file::Source &source;
    uint64_t at = 0, limit, first = 0;
    straw_v10::Bytes buffer;
    explicit NormCursor(straw_file::Source &s) : source(s), limit(s.size) {}
    void need(uint64_t n) const {
        if (at > limit || n > limit - at)
            throw std::runtime_error("truncated normalization header/index");
    }
    uint8_t byte() {
        need(1);
        if (buffer.empty() || at < first || at - first >= buffer.size()) {
            first = at;
            buffer = source.read(at, std::min<uint64_t>(65536, limit - at));
        }
        return buffer[at++ - first];
    }
    uint32_t word() {
        uint32_t out = 0;
        for (unsigned i = 0; i < 4; ++i) out |= uint32_t(byte()) << (8*i);
        return out;
    }
    uint64_t wide() {
        uint64_t lo = word(); return lo | (uint64_t(word()) << 32);
    }
    std::string string() {
        std::string out;
        for (uint8_t b; (b = byte());) {
            if (out.size() == 1048576) throw std::runtime_error("normalization name exceeds limit");
            out.push_back(char(b));
        }
        return out;
    }
    uint32_t count(uint64_t minimum) {
        uint32_t n = word();
        if (n > (limit-at)/minimum) throw std::runtime_error("invalid normalization metadata count");
        return n;
    }
    void skip(uint64_t n) { need(n); at += n; }
};
void normCheck(bool ok, const std::string &message) {
    if (!ok) throw std::runtime_error(message);
}
}

struct StrawNormalizationVectors::Impl {
    straw_file::Source source;
    std::unique_ptr<straw_v10::File> v10;
    std::vector<StrawNormalizationVector> entries;
    std::vector<uint64_t> positions;
    explicit Impl(const std::string &path) : source(path) {
        auto magic = source.read(0, 8);
        straw_v10::Cursor probe(magic);
        probe.magic("HIC\0");
        uint32_t version = probe.word();
        if (version == 10) {
            v10.reset(new straw_v10::File(path));
            entries = v10->normalizationEntries();
            return;
        }
        normCheck(version == 9, "dump-norms supports .hic versions 9 and 10");
        NormCursor c(source); c.at = 8;
        uint64_t footer = c.wide(); c.string();
        uint64_t nvi = c.wide(), nviLength = c.wide();
        uint32_t attrs = c.count(2);
        for (uint32_t i = 0; i < attrs; ++i) { c.string(); c.string(); }
        uint32_t nc = c.count(9);
        std::vector<std::pair<std::string, uint64_t>> chroms;
        for (uint32_t i = 0; i < nc; ++i) {
            auto name = c.string(); uint64_t length = c.wide();
            normCheck(!name.empty() && length <= INT64_MAX, "invalid chromosome in normalization header");
            chroms.emplace_back(std::move(name), length);
        }
        std::array<std::vector<uint32_t>, 2> resolutions;
        for (auto &unit : resolutions) {
            uint32_t nr = c.count(4);
            for (uint32_t i = 0; i < nr; ++i) {
                uint32_t bin = c.word();
                normCheck(bin && bin <= INT32_MAX, "invalid normalization resolution");
                unit.push_back(bin);
            }
        }
        std::vector<uint64_t> fragments(nc);
        if (!resolutions[1].empty()) {
            for (uint32_t i = 0; i < nc; ++i) {
                uint32_t sites = c.count(4);
                fragments[i] = uint64_t(sites)+1;
                c.skip(uint64_t(sites)*4);
            }
        }
        if (!nvi) {
            // V9 files can have zero direct NVI locators. The footer byte count
            // ends after raw expected; normalized expected then precedes NVI.
            c.at = footer;
            uint64_t length = c.wide(); c.skip(length);
            uint32_t ne = c.count(22);
            for (uint32_t i = 0; i < ne; ++i) {
                c.string(); c.string(); c.word();
                uint64_t n = c.wide();
                normCheck(n <= (c.limit-c.at)/4, "invalid normalized expected length");
                c.skip(n*4);
                uint32_t factors = c.count(8); c.skip(uint64_t(factors)*8);
            }
            nvi = c.at;
        }
        normCheck(nvi <= source.size, "normalization index outside file");
        if (nviLength) { source.interval(nvi, nviLength); c.limit = nvi+nviLength; }
        c.at = nvi;
        uint32_t nv = c.count(27);
        using Key = std::tuple<std::string, uint32_t, uint8_t, uint32_t>;
        std::map<Key, std::pair<StrawNormalizationVector, uint64_t>> indexed;
        for (uint32_t i = 0; i < nv; ++i) {
            StrawNormalizationVector entry;
            entry.normalization = c.string(); uint32_t chr = c.word(); entry.unit = c.string();
            uint32_t bin = c.word(); uint64_t pos = c.wide(), length = c.wide();
            normCheck(!entry.normalization.empty() && chr < nc && bin && bin <= INT32_MAX &&
                      (entry.unit == "BP" || entry.unit == "FRAG"), "invalid normalization vector key");
            uint8_t u = entry.unit == "FRAG";
            normCheck(std::find(resolutions[u].begin(), resolutions[u].end(), bin) != resolutions[u].end(),
                      "normalization resolution missing from header");
            source.interval(pos, length);
            normCheck(length >= 8, "truncated normalization vector");
            auto head = source.read(pos, 8);
            straw_v10::Cursor values(head); entry.count = values.wide();
            normCheck(entry.count <= (length-8)/4, "normalization vector exceeds stored length");
            entry.chromosome = chroms[chr].first; entry.resolution = bin;
            uint64_t size = u ? fragments[chr] : chroms[chr].second;
            entry.binCount = size/bin + (size%bin != 0);
            // Historical indexes may contain superseded duplicate keys. Last
            // entry wins, matching the legacy straw lookup and v9 converter.
            Key key{entry.normalization, chr, u, bin};
            indexed[key] = {std::move(entry), pos+8};
        }
        if (nviLength) normCheck(c.at == c.limit, "normalization index length mismatch");
        for (auto &e : indexed) {
            entries.push_back(std::move(e.second.first)); positions.push_back(e.second.second);
        }
    }
    void stream(size_t index, const std::function<void(uint32_t)> &consume) {
        const auto &entry = entries.at(index);
        if (v10) { v10->streamNormalization(entry, consume); return; }
        for (uint64_t begin = 0; begin < entry.count;) {
            uint32_t count = static_cast<uint32_t>(std::min<uint64_t>(65536, entry.count-begin));
            auto bytes = source.read(positions.at(index)+begin*4, uint64_t(count)*4);
            straw_v10::Cursor c(bytes);
            for (uint32_t j = 0; j < count; ++j) consume(c.word());
            begin += count;
        }
    }
};
StrawNormalizationVectors::StrawNormalizationVectors(const std::string &path) : impl(new Impl(path)) {}
StrawNormalizationVectors::~StrawNormalizationVectors() = default;
const std::vector<StrawNormalizationVector> &StrawNormalizationVectors::entries() const { return impl->entries; }
void StrawNormalizationVectors::stream(size_t index, const std::function<void(uint32_t)> &consume) {
    normCheck(bool(consume), "normalization consumer is required");
    impl->stream(index, consume);
}
