#include "straw_v10.h"
#include "v10_binary.h"
#include "file_source.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <mutex>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_set>
#include <zstd.h>

namespace straw_v10 {
namespace {
using straw_file::Source;
struct Locator {
    uint64_t pos = 0, len = 0;
};
struct Resolution {
    uint32_t bin, source;
    uint8_t mode, aggregation;
};
struct Header {
    uint64_t length;
    Locator footer, norm, expected, normExpected;
    std::string genome;
    std::vector<std::pair<std::string, std::string>> attributes;
    std::vector<chromosome> chroms;
    std::array<std::vector<Resolution>, 2> resolutions;
    std::vector<uint64_t> fragments;
    std::vector<std::string> norms;
    uint64_t bins(uint32_t chr, uint8_t unit, uint32_t ri) const {
        uint64_t length = unit ? fragments.at(chr) : uint64_t(chroms.at(chr).length);
        uint32_t r = resolutions.at(unit).at(ri).bin;
        return length / r + (length % r != 0);
    }
};
Locator locator(Cursor &c) {
    Locator l;
    l.pos = c.wide();
    l.len = c.wide();
    require((l.pos == 0) == (l.len == 0), "incomplete locator");
    return l;
}
Header parseHeader(const Bytes &bytes) {
    Cursor c(bytes);
    Header h;
    c.magic("HIC\0");
    require(c.word() == 10, "unsupported file version");
    h.length = c.wide();
    require(h.length == bytes.size() && h.length >= 88, "header length mismatch");
    h.footer = locator(c);
    require(h.footer.pos && h.footer.len, "missing footer");
    h.norm = locator(c);
    h.expected = locator(c);
    h.normExpected = locator(c);
    c.zero(8);
    h.genome = c.str();
    uint32_t n = c.word();
    require(n <= c.left() / 2, "attribute count out of bounds");
    for (uint32_t i = 0; i < n; ++i) {
        auto key = c.str();
        auto value = c.str();
        h.attributes.emplace_back(key, value);
    }
    n = c.word();
    require(n > 0 && n <= c.left() / 10 && n <= INT32_MAX, "chromosome count out of bounds");
    std::set<std::string> names;
    for (uint32_t i = 0; i < n; ++i) {
        auto name = c.str();
        auto length = c.wide();
        require(!name.empty() && names.insert(name).second && length > 0 && length <= INT64_MAX,
                "invalid chromosome");
        h.chroms.push_back({name, static_cast<int32_t>(i), static_cast<int64_t>(length)});
    }
    for (auto &list : h.resolutions) {
        n = c.word();
        require(n <= c.left() / 12, "resolution count out of bounds");
        for (uint32_t i = 0; i < n; ++i) {
            Resolution r;
            r.bin = c.word();
            r.mode = c.byte();
            r.aggregation = c.byte();
            c.zero(2);
            r.source = c.word();
            require(r.bin > 0 && r.bin <= INT32_MAX && (list.empty() || r.bin > list.back().bin),
                    "invalid resolution order");
            require(r.mode <= 1 && r.aggregation <= 1, "unknown resolution enumeration");
            if (!r.mode)
                require(r.source == UINT32_MAX, "materialized source must be absent");
            else
                require(r.aggregation == 1 && r.source < i && list[r.source].mode == 0 &&
                            r.bin % list[r.source].bin == 0,
                        "invalid derived source");
            list.push_back(r);
        }
    }
    h.fragments.assign(h.chroms.size(), 0);
    if (!h.resolutions[1].empty())
        for (size_t i = 0; i < h.chroms.size(); ++i) {
            n = c.word();
            require(n <= c.left() / 8, "fragment-site count out of bounds");
            uint64_t prev = 0;
            for (uint32_t j = 0; j < n; ++j) {
                auto site = c.wide();
                require(site > prev && site < uint64_t(h.chroms[i].length),
                        "invalid fragment site");
                prev = site;
            }
            h.fragments[i] = uint64_t(n) + 1;
        }
    n = c.word();
    require(n <= c.left() / 2, "normalization count out of bounds");
    names.clear();
    for (uint32_t i = 0; i < n; ++i) {
        auto name = c.str();
        require(!name.empty() && name != "NONE" && names.insert(name).second,
                "invalid normalization name");
        h.norms.push_back(name);
    }
    c.done();
    for (uint8_t unit = 0; unit < 2; ++unit)
        for (uint32_t ri = 0; ri < h.resolutions[unit].size(); ++ri)
            for (uint32_t chr = 0; chr < h.chroms.size(); ++chr)
                require(h.bins(chr, unit, ri) <= UINT32_MAX, "chromosome bin count exceeds uint32");
    return h;
}
Bytes decompress(ZSTD_DCtx *decoder, Cursor &c, uint32_t bytes) {
    require(bytes > 0 && bytes <= allocationLimit, "decompressed record exceeds allocation limit");
    c.need(4);
    require(c.p[c.at] == 0x28 && c.p[c.at + 1] == 0xb5 && c.p[c.at + 2] == 0x2f &&
                c.p[c.at + 3] == 0xfd,
            "not a Zstandard data frame");
    require(ZSTD_getDictID_fromFrame(c.p + c.at, c.left()) == 0,
            "Zstandard dictionaries are forbidden");
    size_t frame = ZSTD_findFrameCompressedSize(c.p + c.at, c.left());
    require(!ZSTD_isError(frame) && frame == c.left(), "invalid or concatenated Zstandard frame");
    Bytes result(bytes);
    size_t n = ZSTD_decompressDCtx(decoder, result.data(), result.size(), c.p + c.at, c.left());
    require(!ZSTD_isError(n) && n == bytes, "Zstandard decompression/length/checksum failure");
    return result;
}
uint8_t unitId(const std::string &unit) {
    require(unit == "BP" || unit == "FRAG", "unit must be BP or FRAG");
    return unit == "FRAG";
}
struct Zoom {
    uint8_t unit, mode, aggregation, type, grid;
    uint32_t ri, bin, source, B, columns, blocks;
    uint64_t sum, occupied;
    Locator index;
};
struct BlockEntry {
    uint32_t number, len;
    uint64_t pos;
};
using Pair = std::pair<uint32_t, uint32_t>;
uint32_t distanceDepth(uint64_t d, uint32_t blockBinCount) {
    uint32_t depth = 0;
    using Wide = unsigned __int128;
    Wide lhs = Wide(d) * d, scale = Wide(2) * blockBinCount * blockBinCount;
    while (depth < 32) {
        Wide threshold = (Wide(1) << (depth + 1)) - 1;
        if (threshold * threshold > lhs / scale)
            break;
        ++depth;
    }
    return depth;
}

uint32_t blockNumber(uint32_t x, uint32_t y, const Zoom &z) {
    if (!z.grid)
        return u32(uint64_t(y / z.B) * z.columns + x / z.B);
    uint32_t depth = distanceDepth(uint64_t(y) - x, z.B);
    return u32(uint64_t(depth) * z.columns +
               (uint64_t(x) + y) / (uint64_t(2) * z.B));
}

using BlockRange = std::pair<uint32_t, uint32_t>;

std::vector<BlockRange> candidateBlockRanges(const Zoom &z, uint64_t x0, uint64_t x1,
                                             uint64_t y0, uint64_t y1) {
    std::vector<BlockRange> ranges;
    if (x0 >= x1 || y0 >= y1)
        return ranges;

    auto append = [&](uint64_t lo, uint64_t hi) {
        if (lo > UINT32_MAX)
            return;
        ranges.emplace_back(static_cast<uint32_t>(lo),
                            static_cast<uint32_t>(std::min<uint64_t>(hi, UINT32_MAX)));
    };

    if (!z.grid) {
        uint64_t firstColumn = x0 / z.B, lastColumn = (x1 - 1) / z.B;
        uint64_t firstRow = y0 / z.B, lastRow = (y1 - 1) / z.B;
        ranges.reserve(static_cast<size_t>(lastRow - firstRow + 1));
        for (uint64_t row = firstRow; row <= lastRow; ++row)
            append(row * z.columns + firstColumn, row * z.columns + lastColumn);
        return ranges;
    }

    // Cis blocks use V9's rotated (diagonal, distance) grid.  The old V10
    // reader searched every distance band from zero to the farthest corner,
    // causing far-cis queries to decompress dense near-diagonal blocks.  Map
    // the half-open query rectangle to its exact enclosing range on both
    // rotated axes instead.
    uint64_t lastX = x1 - 1, lastY = y1 - 1;
    uint64_t nearest = 0;
    if (x1 <= y0)
        nearest = y0 - lastX;
    else if (y1 <= x0)
        nearest = x0 - lastY;
    uint64_t farthest = std::max(lastY >= x0 ? lastY - x0 : x0 - lastY,
                                 lastX >= y0 ? lastX - y0 : y0 - lastX);
    uint32_t firstDepth = distanceDepth(nearest, z.B);
    uint32_t lastDepth = distanceDepth(farthest, z.B);
    uint64_t firstPad = (x0 + y0) / (uint64_t(2) * z.B);
    uint64_t lastPad = (lastX + lastY) / (uint64_t(2) * z.B);
    ranges.reserve(uint64_t(lastDepth) - firstDepth + 1);
    for (uint64_t depth = firstDepth; depth <= lastDepth; ++depth)
        append(depth * z.columns + firstPad, depth * z.columns + lastPad);
    return ranges;
}
} // namespace

struct File::Impl {
    Source source;
    std::unique_ptr<ZSTD_DCtx, size_t (*)(ZSTD_DCtx *)> decoder;
    Header h;
    Bytes normIndex;
    std::map<Pair, Locator> matrices;
    std::map<Pair, std::vector<Zoom>> zoomCache;
    // Block indexes are re-read and re-validated on every query otherwise, which
    // a windowed reader (StrawPreparedQuery) pays once per window.
    std::map<uint64_t, std::vector<BlockEntry>> blockIndexCache;
    explicit Impl(const std::string &path)
        : source(path), decoder(ZSTD_createDCtx(), ZSTD_freeDCtx) {
        require(decoder != nullptr, "cannot allocate Zstandard decoder");
        auto prefix = source.read(0, 88);
        Cursor c(prefix);
        c.magic("HIC\0");
        require(c.word() == 10, "not a V10 file");
        auto len = c.wide();
        h = parseHeader(source.read(0, len));
        for (auto loc : {h.footer, h.norm, h.expected, h.normExpected})
            if (loc.len)
                source.interval(loc.pos, loc.len);
        auto data = source.read(h.footer.pos, h.footer.len);
        Cursor f(data);
        f.magic("H10F");
        require(f.word() == 1, "unknown footer version");
        require(f.wide() == h.footer.len, "footer length mismatch");
        uint32_t n = f.word();
        f.zero(4);
        require(uint64_t(n) * 24 == f.left(), "invalid footer count");
        Pair previous{};
        for (uint32_t i = 0; i < n; ++i) {
            uint32_t a = f.word(), b = f.word();
            Pair key{a, b};
            auto loc = locator(f);
            require(a <= b && b < h.chroms.size() && (!i || previous < key),
                    "invalid matrix directory key");
            source.interval(loc.pos, loc.len);
            matrices.emplace(key, loc);
            previous = key;
        }
    }
    uint32_t chromosomeId(const std::string &name) const {
        for (size_t i = 0; i < h.chroms.size(); ++i)
            if (h.chroms[i].name == name)
                return static_cast<uint32_t>(i);
        throw std::runtime_error("V10: unknown chromosome " + name);
    }
    uint32_t resolutionId(uint8_t unit, int32_t bin) const {
        for (size_t i = 0; i < h.resolutions[unit].size(); ++i)
            if (h.resolutions[unit][i].bin == uint32_t(bin))
                return static_cast<uint32_t>(i);
        throw std::runtime_error("V10: unavailable resolution " + std::to_string(bin));
    }
    const std::vector<Zoom> &zooms(Pair key) {
        auto found = zoomCache.find(key);
        if (found != zoomCache.end())
            return found->second;
        std::vector<Zoom> result;
        auto m = matrices.find(key);
        if (m == matrices.end())
            return zoomCache.emplace(key, result).first->second;
        auto data = source.read(m->second.pos, m->second.len);
        Cursor c(data);
        c.magic("H10M");
        require(c.word() == 1, "unknown matrix version");
        require(c.word() == key.first && c.word() == key.second, "matrix key mismatch");
        uint32_t n = c.word();
        c.zero(4);
        require(n == h.resolutions[0].size() + h.resolutions[1].size() &&
                    uint64_t(n) * 76 == c.left(),
                "matrix descriptor count/length mismatch");
        for (uint32_t i = 0; i < n; ++i) {
            Zoom z;
            z.unit = c.byte();
            z.mode = c.byte();
            z.aggregation = c.byte();
            z.type = c.byte();
            z.ri = c.word();
            z.bin = c.word();
            z.source = c.word();
            z.grid = c.byte();
            c.zero(3);
            z.sum = c.wide();
            z.occupied = c.wide();
            c.word();
            c.word();
            z.B = c.word();
            z.columns = c.word();
            z.index = locator(c);
            z.blocks = c.word();
            c.zero(4);
            uint8_t unit = i < h.resolutions[0].size() ? 0 : 1;
            uint32_t ri = i - (unit ? h.resolutions[0].size() : 0);
            const auto &r = h.resolutions[unit][ri];
            require(z.unit == unit && z.ri == ri && z.bin == r.bin && z.mode == r.mode &&
                        z.aggregation == r.aggregation && z.source == r.source,
                    "resolution descriptor mismatch");
            require(z.type <= 1 && z.grid == uint8_t(key.first == key.second) && z.B && z.columns,
                    "invalid matrix geometry/type");
            require(z.columns == (h.bins(key.first, unit, ri) + z.B - 1) / z.B,
                    "invalid block column count");
            if (z.mode)
                require(!z.index.pos && !z.blocks, "derived resolution has storage");
            else if (z.occupied)
                require(z.index.pos && z.blocks, "missing matrix blocks");
            else
                require(!z.index.pos && !z.blocks, "empty matrix has block storage");
            if (z.index.len)
                source.interval(z.index.pos, z.index.len);
            result.push_back(z);
        }
        for (const auto &z : result)
            if (z.mode)
                require(z.type == result[(z.unit ? h.resolutions[0].size() : 0) + z.source].type,
                        "derived value type mismatch");
        return zoomCache.emplace(key, std::move(result)).first->second;
    }
    const Zoom *zoom(Pair key, uint8_t unit, uint32_t ri) {
        const auto &all = zooms(key);
        if (all.empty())
            return nullptr;
        return &all.at((unit ? h.resolutions[0].size() : 0) + ri);
    }
    const std::vector<BlockEntry> &blockIndex(const Zoom &z) {
        static const std::vector<BlockEntry> empty;
        if (!z.index.len)
            return empty;
        auto cached = blockIndexCache.find(z.index.pos);
        if (cached != blockIndexCache.end())
            return cached->second;
        auto data = source.read(z.index.pos, z.index.len);
        Cursor c(data);
        c.magic("H10I");
        require(c.word() == 2 && c.wide() == z.index.len,
                "block index version/length mismatch");
        uint32_t count = c.word();
        c.zero(4);
        require(count == z.blocks && z.index.len == 24 + uint64_t(count) * 16,
                "invalid block index count/length");
        std::vector<BlockEntry> out;
        out.reserve(count);
        for (uint32_t i = 0; i < count; ++i) {
            BlockEntry entry{c.word(), c.word(), c.wide()};
            require(!i || entry.number > out.back().number, "unordered block index");
            require(entry.len > 16, "invalid stored block length");
            source.interval(entry.pos, entry.len);
            out.push_back(entry);
        }
        c.done();
        return blockIndexCache.emplace(z.index.pos, std::move(out)).first->second;
    }
    template <typename Emit>
    void block(Cursor c, uint32_t number, const Zoom &z, Pair key, Emit &&cb) {
        require(c.byte() == 1, "unknown block version");
        uint8_t rep = c.byte(), mode = c.byte(), type = c.byte(), flags = c.byte();
        c.zero(3);
        uint32_t x = c.word(), y = c.word(), w = c.word(), height = c.word();
        uint64_t n = c.wide();
        uint32_t np = c.word(), nv = c.word();
        require(rep <= 2 && mode <= 2 && type == z.type && flags <= 1 && w && height && n,
                "invalid block header");
        uint64_t cells = uint64_t(w) * height, slots = rep == 2 ? cells : n;
        require(n <= cells && slots <= allocationLimit / sizeof(uint64_t) &&
                    uint64_t(np) + nv == c.left(),
                "invalid block stream lengths/size");
        require(x < h.bins(key.first, z.unit, z.ri) && y < h.bins(key.second, z.unit, z.ri),
                "block offsets exceed chromosome");
        Cursor positions = c.take(np), values = c.take(nv);
        if (rep == 0) {
            require(!flags && n <= np, "invalid sparse position stream");
        } else if (rep == 1 || type == 1) {
            require(flags == 1 && np == (cells + 7) / 8, "invalid presence bitmap");
            if (cells % 8)
                require((positions.p[np - 1] >> (cells % 8)) == 0, "nonzero bitmap padding");
            positions.at = positions.size;
        } else
            require(!flags && !np, "dense counts have no presence stream");
        require(rep != 2 || mode == 2, "dense values must be direct");
        auto scalar = [&]() { return type ? uint64_t(values.word()) : values.var(); };
        uint64_t defaultValue = 0;
        std::vector<uint64_t> exceptionOrdinals;
        if (mode == 0) {
            defaultValue = scalar();
        } else if (mode == 1) {
            defaultValue = scalar();
            uint64_t ne = values.var();
            require(ne > 0 && ne < slots && ne <= values.left(), "invalid exception count");
            exceptionOrdinals.reserve(static_cast<size_t>(ne));
            uint64_t prev = 0;
            for (uint64_t i = 0; i < ne; ++i) {
                auto d = values.var();
                require(!i || d, "duplicate exception ordinal");
                auto o = i ? add(prev, d) : d;
                require(o < slots, "exception out of range");
                exceptionOrdinals.push_back(o);
                prev = o;
            }
        } else
            require(slots <= values.left() / (type ? 4 : 1), "truncated values");

        size_t exceptionIndex = 0;
        auto valueAt = [&](uint64_t ordinal) {
            if (mode == 0)
                return defaultValue;
            if (mode == 1) {
                if (exceptionIndex < exceptionOrdinals.size() &&
                    exceptionOrdinals[exceptionIndex] == ordinal) {
                    ++exceptionIndex;
                    auto exception = scalar();
                    require(exception != defaultValue, "exception equals default");
                    return exception;
                }
                return defaultValue;
            }
            return scalar();
        };

        uint64_t emitted = 0;
        auto emit = [&](uint64_t pos, uint64_t v, bool present) {
            if (!present)
                return;
            if (rep != 2 && !type)
                require(v > 0, "sparse/bitmap count must be positive");
            uint32_t bx = u32(uint64_t(x) + pos % w), by = u32(uint64_t(y) + pos / w);
            require(bx < h.bins(key.first, z.unit, z.ri) && by < h.bins(key.second, z.unit, z.ri) &&
                        (key.first != key.second || by >= bx) && blockNumber(bx, by, z) == number,
                    "cell violates block geometry");
            Record r{bx, by, v, type ? asFloat(static_cast<uint32_t>(v)) : 0.0f, type != 0};
            cb(r);
            ++emitted;
        };

        uint64_t consumed = 0;
        if (rep == 0) {
            uint64_t previous = 0;
            for (uint64_t i = 0; i < n; ++i) {
                auto delta = positions.var();
                require(!i || delta, "duplicate sparse cell");
                auto pos = i ? add(previous, delta) : delta;
                require(pos < cells, "sparse position out of bounds");
                emit(pos, valueAt(i), true);
                previous = pos;
                ++consumed;
            }
        } else if (rep == 1) {
            for (uint64_t pos = 0; pos < cells; ++pos) {
                if (!(positions.p[pos / 8] & (1u << (pos % 8))))
                    continue;
                require(consumed < n, "bitmap population mismatch");
                emit(pos, valueAt(consumed), true);
                ++consumed;
            }
            require(consumed == n, "bitmap population mismatch");
        } else {
            for (uint64_t pos = 0; pos < cells; ++pos) {
                uint64_t v = valueAt(pos);
                bool present = type ? (positions.p[pos / 8] & (1u << (pos % 8))) != 0 : v != 0;
                if (type && !present)
                    require(v == 0, "absent dense score must be positive zero");
                emit(pos, v, present);
                ++consumed;
            }
        }
        positions.done();
        require(consumed == slots && exceptionIndex == exceptionOrdinals.size(),
                "value slot count mismatch");
        values.done();
        require(emitted == n, "occupied cell count mismatch");
    }
    template <typename Emit>
    void materialized(Pair key, const Zoom &z, uint64_t x0, uint64_t x1, uint64_t y0, uint64_t y1,
                      Emit &&cb) {
        uint64_t totalBlocks = 0, totalCells = 0, sum = 0;
        const auto &index = blockIndex(z);
        const auto ranges = candidateBlockRanges(z, x0, x1, y0, y1);
        for (const auto &range : ranges) {
            auto entry = std::lower_bound(index.begin(), index.end(), range.first,
                                          [](const BlockEntry &a, uint32_t number) {
                                              return a.number < number;
                                          });
            for (; entry != index.end() && entry->number <= range.second; ++entry) {
                auto bytes = source.read(entry->pos, entry->len);
                Cursor c(bytes);
                c.magic("H10B");
                require(c.byte() == 1 && c.byte() == 1, "unknown block codec/record version");
                c.zero(2);
                uint32_t raw = c.word();
                require(c.word() == entry->number && raw >= 40, "block record/index mismatch");
                auto payload = decompress(decoder.get(), c, raw);
                block(Cursor(payload), entry->number, z, key, [&](const Record &r) {
                    ++totalCells;
                    if (!z.type)
                        sum = add(sum, r.count);
                    cb(r);
                });
                ++totalBlocks;
            }
        }
        if (totalBlocks == index.size())
            require(totalBlocks == z.blocks && totalCells == z.occupied && (z.type || sum == z.sum),
                    "matrix statistics/block count mismatch");
    }
    template <typename Emit>
    void raw(Pair key, uint8_t unit, uint32_t ri, uint64_t x0, uint64_t x1, uint64_t y0,
             uint64_t y1, Emit &&cb) {
        const Zoom *z = zoom(key, unit, ri);
        if (!z)
            return;
        auto inside = [&](uint32_t x, uint32_t y) {
            return (x >= x0 && x < x1 && y >= y0 && y < y1) ||
                   (key.first == key.second && y >= x0 && y < x1 && x >= y0 && x < y1);
        };
        if (!z->mode) {
            materialized(key, *z, x0, x1, y0, y1, [&](const Record &r) {
                if (inside(r.binX, r.binY))
                    cb(r);
            });
            return;
        }
        const Zoom &source = *zoom(key, unit, z->source);
        const uint64_t factor = z->bin / source.bin;
        constexpr size_t limit = 65536;
        struct Split {};
        struct Acc { uint64_t count = 0; double score = 0; };
        uint64_t emitted = 0, total = 0;
        auto publish = [&](const Record &r) {
            ++emitted;
            if (!r.isScore) total = add(total, r.count);
            cb(r);
        };
        // For an unusually large score aggregation factor, even one target
        // pixel may exceed the buffer. Traverse source rows in bounded batches.
        std::function<void(uint64_t,uint64_t,uint64_t,uint64_t,double&)> orderedScore;
        orderedScore = [&](uint64_t a, uint64_t b, uint64_t c, uint64_t d, double &sum) {
            std::vector<Record> records;
            try {
                materialized(key, source, a, b, c, d, [&](const Record &r) {
                    if (r.binX < a || r.binX >= b || r.binY < c || r.binY >= d) return;
                    if (records.size() == limit) throw Split{};
                    records.push_back(r);
                });
            } catch (const Split &) {
                std::vector<Record>().swap(records);
                if (d - c > 1) {
                    const auto mid = c + (d - c) / 2;
                    orderedScore(a, b, c, mid, sum); orderedScore(a, b, mid, d, sum);
                } else {
                    require(b - a > 1, "duplicate source cells");
                    const auto mid = a + (b - a) / 2;
                    orderedScore(a, mid, c, d, sum); orderedScore(mid, b, c, d, sum);
                }
                return;
            }
            std::sort(records.begin(), records.end(), [](const Record &a, const Record &b) {
                return std::tie(a.binY, a.binX) < std::tie(b.binY, b.binX);
            });
            for (const auto &r : records) {
                require(std::isfinite(r.score), "nonfinite derived source score");
                sum += double(r.score);
                require(std::isfinite(sum), "derived score overflow");
            }
        };
        std::function<void(uint64_t,uint64_t,uint64_t,uint64_t)> visit;
        visit = [&](uint64_t a, uint64_t b, uint64_t c, uint64_t d) {
            std::map<std::pair<uint32_t,uint32_t>, Acc> sums;
            std::vector<Record> scores;
            try {
                materialized(key, source, a * factor, std::min(b * factor, h.bins(key.first, unit, source.ri)),
                             c * factor, std::min(d * factor, h.bins(key.second, unit, source.ri)),
                             [&](const Record &r) {
                    const uint32_t x = r.binX / factor, y = r.binY / factor;
                    // Assign a cis cell to exactly one query tile even when
                    // both its canonical and reflected positions are requested.
                    uint32_t px = x, py = y;
                    if (!(x >= x0 && x < x1 && y >= y0 && y < y1)) {
                        if (key.first != key.second) return;
                        std::swap(px, py);
                    }
                    if (px < a || px >= b || py < c || py >= d) return;
                    if (z->type) {
                        if (scores.size() == limit) throw Split{};
                        scores.push_back(r);
                    } else {
                        auto &acc = sums[{y, x}];
                        acc.count = add(acc.count, r.count);
                        if (sums.size() > limit) throw Split{};
                    }
                });
            } catch (const Split &) {
                // Release the parent buffers before descending.
                sums.clear(); std::vector<Record>().swap(scores);
                if (d - c > 1) {
                    const auto mid = c + (d - c) / 2;
                    visit(a, b, c, mid); visit(a, b, mid, d);
                } else if (b - a > 1) {
                    const auto mid = a + (b - a) / 2;
                    visit(a, mid, c, d); visit(mid, b, c, d);
                } else {
                    uint64_t x = a, y = c;
                    if (key.first == key.second && x > y) std::swap(x, y);
                    double sum = 0;
                    orderedScore(x * factor, std::min((x + 1) * factor, h.bins(key.first, unit, source.ri)),
                                 y * factor, std::min((y + 1) * factor, h.bins(key.second, unit, source.ri)), sum);
                    publish(Record{u32(x), u32(y), 0, static_cast<float>(sum), true});
                }
                return;
            }
            std::sort(scores.begin(), scores.end(), [](const Record &a, const Record &b) {
                return std::tie(a.binY, a.binX) < std::tie(b.binY, b.binX);
            });
            for (const auto &r : scores) {
                require(std::isfinite(r.score), "nonfinite derived source score");
                auto &acc = sums[{r.binY / factor, r.binX / factor}];
                acc.score += double(r.score);
                require(std::isfinite(acc.score), "derived score overflow");
            }
            for (const auto &e : sums)
                publish(Record{e.first.second, e.first.first, e.second.count,
                               static_cast<float>(e.second.score), z->type != 0});
        };
        visit(x0, x1, y0, y1);
        if (x0 == 0 && y0 == 0 && x1 == h.bins(key.first, unit, ri) && y1 == h.bins(key.second, unit, ri))
            require(emitted == z->occupied && (z->type || total == z->sum), "derived matrix statistics mismatch");
    }
    std::vector<double> vector(uint8_t kind, uint32_t norm, uint32_t chr, uint8_t unit, uint32_t ri,
                               double &scale, uint64_t begin = 0, uint64_t end = UINT64_MAX,
                               const std::function<void(uint32_t)> &consume = {},
                               std::vector<StrawNormalizationVector> *listing = nullptr) {
        Locator loc = kind == 0 ? h.norm : kind == 1 ? h.expected : h.normExpected;
        if (listing && !loc.len) return {};
        require(loc.len, "requested normalization/expected capability is absent");
        Bytes bytes;
        if (kind == 0 && normIndex.empty()) normIndex = source.read(loc.pos, loc.len);
        if (kind != 0) bytes = source.read(loc.pos, loc.len);
        Cursor c(kind == 0 ? normIndex : bytes);
        c.magic(kind == 0 ? "NVI0" : kind == 1 ? "EVI0" : "NEVI");
        require(c.word() == 1, "unknown vector index version");
        uint32_t n = c.word();
        c.zero(4);
        require(n <= c.left() / 40, "invalid vector entry count");
        std::array<uint32_t, 4> previous{};
        bool have = false, found = false;
        std::vector<double> result;
        for (uint32_t i = 0; i < n; ++i) {
            uint32_t len = c.word();
            require(len >= 40, "invalid vector entry size");
            Cursor e = c.take(len - 4);
            uint32_t ni = kind == 1 ? 0 : e.word(), ci = kind == 0 ? e.word() : 0;
            uint8_t u = e.byte();
            e.zero(3);
            uint32_t r = e.word(), bin = e.word();
            require(u <= 1 && r < h.resolutions[u].size() && bin == h.resolutions[u][r].bin &&
                        (kind == 1 || ni < h.norms.size()) && (kind != 0 || ci < h.chroms.size()),
                    "invalid vector key");
            std::array<uint32_t, 4> key{{ni, ci, u, r}};
            require(!have || previous < key, "unsorted/duplicate vector key");
            have = true;
            previous = key;
            uint64_t count = e.wide();
            uint32_t nominal = e.word(), chunks = e.word();
            uint64_t required = kind == 0 ? h.bins(ci, u, r) : 0;
            if (kind)
                for (uint32_t ch = 0; ch < h.chroms.size(); ++ch)
                    required = std::max(required, h.bins(ch, u, r));
            require(count == required && (!count || (nominal && chunks)) && (count || !chunks),
                    "invalid vector length/chunks");
            if (listing) {
                StrawNormalizationVector entry;
                entry.normalization = h.norms[ni]; entry.chromosome = h.chroms[ci].name;
                entry.unit = u ? "FRAG" : "BP"; entry.resolution = bin;
                entry.count = count; entry.binCount = required;
                listing->push_back(std::move(entry));
            }
            bool match = !listing &&
                u == unit && r == ri && (kind == 1 || ni == norm) && (kind != 0 || ci == chr);
            double factor = 1;
            if (kind) {
                uint32_t ns = e.word();
                e.zero(4);
                require(ns <= e.left() / 8, "invalid scale factor count");
                uint32_t prev = 0;
                for (uint32_t j = 0; j < ns; ++j) {
                    uint32_t ch = e.word(), bits = e.word();
                    require(ch < h.chroms.size() && (!j || ch > prev), "invalid scale factor key");
                    if (ch == chr)
                        factor = asFloat(bits);
                    prev = ch;
                }
            }
            require(uint64_t(chunks) * 32 == e.left(), "vector descriptor length mismatch");
            uint64_t next = 0;
            if (match) {
                end = std::min(end, count);
                require(begin <= end && (consume || end - begin <= allocationLimit / sizeof(double)),
                        "vector range exceeds allocation limit");
                if (!consume) result.reserve(end - begin);
                found = true;
                scale = factor;
            }
            for (uint32_t j = 0; j < chunks; ++j) {
                uint64_t first = e.wide();
                uint32_t nc = e.word();
                uint8_t transform = e.byte(), codec = e.byte();
                e.zero(2);
                uint64_t pos = e.wide();
                uint32_t stored = e.word(), raw = e.word();
                require(first == next && nc && uint64_t(nc) * 4 == raw && transform <= 2 &&
                            codec == 1 && stored > 16,
                        "invalid vector chunk descriptor");
                next = add(next, nc);
                require(next <= count, "vector chunk exceeds length");
                source.interval(pos, stored);
                if (!match || next <= begin || first >= end)
                    continue;
                auto chunkBytes = source.read(pos, stored);
                Cursor chunk(chunkBytes);
                chunk.magic("H10V");
                require(chunk.byte() == codec && chunk.byte() == transform,
                        "vector chunk codec/transform mismatch");
                chunk.zero(2);
                require(chunk.word() == raw && chunk.word() == nc, "vector chunk size mismatch");
                auto data = decompress(decoder.get(), chunk, raw);
                Cursor values(data);
                uint32_t prevBits = 0;
                for (uint32_t k = 0; k < nc; ++k) {
                    uint32_t bits = 0;
                    if (transform == 1)
                        for (unsigned lane = 0; lane < 4; ++lane)
                            bits |= uint32_t(data[uint64_t(lane) * nc + k]) << (8 * lane);
                    else {
                        bits = values.word();
                        if (transform == 2 && k)
                            bits ^= prevBits;
                    }
                    prevBits = bits;
                    if (first + k >= begin && first + k < end) {
                        if (consume) consume(bits);
                        else result.push_back(asFloat(bits));
                    }
                }
            }
            require(next == count, "incomplete vector coverage");
            e.done();
        }
        c.done();
        require(found || listing, "requested normalization/expected capability is absent");
        return result;
    }
    uint32_t normId(const std::string &norm) const {
        auto it = std::find(h.norms.begin(), h.norms.end(), norm);
        require(it != h.norms.end(), "unknown normalization " + norm);
        return static_cast<uint32_t>(it - h.norms.begin());
    }
    struct Region {
        uint32_t chr;
        uint64_t first, last, origin, end;
    };
    Region region(const std::string &location, uint8_t unit, uint32_t ri) const {
        auto colon = location.find(':');
        uint32_t chr = chromosomeId(location.substr(0, colon));
        uint64_t length = unit ? h.fragments[chr] : uint64_t(h.chroms[chr].length);
        uint64_t a = 0, b = length;
        if (colon != std::string::npos) {
            auto second = location.find(':', colon + 1);
            require(second != std::string::npos, "region needs start:end");
            auto number = [](const std::string &s) {
                require(!s.empty() && s.find_first_not_of("0123456789") == std::string::npos,
                        "invalid region coordinate");
                return std::stoull(s);
            };
            a = number(location.substr(colon + 1, second - colon - 1));
            b = number(location.substr(second + 1));
            // Straw's external genomic bounds are inclusive, for every version.
            // Internal bin rectangles remain half-open.
            require(a <= b && b <= length, "region out of bounds");
        }
        uint64_t bin = h.resolutions[unit][ri].bin;
        const uint64_t bins = h.bins(chr, unit, ri);
        const uint64_t first = std::min(bins, a / bin + (a % bin != 0));
        return {chr, first, std::max(first, std::min(bins, b / bin + 1)), a / bin, b / bin + 1};
    }
};

bool isV10(const std::string &path) {
    Source s(path);
    auto data = s.read(0, 8);
    Cursor c(data);
    c.magic("HIC\0");
    auto v = c.word();
    require(v >= 6 && v <= 10, "unsupported .hic version");
    return v == 10;
}
File::File(const std::string &path) : impl(new Impl(path)) {}
File::~File() = default;
std::vector<chromosome> File::chromosomes() const {
    return impl->h.chroms;
}
std::vector<int32_t> File::resolutions(const std::string &unit) const {
    std::vector<int32_t> out;
    for (auto r : impl->h.resolutions[unitId(unit)])
        out.push_back(static_cast<int32_t>(r.bin));
    return out;
}
std::vector<int32_t> File::derivedResolutions(const std::string &unit) const {
    std::vector<int32_t> out;
    for (auto r : impl->h.resolutions[unitId(unit)])
        if (r.mode) out.push_back(static_cast<int32_t>(r.bin));
    return out;
}
std::string File::genome() const { return impl->h.genome; }
std::vector<std::string> File::normalizations() const { return impl->h.norms; }
std::vector<std::pair<std::string, std::string>> File::attributes() const { return impl->h.attributes; }
bool File::allOverview(AllOverview &result) const {
    std::map<std::string, std::string> attrs;
    for (const auto &a : impl->h.attributes) attrs[a.first] = a.second;
    if (attrs["hicAll.version"] != "1" && attrs["hicAll.version"] != "2") return false;
    auto number = [&](const std::string &key) {
        const auto &text = attrs["hicAll." + key];
        require(!text.empty() && text.find_first_not_of("0123456789") == std::string::npos,
                "invalid scaled All " + key);
        uint64_t n = std::stoull(text);
        require(n > 0, "invalid scaled All " + key); return n;
    };
    result = {};
    result.version = u32(number("version"));
    result.finestResolution = u32(number("finestResolution"));
    result.sourceResolution = u32(number("sourceResolution"));
    result.groupFactor = number("groupFactor");
    result.binCount = u32(number("binCount"));
    require(result.groupFactor <= UINT64_MAX / result.sourceResolution, "All scale overflow");
    result.storageStride = result.version == 2 ? number("storageStride") : 1;
    result.coordinateScale = result.version == 2 ? number("coordinateScale") : 0;
    require(result.storageStride <= UINT64_MAX / result.finestResolution, "All stride overflow");
    result.displayResolution = result.storageStride * result.finestResolution;
    require(result.version != 2 || result.coordinateScale <= UINT64_MAX / result.displayResolution, "All physical scale overflow");
    result.effectiveResolution = result.version == 2 ? result.coordinateScale * result.displayResolution
        : result.groupFactor * result.sourceResolution;
    result.chromosomeIndex = UINT32_MAX;
    uint64_t total = 0;
    std::istringstream bins(attrs["hicAll.chromosomeBins"]);
    for (const auto &c : impl->h.chroms) {
        std::string name = c.name;
        std::transform(name.begin(), name.end(), name.begin(), [](unsigned char v) { return std::tolower(v); });
        if (name == "all") { result.chromosomeIndex = c.index; continue; }
        std::string word;
        require(bool(std::getline(bins, word, ',')) && !word.empty() &&
                word.find_first_not_of("0123456789") == std::string::npos, "invalid All chromosome layout");
        uint64_t n = std::stoull(word), length = c.length;
        require(n == length / result.sourceResolution + (length % result.sourceResolution != 0),
                "inconsistent All chromosome bins");
        result.chromosomeIndices.push_back(c.index);
        result.sourceOffsets.push_back(total); result.chromosomeBins.push_back(n);
        result.chromosomeLengths.push_back(length);
        uint64_t amount = result.version == 2 ? length : n;
        require(amount <= UINT64_MAX - total, "All layout overflow"); total += amount;
    }
    if (result.chromosomeIndex != UINT32_MAX) result.syntheticLength = impl->h.chroms[result.chromosomeIndex].length;
    bool consistent = result.version == 2
        ? result.syntheticLength == std::max<uint64_t>(1, total / result.coordinateScale) &&
          result.binCount == result.syntheticLength / result.displayResolution + (result.syntheticLength % result.displayResolution != 0)
        : result.syntheticLength == uint64_t(result.binCount) * result.finestResolution &&
          total / result.groupFactor + (total % result.groupFactor != 0) == result.binCount;
    std::string trailing;
    require(!std::getline(bins, trailing, ',') && result.chromosomeIndex != UINT32_MAX &&
            !impl->h.resolutions[0].empty() && result.finestResolution == impl->h.resolutions[0][0].bin &&
            consistent,
            "inconsistent scaled All overview");
    return true;
}
bool File::hasMatrix(const std::string &chr1, const std::string &chr2,
                     const std::string &unit, int32_t resolution) {
    auto &p = *impl;
    uint32_t a = p.chromosomeId(chr1), b = p.chromosomeId(chr2);
    if (a > b) std::swap(a, b);
    uint8_t u = unitId(unit);
    uint32_t ri = p.resolutionId(u, resolution);
    return p.zoom({a, b}, u, ri) != nullptr;
}
void File::raw(const std::string &chr1, const std::string &chr2, const std::string &unit,
               int32_t resolution, uint64_t x0, uint64_t x1, uint64_t y0, uint64_t y1,
               const Callback &cb) {
    auto &p = *impl;
    uint32_t a = p.chromosomeId(chr1), b = p.chromosomeId(chr2);
    uint8_t u = unitId(unit);
    uint32_t ri = p.resolutionId(u, resolution);
    require(x0 <= x1 && y0 <= y1, "reversed bin interval");
    x1 = std::min(x1, p.h.bins(a, u, ri));
    y1 = std::min(y1, p.h.bins(b, u, ri));
    if (x0 >= x1 || y0 >= y1)
        return;
    bool transpose = a > b;
    if (transpose) {
        std::swap(a, b);
        std::swap(x0, y0);
        std::swap(x1, y1);
    }
    p.raw({a, b}, u, ri, x0, x1, y0, y1, [&](Record r) {
        if (transpose)
            std::swap(r.binX, r.binY);
        cb(r);
    });
}
void File::streamRaw(const std::string &a, const std::string &b, const std::string &unit,
                     int32_t resolution, const Callback &cb) {
    auto &p = *impl;
    auto u = unitId(unit);
    auto ri = p.resolutionId(u, resolution);
    auto x = p.region(a, u, ri), y = p.region(b, u, ri);
    raw(p.h.chroms[x.chr].name, p.h.chroms[y.chr].name, unit, resolution, x.first, x.last, y.first,
        y.last, [&](Record r) {
            if (x.chr == y.chr &&
                !(r.binX >= x.first && r.binX < x.last && r.binY >= y.first && r.binY < y.last))
                std::swap(r.binX, r.binY);
            cb(r);
        });
}
std::vector<double> File::normalization(const std::string &chr, const std::string &unit,
                                        int32_t resolution, const std::string &norm) {
    auto &p = *impl;
    auto ch = p.chromosomeId(chr);
    auto u = unitId(unit);
    auto ri = p.resolutionId(u, resolution);
    if (norm == "NONE") {
        require(p.h.bins(ch, u, ri) <= allocationLimit / sizeof(double),
                "vector exceeds allocation limit");
        return std::vector<double>(p.h.bins(ch, u, ri), 1);
    }
    double scale = 1;
    return p.vector(0, p.normId(norm), ch, u, ri, scale);
}
std::vector<StrawNormalizationVector> File::normalizationEntries() {
    std::vector<StrawNormalizationVector> result;
    double scale = 1;
    impl->vector(0, 0, 0, 0, 0, scale, 0, UINT64_MAX, {}, &result);
    return result;
}
void File::streamNormalization(const StrawNormalizationVector &entry,
                               const std::function<void(uint32_t)> &consume) {
    require(bool(consume), "normalization consumer is required");
    auto &p = *impl;
    auto u = unitId(entry.unit);
    double scale = 1;
    p.vector(0, p.normId(entry.normalization), p.chromosomeId(entry.chromosome), u,
             p.resolutionId(u, entry.resolution), scale, 0, UINT64_MAX, consume);
}
std::vector<double> File::expected(const std::string &chr, const std::string &unit,
                                   int32_t resolution, const std::string &norm) {
    auto &p = *impl;
    auto ch = p.chromosomeId(chr);
    auto u = unitId(unit);
    auto ri = p.resolutionId(u, resolution);
    double scale = 1;
    auto values = p.vector(norm == "NONE" ? 1 : 2,
                           norm == "NONE" ? 0 : p.normId(norm), ch, u, ri, scale);
    require(scale && std::isfinite(scale), "invalid expected-vector scale");
    for (auto &value : values) value /= scale;
    return values;
}
void File::stream(const std::string &matrixType, const std::string &norm,
                  const std::string &chr1loc, const std::string &chr2loc, const std::string &unit,
                  int32_t resolution, const StrawRecordCallback &cb) {
    require(matrixType == "observed" || matrixType == "oe" || matrixType == "expected",
            "unknown matrix type");
    auto &p = *impl;
    uint8_t u = unitId(unit);
    uint32_t ri = p.resolutionId(u, resolution);
    auto x = p.region(chr1loc, u, ri), y = p.region(chr2loc, u, ri);
    if (x.first == x.last || y.first == y.last)
        return;
    std::vector<double> n1, n2, expected;
    double scale = 1;
    double transExpected = 0;
    uint64_t expectedBegin = 0;
    if (norm != "NONE" && matrixType != "expected") {
        auto ni = p.normId(norm);
        n1 = p.vector(0, ni, x.chr, u, ri, scale, x.first, x.last);
        n2 = p.vector(0, ni, y.chr, u, ri, scale, y.first, y.last);
    }
    if (matrixType != "observed" && x.chr != y.chr) {
        const auto *z = p.zoom({std::min(x.chr, y.chr), std::max(x.chr, y.chr)}, u, ri);
        if (!z) return;
        double sum = double(z->sum);
        if (z->type) std::memcpy(&sum, &z->sum, sizeof(sum));
        // As in V9, use the raw whole-matrix mean even for normalized O/E.
        transExpected = (sum / p.h.bins(x.chr, u, ri)) / p.h.bins(y.chr, u, ri);
    }
    if (matrixType != "observed" && x.chr == y.chr) {
        expectedBegin = x.last <= y.first   ? y.first - x.last + 1
                        : y.last <= x.first ? x.first - y.last + 1
                                            : 0;
        uint64_t expectedEnd = std::max(x.last > y.first ? x.last - y.first : y.first - x.last + 2,
                                        y.last > x.first ? y.last - x.first : x.first - y.last + 2);
        expected = p.vector(norm == "NONE" ? 1 : 2, norm == "NONE" ? 0 : p.normId(norm), x.chr, u,
                            ri, scale, expectedBegin, expectedEnd);
    }
    raw(p.h.chroms[x.chr].name, p.h.chroms[y.chr].name, unit, resolution, x.first, x.last, y.first,
        y.last, [&](Record r) {
            if (x.chr == y.chr &&
                !(r.binX >= x.first && r.binX < x.last && r.binY >= y.first && r.binY < y.last))
                std::swap(r.binX, r.binY);
            double value = r.isScore ? double(r.score) : double(r.count);
            if (norm != "NONE" && matrixType != "expected") {
                double a = n1.at(r.binX - x.first), b = n2.at(r.binY - y.first);
                if (!a || !b || !std::isfinite(a) || !std::isfinite(b))
                    return;
                value /= a * b;
            }
            if (matrixType != "observed") {
                double e = transExpected;
                if (x.chr == y.chr) {
                    uint64_t d = r.binX > r.binY ? r.binX - r.binY : r.binY - r.binX;
                    if (d < expectedBegin || d - expectedBegin >= expected.size() || !scale ||
                        !std::isfinite(scale))
                        return;
                    e = expected[d - expectedBegin] / scale;
                }
                if (!e || !std::isfinite(e))
                    return;
                value = matrixType == "oe" ? value / e : e;
            }
            uint64_t bx = uint64_t(r.binX) * resolution, by = uint64_t(r.binY) * resolution;
            require(bx <= INT32_MAX && by <= INT32_MAX,
                    "legacy contactRecord coordinate overflow; use File::raw");
            cb({static_cast<int32_t>(bx), static_cast<int32_t>(by), static_cast<float>(value)});
        });
}
std::vector<std::vector<float>> File::matrix(const std::string &type, const std::string &norm,
                                             const std::string &a, const std::string &b,
                                             const std::string &unit, int32_t resolution) {
    auto u = unitId(unit);
    auto ri = impl->resolutionId(u, resolution);
    auto x = impl->region(a, u, ri), y = impl->region(b, u, ri);
    uint64_t rows = x.end - x.origin, cols = y.end - y.origin;
    require(rows <= allocationLimit / sizeof(std::vector<float>) &&
                (!cols || rows <= allocationLimit / sizeof(float) / cols),
            "dense matrix exceeds allocation limit");
    std::vector<std::vector<float>> result(rows, std::vector<float>(cols, 0));
    stream(type, norm, a, b, unit, resolution, [&](const contactRecord &r) {
        uint64_t bx = r.binX / resolution, by = r.binY / resolution;
        result.at(bx - x.origin).at(by - y.origin) = r.counts;
        if (x.chr == y.chr && bx != by && by >= x.first && by < x.last && bx >= y.first &&
            bx < y.last)
            result[by - x.origin][bx - y.origin] = r.counts;
    });
    return result;
}
int64_t File::countRecords(int32_t resolution, bool interOnly, bool printByChromosome) {
    uint32_t ri = impl->resolutionId(0, resolution);
    uint64_t total = 0;
    for (uint32_t a = 0; a < impl->h.chroms.size(); ++a) {
        if (impl->h.chroms[a].name == "All" || impl->h.chroms[a].name == "ALL")
            continue;
        for (uint32_t b = a; b < impl->h.chroms.size(); ++b) {
            if ((interOnly && a == b) || (printByChromosome && a != b))
                continue;
            const Zoom *z = impl->zoom({a, b}, 0, ri);
            uint64_t n = z ? z->occupied : 0;
            total = add(total, n);
            if (printByChromosome)
                std::cout << impl->h.chroms[a].name << " " << n << " " << (n * 12 / 1000000000)
                          << " GB\n";
        }
    }
    require(total <= INT64_MAX, "record count exceeds legacy API");
    return printByChromosome ? 0 : static_cast<int64_t>(total);
}

std::vector<std::pair<std::string, int64_t>> File::countRecordsByChromosome(int32_t resolution) {
    uint32_t ri = impl->resolutionId(0, resolution);
    std::vector<std::pair<std::string, int64_t>> result;
    for (uint32_t a = 0; a < impl->h.chroms.size(); ++a) {
        if (impl->h.chroms[a].name == "All" || impl->h.chroms[a].name == "ALL")
            continue;
        const Zoom *z = impl->zoom({a, a}, 0, ri);
        uint64_t n = z ? z->occupied : 0;
        require(n <= static_cast<uint64_t>(INT64_MAX), "record count exceeds legacy API");
        result.emplace_back(impl->h.chroms[a].name, static_cast<int64_t>(n));
    }
    return result;
}
void readHeader(std::istream &input, int64_t &master, std::string &genome, int32_t &count,
                int64_t &nvi, int64_t &nviLength, std::map<std::string, chromosome> &chromosomes) {
    input.clear();
    input.seekg(0);
    Bytes prefix(16);
    input.read(reinterpret_cast<char *>(prefix.data()), 16);
    require(bool(input), "short header");
    Cursor c(prefix);
    c.magic("HIC\0");
    require(c.word() == 10, "not V10");
    auto n = c.wide();
    require(n >= 88 && n <= allocationLimit, "invalid header length");
    Bytes bytes(n);
    input.seekg(0);
    input.read(reinterpret_cast<char *>(bytes.data()), n);
    require(bool(input), "short header");
    auto h = parseHeader(bytes);
    require(h.footer.pos <= INT64_MAX && h.norm.pos <= INT64_MAX && h.norm.len <= INT64_MAX,
            "header offset exceeds legacy API");
    master = h.footer.pos;
    genome = h.genome;
    count = static_cast<int32_t>(h.chroms.size());
    nvi = h.norm.pos;
    nviLength = h.norm.len;
    chromosomes.clear();
    for (auto chr : h.chroms)
        chromosomes.emplace(chr.name, chr);
}
} // namespace straw_v10

#include "hic_slice.h"
#include <zlib.h>
// Keep the legacy slice wire format while applying the requested V10 query.
void dumpV10(const std::string &type, const std::string &norm, const std::string &path,
             const std::string &unit, int32_t resolution, const std::string &output,
             bool compressed, ContactFilter filter) {
    straw_v10::File file(path);
    FILE *plain = nullptr;
    gzFile gzip = nullptr;
    if (compressed)
        gzip = gzopen(output.c_str(), "wb");
    else
        plain = std::fopen(output.c_str(), "wb");
    straw_v10::require(compressed ? gzip != nullptr : plain != nullptr, "cannot open slice output");
    auto write = [&](const void *data, size_t size) {
        bool ok = compressed
                      ? gzwrite(gzip, data, static_cast<unsigned>(size)) == static_cast<int>(size)
                      : std::fwrite(data, 1, size, plain) == size;
        straw_v10::require(ok, "slice write failure");
    };
    try {
        auto chromosomes = file.chromosomes();
        std::map<std::string, int16_t> keys;
        for (const auto &c : chromosomes)
            if (c.name != "All" && c.name != "ALL") {
                straw_v10::require(keys.size() < INT16_MAX, "too many slice chromosomes");
                keys[c.name] = static_cast<int16_t>(keys.size());
            }
        int32_t n = static_cast<int32_t>(keys.size());
        write("HICSLICE", 8);
        write(&resolution, 4);
        write(&n, 4);
        for (const auto &k : keys) {
            int32_t size = static_cast<int32_t>(k.first.size());
            write(&size, 4);
            write(k.first.data(), size);
            write(&k.second, 2);
        }
        for (const auto &a : chromosomes)
            for (const auto &b : chromosomes) {
                if (!keys.count(a.name) || !keys.count(b.name) || a.index > b.index)
                    continue;
                if (filter == ContactFilter::INTER && a.index == b.index)
                    continue;
                if (filter != ContactFilter::ALL && filter != ContactFilter::INTER &&
                    a.index != b.index)
                    continue;
                file.stream(type, norm, a.name, b.name, unit, resolution,
                            [&](const contactRecord &r) {
                                int64_t d = std::abs(int64_t(r.binY) - r.binX);
                                if (filter == ContactFilter::INTRA_SHORT && d >= 5000000)
                                    return;
                                if (filter == ContactFilter::INTRA_LONG && d <= 5000000)
                                    return;
                                CompressedContactRecord record{};
                                record.chr1Key = keys.at(a.name);
                                record.chr2Key = keys.at(b.name);
                                record.binX = r.binX / resolution;
                                record.binY = r.binY / resolution;
                                record.value = r.counts;
                                write(&record, sizeof(record));
                            });
            }
    } catch (...) {
        if (gzip)
            gzclose(gzip);
        if (plain)
            std::fclose(plain);
        throw;
    }
    int status = compressed ? gzclose(gzip) : std::fclose(plain);
    straw_v10::require(status == 0, "cannot finish slice output");
}
