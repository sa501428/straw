#include "straw_v10.h"
#include "v10_binary.h"
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
struct Source {
    std::ifstream file;
    CURL *curl = nullptr;
    uint64_t size = 0, requested = 0, responseStart = UINT64_MAX, responseEnd = 0;
    Bytes response;
    static size_t write(char *p, size_t a, size_t b, void *context) {
        auto &s = *static_cast<Source *>(context);
        if (b && a > SIZE_MAX / b)
            return 0;
        size_t n = a * b;
        if (n > s.requested - s.response.size())
            return 0;
        try {
            s.response.insert(s.response.end(), p, p + n);
        } catch (...) {
            return 0;
        }
        return n;
    }
    static size_t header(char *p, size_t a, size_t b, void *context) {
        auto &s = *static_cast<Source *>(context);
        size_t n = a * b;
        try {
            std::string line(p, n), lower = line;
            std::transform(lower.begin(), lower.end(), lower.begin(),
                           [](unsigned char c) { return std::tolower(c); });
            if (lower.compare(0, 14, "content-range:") == 0) {
                unsigned long long first, last, total;
                if (std::sscanf(lower.c_str(), "content-range: bytes %llu-%llu/%llu", &first, &last,
                                &total) == 3) {
                    s.responseStart = first;
                    s.responseEnd = last;
                    s.size = total;
                }
            }
        } catch (...) {
            return 0;
        }
        return n;
    }
    explicit Source(const std::string &path) {
        if (path.compare(0, 7, "http://") == 0 || path.compare(0, 8, "https://") == 0) {
            static std::once_flag once;
            std::call_once(once, [] {
                require(curl_global_init(CURL_GLOBAL_DEFAULT) == CURLE_OK,
                        "curl initialization failed");
            });
            curl = curl_easy_init();
            require(curl != nullptr, "curl initialization failed");
            curl_easy_setopt(curl, CURLOPT_URL, path.c_str());
            curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
            curl_easy_setopt(curl, CURLOPT_MAXREDIRS, 10L);
            curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L);
            curl_easy_setopt(curl, CURLOPT_CONNECTTIMEOUT, 30L);
            curl_easy_setopt(curl, CURLOPT_LOW_SPEED_LIMIT, 1L);
            curl_easy_setopt(curl, CURLOPT_LOW_SPEED_TIME, 60L);
            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write);
            curl_easy_setopt(curl, CURLOPT_WRITEDATA, this);
            curl_easy_setopt(curl, CURLOPT_HEADERFUNCTION, header);
            curl_easy_setopt(curl, CURLOPT_HEADERDATA, this);
            curl_easy_setopt(curl, CURLOPT_USERAGENT, "straw-v10");
        } else {
            file.open(path, std::ios::binary);
            require(bool(file), "cannot open " + path);
            file.seekg(0, std::ios::end);
            auto end = file.tellg();
            require(end >= 0, "cannot stat input");
            size = uint64_t(end);
        }
    }
    ~Source() {
        if (curl)
            curl_easy_cleanup(curl);
    }
    void interval(uint64_t pos, uint64_t len) const {
        require(len > 0 && pos <= size && len <= size - pos, "file interval out of bounds");
    }
    Bytes read(uint64_t pos, uint64_t len) {
        require(len > 0 && len <= allocationLimit, "record exceeds allocation limit");
        if (size)
            interval(pos, len);
        if (!curl) {
            Bytes out(static_cast<size_t>(len));
            file.clear();
            file.seekg(static_cast<std::streamoff>(pos));
            file.read(reinterpret_cast<char *>(out.data()), static_cast<std::streamsize>(len));
            require(bool(file), "short file read");
            return out;
        }
        requested = len;
        response.clear();
        responseStart = UINT64_MAX;
        std::string range = std::to_string(pos) + "-" + std::to_string(add(pos, len) - 1);
        curl_easy_setopt(curl, CURLOPT_RANGE, range.c_str());
        CURLcode rc = curl_easy_perform(curl);
        long status = 0;
        curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &status);
        require(rc == CURLE_OK, std::string("HTTP range read failed: ") + curl_easy_strerror(rc));
        require(status == 206 && responseStart == pos && responseEnd == pos + len - 1 &&
                    response.size() == len,
                "server must return the exact requested HTTP byte range (206)");
        interval(pos, len);
        return std::move(response);
    }
};
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
    auto requiredSource = [](uint32_t bin) -> uint32_t {
        if (bin == 20 || bin == 50) return 10;
        if (bin == 200 || bin == 500) return 100;
        if (bin == 2000) return 1000;
        return 0;
    };
    for (uint32_t i = 0; i < h.resolutions[0].size(); ++i) {
        const auto &r = h.resolutions[0][i];
        if (uint32_t sourceBin = requiredSource(r.bin)) {
            auto source = std::find_if(h.resolutions[0].begin(), h.resolutions[0].end(),
                                       [&](const Resolution &s) { return s.bin == sourceBin; });
            require(source != h.resolutions[0].end() && r.mode &&
                        r.source == uint32_t(source - h.resolutions[0].begin()) && !source->mode,
                    "mandatory BP derivation policy is not satisfied");
        }
        require(r.bin != 500000 || !r.mode, "500 kb must be materialized");
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
        const Zoom &s = *zoom(key, unit, z->source);
        uint64_t factor = z->bin / s.bin;
        struct Acc {
            uint64_t count = 0;
            double score = 0;
        };
        std::map<std::pair<uint32_t, uint32_t>, Acc> sums;
        // The source cells used to be buffered whole in a second std::map purely
        // so duplicates could be detected. A hash set of the packed coordinate
        // does the same job without holding a Record per source cell.
        std::unordered_set<uint64_t> seen;
        // Integer counts are order-independent, so they accumulate as cells
        // arrive. Float scores are not, so they keep the previous behaviour of
        // being summed in ascending (binY, binX) source order.
        std::vector<std::pair<std::pair<uint32_t, uint32_t>, float>> scores;
        materialized(key, s, x0 * factor, std::min(x1 * factor, h.bins(key.first, unit, s.ri)),
                     y0 * factor, std::min(y1 * factor, h.bins(key.second, unit, s.ri)),
                     [&](const Record &r) {
                         uint32_t x = r.binX / factor, y = r.binY / factor;
                         if (!inside(x, y))
                             return;
                         require(seen.insert((uint64_t(r.binY) << 32) | r.binX).second,
                                 "duplicate source cell");
                         if (!z->type) {
                             auto &a = sums[{y, x}];
                             a.count = add(a.count, r.count);
                         } else {
                             scores.emplace_back(std::make_pair(r.binY, r.binX), r.score);
                         }
                     });
        if (z->type) {
            std::sort(scores.begin(), scores.end(),
                      [](const std::pair<std::pair<uint32_t, uint32_t>, float> &a,
                         const std::pair<std::pair<uint32_t, uint32_t>, float> &b) {
                          return a.first < b.first;
                      });
            for (const auto &entry : scores) {
                auto &a = sums[{entry.first.first / factor, entry.first.second / factor}];
                require(std::isfinite(entry.second), "nonfinite derived source score");
                a.score += double(entry.second);
                require(std::isfinite(a.score), "derived score overflow");
            }
        }
        if (x0 == 0 && y0 == 0 && x1 == h.bins(key.first, unit, ri) &&
            y1 == h.bins(key.second, unit, ri)) {
            require(sums.size() == z->occupied, "derived occupied count mismatch");
            if (!z->type) {
                uint64_t total = 0;
                for (const auto &e : sums)
                    total = add(total, e.second.count);
                require(total == z->sum, "derived sum mismatch");
            }
        }
        for (const auto &e : sums)
            cb(Record{e.first.second, e.first.first, e.second.count,
                      static_cast<float>(e.second.score), z->type != 0});
    }
    std::vector<double> vector(uint8_t kind, uint32_t norm, uint32_t chr, uint8_t unit, uint32_t ri,
                               double &scale, uint64_t begin = 0, uint64_t end = UINT64_MAX) {
        Locator loc = kind == 0 ? h.norm : kind == 1 ? h.expected : h.normExpected;
        require(loc.len, "requested normalization/expected capability is absent");
        auto bytes = source.read(loc.pos, loc.len);
        Cursor c(bytes);
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
            bool match =
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
                require(begin <= end && end - begin <= allocationLimit / sizeof(double),
                        "vector range exceeds allocation limit");
                result.reserve(end - begin);
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
                    if (first + k >= begin && first + k < end)
                        result.push_back(asFloat(bits));
                }
            }
            require(next == count, "incomplete vector coverage");
            e.done();
        }
        c.done();
        require(found, "requested normalization/expected capability is absent");
        return result;
    }
    uint32_t normId(const std::string &norm) const {
        auto it = std::find(h.norms.begin(), h.norms.end(), norm);
        require(it != h.norms.end(), "unknown normalization " + norm);
        return static_cast<uint32_t>(it - h.norms.begin());
    }
    struct Region {
        uint32_t chr;
        uint64_t first, last;
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
            // Consume the caller's numeric values directly. V10 does not
            // classify or convert their coordinate origin; this API's end is
            // exclusive and may equal the chromosome length.
            require(a <= b && b <= length, "region out of bounds");
        }
        uint64_t bin = h.resolutions[unit][ri].bin;
        return {chr, a / bin, a == b ? a / bin : b / bin + (b % bin != 0)};
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
    uint64_t expectedBegin = 0;
    if (norm != "NONE" && matrixType != "expected") {
        auto ni = p.normId(norm);
        n1 = p.vector(0, ni, x.chr, u, ri, scale, x.first, x.last);
        n2 = p.vector(0, ni, y.chr, u, ri, scale, y.first, y.last);
    }
    if (matrixType != "observed") {
        require(x.chr == y.chr, "V10 expected values are defined only for cis matrices");
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
                uint64_t d = r.binX > r.binY ? r.binX - r.binY : r.binY - r.binX;
                if (d < expectedBegin || d - expectedBegin >= expected.size() || !scale ||
                    !std::isfinite(scale))
                    return;
                double e = expected[d - expectedBegin] / scale;
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
    uint64_t rows = x.last - x.first, cols = y.last - y.first;
    require(rows <= allocationLimit / sizeof(std::vector<float>) &&
                (!cols || rows <= allocationLimit / sizeof(float) / cols),
            "dense matrix exceeds allocation limit");
    std::vector<std::vector<float>> result(rows, std::vector<float>(cols, 0));
    stream(type, norm, a, b, unit, resolution, [&](const contactRecord &r) {
        uint64_t bx = r.binX / resolution, by = r.binY / resolution;
        result.at(bx - x.first).at(by - y.first) = r.counts;
        if (x.chr == y.chr && bx != by && by >= x.first && by < x.last && bx >= y.first &&
            bx < y.last)
            result[by - x.first][bx - y.first] = r.counts;
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
