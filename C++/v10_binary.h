#pragma once
// Bounds-checked primitives used only by the V10 reader.
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace straw_v10 {
using Bytes = std::vector<uint8_t>;
inline void require(bool ok, const std::string &message) {
    if (!ok)
        throw std::runtime_error("V10: " + message);
}
inline uint64_t add(uint64_t a, uint64_t b) {
    require(b <= UINT64_MAX - a, "integer overflow");
    return a + b;
}
inline uint32_t u32(uint64_t x) {
    require(x <= UINT32_MAX, "value exceeds uint32");
    return static_cast<uint32_t>(x);
}
inline float asFloat(uint32_t bits) {
    float f;
    std::memcpy(&f, &bits, 4);
    return f;
}
inline double asDouble(uint64_t bits) {
    double f;
    std::memcpy(&f, &bits, 8);
    return f;
}
// Per-record allocation ceiling; independent of the size of the .hic file.
constexpr uint64_t allocationLimit = 512ULL * 1024 * 1024;
struct Cursor {
    const uint8_t *p;
    size_t size, at = 0;
    explicit Cursor(const Bytes &b) : p(b.data()), size(b.size()) {}
    Cursor(const uint8_t *b, size_t n) : p(b), size(n) {}
    size_t left() const {
        return size - at;
    }
    void need(uint64_t n) const {
        require(n <= left(), "truncated record");
    }
    uint64_t integer(unsigned n) {
        need(n);
        uint64_t v = 0;
        for (unsigned i = 0; i < n; ++i)
            v |= uint64_t(p[at++]) << (8 * i);
        return v;
    }
    uint8_t byte() {
        return static_cast<uint8_t>(integer(1));
    }
    uint32_t word() {
        return static_cast<uint32_t>(integer(4));
    }
    uint64_t wide() {
        return integer(8);
    }
    uint64_t var() {
        uint64_t v = 0;
        for (unsigned i = 0; i < 10; ++i) {
            uint8_t b = byte();
            require(i < 9 || b <= 1, "overflowing ULEB128");
            v |= uint64_t(b & 127) << (i * 7);
            if (!(b & 128)) {
                require(i == 0 || b != 0, "non-canonical ULEB128");
                return v;
            }
        }
        throw std::runtime_error("V10: unterminated ULEB128");
    }
    void zero(size_t n) {
        while (n--)
            require(byte() == 0, "nonzero reserved field");
    }
    void magic(const char *s) {
        need(4);
        require(std::memcmp(p + at, s, 4) == 0, "bad record magic");
        at += 4;
    }
    std::string str() {
        size_t begin = at;
        while (at < size && p[at])
            ++at;
        require(at < size && at - begin <= 1024 * 1024, "invalid string length or terminator");
        std::string s(reinterpret_cast<const char *>(p + begin), at - begin);
        ++at;
        // Validate UTF-8, including overlong forms, surrogates and > U+10ffff.
        for (size_t i = 0; i < s.size();) {
            uint8_t c = static_cast<uint8_t>(s[i++]);
            if (c < 128)
                continue;
            unsigned n = c >= 0xc2 && c <= 0xdf   ? 1
                         : c >= 0xe0 && c <= 0xef ? 2
                         : c >= 0xf0 && c <= 0xf4 ? 3
                                                  : 0;
            require(n && i + n <= s.size(), "invalid UTF-8");
            uint32_t v = c & ((1u << (6 - n)) - 1);
            for (unsigned j = 0; j < n; ++j) {
                uint8_t d = static_cast<uint8_t>(s[i++]);
                require((d & 0xc0) == 0x80, "invalid UTF-8");
                v = (v << 6) | (d & 63);
            }
            require(v >= (n == 1   ? 128u
                          : n == 2 ? 2048u
                                   : 65536u) &&
                        v <= 0x10ffff && !(v >= 0xd800 && v <= 0xdfff),
                    "invalid UTF-8");
        }
        return s;
    }
    Cursor take(uint64_t n) {
        need(n);
        Cursor c(p + at, static_cast<size_t>(n));
        at += n;
        return c;
    }
    void done() const {
        require(at == size, "trailing bytes in record");
    }
};
} // namespace straw_v10
