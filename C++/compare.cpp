#include "compare.h"
#include "straw.h"
#include "straw_v10.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Options {
    std::string first, second;
    int32_t exhaustiveAt = 100000;
    size_t samples = 4;
    uint64_t windowBins = 256;
    uint64_t seed = 1;
    double absTol = 0.0;
    double relTol = 1e-6;
    size_t maxErrors = 20;
    bool exhaustiveAll = false;
    std::vector<std::string> norms{"VC", "VC_SQRT", "KR"};
};

struct Summary {
    uint64_t metadata = 0, matrices = 0, windows = 0, cells = 0;
    uint64_t vectors = 0, vectorValues = 0, differences = 0;
};

using Cell = std::pair<uint64_t, uint64_t>;
using Cells = std::map<Cell, double>;

bool closeEnough(double a, double b, double absTol, double relTol) {
    if (std::isnan(a) || std::isnan(b)) return std::isnan(a) && std::isnan(b);
    if (std::isinf(a) || std::isinf(b)) return a == b;
    double d = std::fabs(a - b);
    return d <= absTol || d <= relTol * std::max(std::fabs(a), std::fabs(b));
}

void difference(Summary &s, const Options &o, const std::string &message) {
    ++s.differences;
    if (s.differences <= o.maxErrors) std::cerr << "DIFF: " << message << '\n';
    else if (s.differences == o.maxErrors + 1)
        std::cerr << "DIFF: further differences suppressed\n";
}

std::string location(const std::string &chr, uint64_t begin, uint64_t end) {
    return chr + ":" + std::to_string(begin) + ":" + std::to_string(end);
}

Cells rawCells(const std::string &path, const std::string &chr1, const std::string &chr2,
               int32_t resolution, uint64_t x0, uint64_t x1, uint64_t y0, uint64_t y1,
               bool v10, bool transposeLegacy) {
    Cells out;
    auto add = [&](uint64_t x, uint64_t y, double value) {
        // Legacy queries can include a bin touching the inclusive region end.
        // Filtering here gives both readers identical half-open semantics.
        bool direct = x >= x0 && x < x1 && y >= y0 && y < y1;
        bool reflected = chr1 == chr2 && y >= x0 && y < x1 && x >= y0 && x < y1;
        if (!direct && !reflected) return;
        if (!direct) std::swap(x, y);
        out[{x / static_cast<uint64_t>(resolution), y / static_cast<uint64_t>(resolution)}] += value;
    };
    const std::string a = location(chr1, x0, x1), b = location(chr2, y0, y1);
    if (v10) {
        straw_v10::File(path).streamRaw(a, b, "BP", resolution,
            [&](const straw_v10::Record &r) {
                double value = r.isScore ? static_cast<double>(r.score) : static_cast<double>(r.count);
                add(uint64_t(r.binX) * resolution, uint64_t(r.binY) * resolution, value);
            });
    } else {
        strawStream("observed", "NONE", path, a, b, "BP", resolution,
            [&](const contactRecord &r) {
                uint64_t x = r.binX, y = r.binY;
                // The V6-V9 reader emits trans records in file chromosome-index
                // order even when the caller requested the opposite order. V10
                // emits them in request order. Put legacy records back into the
                // requested orientation before filtering and comparing them.
                if (transposeLegacy) std::swap(x, y);
                add(x, y, r.counts);
            });
    }
    return out;
}

void compareCells(const Cells &a, const Cells &b, const Options &o, Summary &s,
                  const std::string &context) {
    auto i = a.begin(), j = b.begin();
    while (i != a.end() || j != b.end()) {
        Cell key;
        double av = 0, bv = 0;
        if (j == b.end() || (i != a.end() && i->first < j->first)) {
            key = i->first; av = i->second; ++i;
        } else if (i == a.end() || j->first < i->first) {
            key = j->first; bv = j->second; ++j;
        } else {
            key = i->first; av = i->second; bv = j->second; ++i; ++j;
        }
        ++s.cells;
        if (!closeEnough(av, bv, o.absTol, o.relTol)) {
            std::ostringstream m;
            m << context << " bin(" << key.first << ',' << key.second << ") "
              << std::setprecision(17) << av << " != " << bv;
            difference(s, o, m.str());
        }
    }
}

void compareVector(const std::string &kind, const std::string &chr, int32_t resolution,
                   const std::string &norm, const std::vector<double> &a,
                   const std::vector<double> &b, bool exhaustive, bool firstV10, bool secondV10,
                   size_t canonicalLength, const Options &o, Summary &s) {
    ++s.vectors;
    auto compatibleLengths = [&]() {
        if (a.size() == b.size()) return true;
        if (firstV10 == secondV10) return false;
        size_t v10 = firstV10 ? a.size() : b.size();
        size_t legacy = firstV10 ? b.size() : a.size();
        if (v10 != canonicalLength) return false;
        // V9 normalization vectors commonly contain floor(length/bin)+1
        // entries, while V10 contains exactly ceil(length/bin). The only
        // difference is an unused trailing V9 entry when the length divides
        // the bin size exactly.
        if (kind == "normalization")
            return legacy == canonicalLength + 1;
        // V9 expected vectors commonly contain floor(maxLength/bin) entries;
        // V10 covers every possible distance with ceil(maxLength/bin).
        return canonicalLength && legacy + 1 == canonicalLength;
    };
    if (!compatibleLengths()) {
        difference(s, o, kind + " " + norm + " " + chr + " @" + std::to_string(resolution) +
                          " length " + std::to_string(a.size()) + " != " + std::to_string(b.size()));
    }
    size_t n = std::min(a.size(), b.size());
    std::vector<size_t> indices;
    if (exhaustive || n <= o.windowBins * std::max<size_t>(1, o.samples)) {
        indices.resize(n);
        for (size_t i = 0; i < n; ++i) indices[i] = i;
    } else if (n) {
        std::mt19937_64 rng(o.seed ^ uint64_t(resolution) ^ std::hash<std::string>{}(chr + kind + norm));
        indices = {0, n - 1, n / 2};
        for (size_t k = 0; k < o.samples * o.windowBins; ++k) indices.push_back(rng() % n);
        std::sort(indices.begin(), indices.end());
        indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
    }
    for (size_t i : indices) {
        ++s.vectorValues;
        if (!closeEnough(a[i], b[i], o.absTol, o.relTol)) {
            std::ostringstream m;
            m << kind << ' ' << norm << ' ' << chr << " @" << resolution << " index " << i
              << ' ' << std::setprecision(17) << a[i] << " != " << b[i];
            difference(s, o, m.str());
        }
    }
}

std::map<std::string, chromosome> chromosomeMap(const std::string &path) {
    std::map<std::string, chromosome> out;
    for (const auto &c : getChromosomesForFile(path))
        if (c.name != "ALL" && c.name != "All" && c.name != "all") out[c.name] = c;
    return out;
}

Options parse(int argc, char **argv) {
    if (argc < 4) throw std::runtime_error("usage: straw compare <first.hic> <second.hic> [--all] [--exhaustive-at BP] [--samples N] [--window-bins N] [--seed N] [--abs-tol X] [--rel-tol X] [--norm NAME] [--max-errors N]");
    Options o; o.first = argv[2]; o.second = argv[3];
    bool customNorm = false;
    for (int i = 4; i < argc; ++i) {
        std::string a = argv[i];
        auto value = [&]() -> std::string {
            if (++i >= argc) throw std::runtime_error("missing value after " + a);
            return argv[i];
        };
        if (a == "--all") o.exhaustiveAll = true;
        else if (a == "--exhaustive-at") o.exhaustiveAt = std::stoi(value());
        else if (a == "--samples") o.samples = std::stoull(value());
        else if (a == "--window-bins") o.windowBins = std::stoull(value());
        else if (a == "--seed") o.seed = std::stoull(value());
        else if (a == "--abs-tol") o.absTol = std::stod(value());
        else if (a == "--rel-tol") o.relTol = std::stod(value());
        else if (a == "--max-errors") o.maxErrors = std::stoull(value());
        else if (a == "--norm") {
            if (!customNorm) { o.norms.clear(); customNorm = true; }
            o.norms.push_back(value());
        } else throw std::runtime_error("unknown compare option: " + a);
    }
    if (o.exhaustiveAt <= 0 || !o.samples || !o.windowBins || o.absTol < 0 || o.relTol < 0)
        throw std::runtime_error("compare sizes and tolerances must be positive (absolute tolerance may be zero)");
    return o;
}

} // namespace

int compareMain(int argc, char *argv[]) {
    try {
        Options o = parse(argc, argv);
        Summary s;
        bool firstV10 = straw_v10::isV10(o.first), secondV10 = straw_v10::isV10(o.second);
        auto ca = chromosomeMap(o.first), cb = chromosomeMap(o.second);
        std::set<std::string> names;
        for (const auto &x : ca) names.insert(x.first);
        for (const auto &x : cb) names.insert(x.first);
        for (const auto &name : names) {
            ++s.metadata;
            if (!ca.count(name)) difference(s, o, "chromosome " + name + " only in second file");
            else if (!cb.count(name)) difference(s, o, "chromosome " + name + " only in first file");
            else if (ca[name].length != cb[name].length)
                difference(s, o, "chromosome " + name + " length " + std::to_string(ca[name].length) +
                                  " != " + std::to_string(cb[name].length));
        }
        std::vector<std::string> commonChroms;
        for (const auto &x : ca) if (cb.count(x.first)) commonChroms.push_back(x.first);

        auto ra = getResolutionsForFile(o.first), rb = getResolutionsForFile(o.second);
        std::set<int32_t> sa(ra.begin(), ra.end()), sb(rb.begin(), rb.end());
        for (int32_t r : sa) if (!sb.count(r)) difference(s, o, "resolution " + std::to_string(r) + " only in first file");
        for (int32_t r : sb) if (!sa.count(r)) difference(s, o, "resolution " + std::to_string(r) + " only in second file");
        std::vector<int32_t> resolutions;
        for (int32_t r : sa) if (sb.count(r)) resolutions.push_back(r);
        if (resolutions.empty()) throw std::runtime_error("files have no shared BP resolution");
        int32_t coarsest = *std::max_element(resolutions.begin(), resolutions.end());
        const auto &v10Chroms = firstV10 ? ca : cb;

        std::cout << "Comparing " << commonChroms.size() << " chromosomes at " << resolutions.size()
                  << " shared BP resolutions (seed " << o.seed << ")\n";
        for (int32_t resolution : resolutions) {
            bool exhaustive = o.exhaustiveAll || resolution >= o.exhaustiveAt || resolution == coarsest;
            std::cout << "  " << resolution << " BP: " << (exhaustive ? "exhaustive" : "sampled") << '\n';
            for (size_t i = 0; i < commonChroms.size(); ++i) {
                for (size_t j = i; j < commonChroms.size(); ++j) {
                    const auto &x = ca[commonChroms[i]], &y = ca[commonChroms[j]];
                    uint64_t xlen = std::min<uint64_t>(x.length, cb[x.name].length);
                    uint64_t ylen = std::min<uint64_t>(y.length, cb[y.name].length);
                    ++s.matrices;
                    auto check = [&](uint64_t x0, uint64_t x1, uint64_t y0, uint64_t y1) {
                        ++s.windows;
                        std::string context = x.name + "/" + y.name + " @" + std::to_string(resolution) +
                            " [" + std::to_string(x0) + "," + std::to_string(x1) + ")x[" +
                            std::to_string(y0) + "," + std::to_string(y1) + ")";
                        compareCells(rawCells(o.first, x.name, y.name, resolution, x0, x1, y0, y1,
                                              firstV10,
                                              !firstV10 && ca[x.name].index > ca[y.name].index),
                                     rawCells(o.second, x.name, y.name, resolution, x0, x1, y0, y1,
                                              secondV10,
                                              !secondV10 && cb[x.name].index > cb[y.name].index),
                                     o, s, context);
                    };
                    if (exhaustive) check(0, xlen, 0, ylen);
                    else {
                        uint64_t wx = std::min<uint64_t>(xlen, o.windowBins * uint64_t(resolution));
                        uint64_t wy = std::min<uint64_t>(ylen, o.windowBins * uint64_t(resolution));
                        std::mt19937_64 rng(o.seed ^ uint64_t(resolution) ^
                            std::hash<std::string>{}(x.name + "\0" + y.name));
                        for (size_t k = 0; k < o.samples; ++k) {
                            uint64_t x0 = xlen > wx ? (rng() % ((xlen - wx) / resolution + 1)) * resolution : 0;
                            uint64_t y0 = ylen > wy ? (rng() % ((ylen - wy) / resolution + 1)) * resolution : 0;
                            check(x0, std::min(xlen, x0 + wx), y0, std::min(ylen, y0 + wy));
                        }
                    }
                }
            }
            for (const auto &chr : commonChroms) {
                std::vector<double> a, b;
                size_t expectedLength = 0;
                for (const auto &entry : v10Chroms)
                    if (entry.first != "ALL" && entry.first != "All" && entry.first != "all")
                        expectedLength = std::max<size_t>(
                            expectedLength,
                            entry.second.length / uint64_t(resolution) +
                                (entry.second.length % uint64_t(resolution) != 0));
                size_t normalizationLength =
                    v10Chroms.at(chr).length / uint64_t(resolution) +
                    (v10Chroms.at(chr).length % uint64_t(resolution) != 0);
                bool aa = getExpectedVectorForFile(o.first, chr, resolution, "NONE", a);
                bool bb = getExpectedVectorForFile(o.second, chr, resolution, "NONE", b);
                if (aa != bb) difference(s, o, "expected NONE availability differs for " + chr + " @" + std::to_string(resolution));
                else if (aa) compareVector("expected", chr, resolution, "NONE", a, b, exhaustive,
                                           firstV10, secondV10, expectedLength, o, s);
                for (const auto &norm : o.norms) {
                    aa = getNormalizationVectorForFile(o.first, chr, resolution, norm, a);
                    bb = getNormalizationVectorForFile(o.second, chr, resolution, norm, b);
                    if (aa != bb) difference(s, o, "normalization " + norm + " availability differs for " + chr + " @" + std::to_string(resolution));
                    else if (aa) compareVector("normalization", chr, resolution, norm, a, b, exhaustive,
                                               firstV10, secondV10, normalizationLength, o, s);
                    aa = getExpectedVectorForFile(o.first, chr, resolution, norm, a);
                    bb = getExpectedVectorForFile(o.second, chr, resolution, norm, b);
                    if (aa != bb) difference(s, o, "normalized expected " + norm + " availability differs for " + chr + " @" + std::to_string(resolution));
                    else if (aa) compareVector("expected", chr, resolution, norm, a, b, exhaustive,
                                               firstV10, secondV10, expectedLength, o, s);
                }
            }
        }
        std::cout << "Checked " << s.matrices << " chromosome-pair/resolution matrices, " << s.windows
                  << " regions, " << s.cells << " nonzero-cell union entries, " << s.vectors
                  << " vectors, and " << s.vectorValues << " vector values.\n";
        if (s.differences) {
            std::cout << "RESULT: DIFFERENT (" << s.differences << " differences)\n";
            return 1;
        }
        std::cout << "RESULT: IDENTICAL WITHIN TOLERANCE\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "compare error: " << e.what() << '\n';
        return 2;
    }
}
