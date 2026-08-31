#include "compare.h"
#include "straw.h"
#include "straw_v10.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

double since(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

struct Options {
    std::string first, second;
    int32_t exhaustiveAt = 100000;
    size_t samples = 4;
    uint64_t windowBins = 256;
    uint64_t minWindowBins = 0;   // 0 means "use windowBins"
    uint64_t maxWindowBins = 0;   // 0 means "use windowBins"
    bool varyDistance = false;
    uint64_t minDistance = 0;
    uint64_t maxDistance = 0;     // 0 means "up to the chromosome length"
    size_t sampleResolutions = 0; // 0 means "every shared resolution"
    std::vector<int32_t> onlyResolutions;
    size_t randomRegions = 0;
    double intraFraction = 0.5;
    size_t stratified = 0;
    uint64_t nearMax = 1000000;
    uint64_t midMax = 10000000;
    bool skipVectors = false;
    size_t repeat = 1;
    uint64_t seed = 1;
    double absTol = 0.0;
    double relTol = 1e-6;
    size_t maxErrors = 20;
    bool exhaustiveAll = false;
    bool timing = true;
    std::string timingCsv;
    std::vector<std::string> norms{"VC", "VC_SQRT", "KR"};

    uint64_t windowLow() const {
        return std::max<uint64_t>(1, minWindowBins ? minWindowBins : windowBins);
    }
    uint64_t windowHigh() const {
        return std::max(windowLow(), maxWindowBins ? maxWindowBins : windowBins);
    }
    bool varySize() const { return windowHigh() > windowLow(); }
};

struct Summary {
    uint64_t metadata = 0, matrices = 0, windows = 0, cells = 0;
    uint64_t vectors = 0, vectorValues = 0, differences = 0;
};

using Cell = std::pair<uint64_t, uint64_t>;
using Cells = std::map<Cell, double>;

// ---------------------------------------------------------------------------
// Timing
// ---------------------------------------------------------------------------

struct Stat {
    std::vector<double> values;
    double total = 0;
    uint64_t records = 0;

    void add(double seconds, uint64_t recordCount) {
        values.push_back(seconds);
        total += seconds;
        records += recordCount;
    }
    bool empty() const { return values.empty(); }
    double mean() const { return values.empty() ? 0.0 : total / double(values.size()); }
    double quantile(double q) const {
        if (values.empty()) return 0.0;
        std::vector<double> sorted = values;
        std::sort(sorted.begin(), sorted.end());
        size_t i = size_t(q * double(sorted.size() - 1) + 0.5);
        return sorted[std::min(i, sorted.size() - 1)];
    }
};

struct FileTiming {
    Stat region, normalization, expected;
    std::map<int32_t, Stat> byResolution;
    std::map<int, Stat> bySize;      // key: floor(log2(width in bins))
    std::map<int, Stat> byDistance;  // key: floor(log10(bp)), -1 inter, -2 diagonal
    std::map<std::pair<int32_t, int>, Stat> byStratum;  // key: (resolution, stratum)
};

// Strata for the balanced sweep: three intra-chromosomal distance bands plus
// inter-chromosomal, so a benchmark covers near-diagonal blocks, mid-range
// blocks, the sparse far corner, and trans contacts in equal measure.
enum Stratum { STRATUM_NONE = -1, STRATUM_NEAR = 0, STRATUM_MID = 1, STRATUM_FAR = 2,
               STRATUM_INTER = 3, STRATUM_COUNT = 4 };

const char *stratumName(int s) {
    switch (s) {
        case STRATUM_NEAR: return "near";
        case STRATUM_MID: return "mid";
        case STRATUM_FAR: return "far";
        case STRATUM_INTER: return "inter";
        default: return "-";
    }
}

// Half-open [lo, hi) diagonal-offset band in base pairs; hi == 0 means "up to
// the chromosome length".
struct Band {
    bool active = false;
    uint64_t lo = 0, hi = 0;
};

Band stratumBand(int s, const Options &o) {
    Band b;
    b.active = true;
    if (s == STRATUM_NEAR) { b.lo = 0; b.hi = o.nearMax; }
    else if (s == STRATUM_MID) { b.lo = o.nearMax; b.hi = o.midMax; }
    else { b.lo = o.midMax; b.hi = 0; }
    return b;
}

struct Region {
    std::string chr1, chr2;
    int32_t resolution = 0;
    uint64_t x0 = 0, x1 = 0, y0 = 0, y1 = 0;
    uint64_t widthBins = 0;
    int64_t distance = -1;  // -1 for inter-chromosomal pairs
    int stratum = STRATUM_NONE;
};

int sizeBucket(uint64_t widthBins) {
    int b = 0;
    while ((uint64_t(1) << (b + 1)) <= widthBins && b < 62) ++b;
    return b;
}

int distanceBucket(int64_t distance) {
    if (distance < 0) return -1;
    if (distance == 0) return -2;
    return int(std::floor(std::log10(double(distance))));
}

std::string sizeBucketLabel(int b) {
    std::ostringstream o;
    o << (uint64_t(1) << b) << '-' << ((uint64_t(1) << (b + 1)) - 1) << " bins";
    return o.str();
}

std::string distanceBucketLabel(int b) {
    if (b == -1) return "inter-chromosomal";
    if (b == -2) return "diagonal (d=0)";
    std::ostringstream o;
    o << "1e" << b << "-1e" << (b + 1) << " bp";
    return o.str();
}

// Collects per-query wall-clock timings and optionally writes one CSV row per
// query so that size/distance/resolution sweeps can be analysed offline.
struct Bench {
    const Options *o = nullptr;
    FileTiming files[2];
    std::set<int32_t> derived[2];  // resolutions computed on the fly
    std::ofstream csv;
    double wall = 0;

    void open(const Options &options) {
        o = &options;
        if (options.timingCsv.empty()) return;
        csv.open(options.timingCsv);
        if (!csv) throw std::runtime_error("cannot write timing CSV: " + options.timingCsv);
        csv << "type,file,path,chrom1,chrom2,resolution,derived,stratum,norm,x0,x1,y0,y1,"
               "width_bins,distance,records,seconds,iteration\n";
        csv << std::setprecision(9);
    }

    bool isDerived(int which, int32_t resolution) const {
        return derived[which].count(resolution) != 0;
    }
    // A resolution is only interesting as "derived" if it costs one of the two
    // readers extra work, so label it when either file computes it on the fly.
    std::string resolutionLabel(int32_t resolution) const {
        std::string label = std::to_string(resolution) + " BP";
        bool a = isDerived(0, resolution), b = isDerived(1, resolution);
        if (a && b) label += " [derived]";
        else if (a) label += " [derived:1st]";
        else if (b) label += " [derived:2nd]";
        return label;
    }

    const std::string &path(int which) const { return which ? o->second : o->first; }

    void region(int which, const Region &r, double seconds, uint64_t records, size_t iteration) {
        FileTiming &f = files[which];
        f.region.add(seconds, records);
        f.byResolution[r.resolution].add(seconds, records);
        f.bySize[sizeBucket(r.widthBins)].add(seconds, records);
        f.byDistance[distanceBucket(r.distance)].add(seconds, records);
        if (r.stratum != STRATUM_NONE)
            f.byStratum[{r.resolution, r.stratum}].add(seconds, records);
        if (!csv) return;
        csv << "region," << (which ? "second" : "first") << ',' << path(which) << ',' << r.chr1
            << ',' << r.chr2 << ',' << r.resolution << ',' << (isDerived(which, r.resolution) ? 1 : 0)
            << ',' << stratumName(r.stratum) << ",NONE," << r.x0 << ',' << r.x1 << ','
            << r.y0 << ',' << r.y1 << ',' << r.widthBins << ',' << r.distance << ',' << records
            << ',' << seconds << ',' << iteration << '\n';
    }

    void vector(int which, const std::string &kind, const std::string &chr, int32_t resolution,
                const std::string &norm, double seconds, uint64_t values) {
        FileTiming &f = files[which];
        (kind == "normalization" ? f.normalization : f.expected).add(seconds, values);
        if (!csv) return;
        csv << kind << ',' << (which ? "second" : "first") << ',' << path(which) << ',' << chr
            << ",," << resolution << ',' << (isDerived(which, resolution) ? 1 : 0) << ",-," << norm
            << ",,,,,,," << values << ',' << seconds << ",0\n";
    }
};

void statRow(const std::string &label, const std::string &file, const Stat &s,
             const std::string &unit) {
    if (s.empty()) return;
    std::cout << "    " << std::left << std::setw(24) << label << std::setw(8) << file << std::right
              << std::setw(7) << s.values.size() << std::fixed << std::setprecision(3)
              << std::setw(11) << s.total << std::setw(10) << s.mean() * 1e3 << std::setw(10)
              << s.quantile(0.5) * 1e3 << std::setw(10) << s.quantile(0.95) * 1e3 << std::setw(10)
              << s.quantile(1.0) * 1e3;
    if (s.total > 0)
        std::cout << std::setw(14) << std::setprecision(0) << double(s.records) / s.total;
    else
        std::cout << std::setw(14) << '-';
    std::cout << "  " << unit << '\n';
    std::cout << std::setprecision(6);
}

void statHeader(const std::string &title) {
    std::cout << "  " << title << '\n'
              << "    " << std::left << std::setw(24) << "bucket" << std::setw(8) << "file"
              << std::right << std::setw(7) << "n" << std::setw(11) << "total(s)" << std::setw(10)
              << "mean(ms)" << std::setw(10) << "p50(ms)" << std::setw(10) << "p95(ms)"
              << std::setw(10) << "max(ms)" << std::setw(14) << "rate/s" << '\n';
}

void report(const Bench &bench, const Options &o) {
    std::cout << "\nTiming (wall clock per query; each query re-opens the file and parses its "
                 "header)\n";
    std::cout << "  first  = " << o.first << "\n  second = " << o.second << '\n';
    statHeader("Totals by operation");
    const char *names[2] = {"first", "second"};
    for (int w = 0; w < 2; ++w) statRow("region read", names[w], bench.files[w].region, "records/s");
    for (int w = 0; w < 2; ++w)
        statRow("normalization vector", names[w], bench.files[w].normalization, "values/s");
    for (int w = 0; w < 2; ++w)
        statRow("expected vector", names[w], bench.files[w].expected, "values/s");

    double a = bench.files[0].region.total, b = bench.files[1].region.total;
    if (a > 0 && b > 0) {
        std::cout << "  Region-read total: first " << std::fixed << std::setprecision(3) << a
                  << "s vs second " << b << "s ("
                  << (a < b ? "first " : "second ") << std::setprecision(2)
                  << (a < b ? b / a : a / b) << "x faster)\n"
                  << std::setprecision(6);
    }

    std::set<int32_t> resolutions;
    std::set<int> sizes, distances;
    for (int w = 0; w < 2; ++w) {
        for (const auto &e : bench.files[w].byResolution) resolutions.insert(e.first);
        for (const auto &e : bench.files[w].bySize) sizes.insert(e.first);
        for (const auto &e : bench.files[w].byDistance) distances.insert(e.first);
    }
    if (resolutions.size() > 1 || sizes.size() > 1 || distances.size() > 1) std::cout << '\n';
    if (resolutions.size() > 1) {
        statHeader("Region reads by resolution");
        for (int32_t r : resolutions)
            for (int w = 0; w < 2; ++w) {
                auto it = bench.files[w].byResolution.find(r);
                if (it != bench.files[w].byResolution.end())
                    statRow(bench.resolutionLabel(r), names[w], it->second, "records/s");
            }
    }
    if (sizes.size() > 1) {
        statHeader("Region reads by size");
        for (int s : sizes)
            for (int w = 0; w < 2; ++w) {
                auto it = bench.files[w].bySize.find(s);
                if (it != bench.files[w].bySize.end())
                    statRow(sizeBucketLabel(s), names[w], it->second, "records/s");
            }
    }
    if (distances.size() > 1) {
        statHeader("Region reads by diagonal distance");
        for (int d : distances)
            for (int w = 0; w < 2; ++w) {
                auto it = bench.files[w].byDistance.find(d);
                if (it != bench.files[w].byDistance.end())
                    statRow(distanceBucketLabel(d), names[w], it->second, "records/s");
            }
    }
    if (!bench.files[0].byStratum.empty()) {
        std::cout << "\n  Stratified head-to-head (median region-read latency)\n"
                  << "    " << std::left << std::setw(26) << "resolution" << std::setw(8)
                  << "stratum" << std::right << std::setw(7) << "n" << std::setw(12) << "first(ms)"
                  << std::setw(12) << "second(ms)" << std::setw(9) << "ratio" << "  winner\n";
        for (const auto &entry : bench.files[0].byStratum) {
            auto other = bench.files[1].byStratum.find(entry.first);
            if (other == bench.files[1].byStratum.end()) continue;
            double a = entry.second.quantile(0.5), b = other->second.quantile(0.5);
            std::cout << "    " << std::left << std::setw(26)
                      << bench.resolutionLabel(entry.first.first) << std::setw(8)
                      << stratumName(entry.first.second) << std::right << std::setw(7)
                      << entry.second.values.size() << std::fixed << std::setprecision(3)
                      << std::setw(12) << a * 1e3 << std::setw(12) << b * 1e3;
            if (a > 0 && b > 0)
                std::cout << std::setw(8) << std::setprecision(2) << (a < b ? b / a : a / b) << 'x'
                          << "  " << (a < b ? "first" : "second");
            else
                std::cout << std::setw(9) << '-' << "  -";
            std::cout << '\n' << std::setprecision(3);
        }
        std::cout << std::setprecision(6);
    }
    std::cout << "  Total elapsed: " << std::fixed << std::setprecision(3) << bench.wall << "s\n"
              << std::setprecision(6);
    if (!o.timingCsv.empty()) std::cout << "  Per-query timings written to " << o.timingCsv << '\n';
}

// ---------------------------------------------------------------------------
// Sampling
// ---------------------------------------------------------------------------

double unit01(std::mt19937_64 &rng) {
    return double(rng() >> 11) * (1.0 / 9007199254740992.0);
}

// Log-uniform draw so that a sweep covers small and large values evenly on a
// log scale instead of being dominated by the top of the range.
uint64_t logUniform(std::mt19937_64 &rng, uint64_t lo, uint64_t hi) {
    if (hi <= lo) return lo;
    if (lo == 0) return logUniform(rng, 1, hi + 1) - 1;
    double l = std::log(double(lo)), h = std::log(double(hi));
    double v = std::exp(l + (h - l) * unit01(rng));
    uint64_t out = uint64_t(v + 0.5);
    return std::min(hi, std::max(lo, out));
}

// A chromosome can host a band only if some region still fits beyond its lower
// diagonal offset.
bool bandFits(uint64_t length, int32_t resolution, const Band &band) {
    if (!band.active || length <= band.lo) return !band.active || band.lo == 0;
    return (length - band.lo) / uint64_t(resolution) >= 1;
}

Region sampleRegion(std::mt19937_64 &rng, const Options &o, const std::string &chr1,
                    const std::string &chr2, int32_t resolution, uint64_t xlen, uint64_t ylen,
                    const Band &band = Band()) {
    const uint64_t res = uint64_t(resolution);
    const bool intra = chr1 == chr2;
    const bool useBand = intra && band.active;

    // A wide region pushes both corners apart, so cap the width at whatever
    // still leaves room for the band's lower offset.
    uint64_t widthHi = o.windowHigh();
    if (useBand && band.lo)
        widthHi = std::min(widthHi, xlen > band.lo ? (xlen - band.lo) / res : 0);
    widthHi = std::max<uint64_t>(1, widthHi);
    uint64_t widthLo = std::min(o.windowLow(), widthHi);
    const uint64_t widthBins = logUniform(rng, widthLo, widthHi);
    const uint64_t wx = std::min<uint64_t>(xlen, widthBins * res);
    const uint64_t wy = std::min<uint64_t>(ylen, widthBins * res);

    Region r;
    r.chr1 = chr1;
    r.chr2 = chr2;
    r.resolution = resolution;
    if (useBand) {
        // Both corners must stay inside the chromosome, so the reachable
        // diagonal offset shrinks as the window grows. The offset is drawn in
        // bins rather than base pairs: a log-uniform draw over base pairs would
        // collapse onto the diagonal once rounded to a coarse resolution.
        uint64_t spanBins = (xlen > wx ? xlen - wx : 0) / res;
        uint64_t hiBins = band.hi ? std::min<uint64_t>(band.hi / res, spanBins) : spanBins;
        uint64_t loBins = std::min<uint64_t>((band.lo + res - 1) / res, hiBins);
        uint64_t d = logUniform(rng, loBins, hiBins) * res;
        uint64_t slack = spanBins * res > d ? spanBins * res - d : 0;
        r.x0 = slack ? (rng() % (slack / res + 1)) * res : 0;
        r.y0 = r.x0 + d;
        r.distance = int64_t(d);
    } else {
        r.x0 = xlen > wx ? (rng() % ((xlen - wx) / res + 1)) * res : 0;
        r.y0 = ylen > wy ? (rng() % ((ylen - wy) / res + 1)) * res : 0;
        r.distance = intra ? std::llabs(int64_t(r.y0) - int64_t(r.x0)) : -1;
    }
    r.x1 = std::min(xlen, r.x0 + wx);
    r.y1 = std::min(ylen, r.y0 + wy);
    r.widthBins = std::max<uint64_t>(1, (r.x1 - r.x0 + res - 1) / res);
    return r;
}

// ---------------------------------------------------------------------------
// Comparison
// ---------------------------------------------------------------------------

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

struct Read {
    Cells cells;
    double seconds = 0;
    uint64_t records = 0;
};

Read rawCells(const std::string &path, const Region &r) {
    Read out;
    const int32_t resolution = r.resolution;
    auto add = [&](uint64_t x, uint64_t y, double value) {
        ++out.records;
        // Legacy queries can include a bin touching the inclusive region end.
        // Filtering here gives both readers identical half-open semantics.
        bool direct = x >= r.x0 && x < r.x1 && y >= r.y0 && y < r.y1;
        bool reflected = r.chr1 == r.chr2 && y >= r.x0 && y < r.x1 && x >= r.y0 && x < r.y1;
        if (!direct && !reflected) return;
        if (!direct) std::swap(x, y);
        out.cells[{x / static_cast<uint64_t>(resolution), y / static_cast<uint64_t>(resolution)}] +=
            value;
    };
    const std::string a = location(r.chr1, r.x0, r.x1), b = location(r.chr2, r.y0, r.y1);
    auto start = Clock::now();
    if (straw_v10::isV10(path)) {
        straw_v10::File(path).streamRaw(a, b, "BP", resolution,
            [&](const straw_v10::Record &record) {
                double value = record.isScore ? static_cast<double>(record.score)
                                              : static_cast<double>(record.count);
                add(uint64_t(record.binX) * resolution, uint64_t(record.binY) * resolution, value);
            });
    } else {
        strawStream("observed", "NONE", path, a, b, "BP", resolution,
            [&](const contactRecord &record) { add(record.binX, record.binY, record.counts); });
    }
    out.seconds = since(start);
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
    if (exhaustive || n <= o.windowHigh() * std::max<size_t>(1, o.samples)) {
        indices.resize(n);
        for (size_t i = 0; i < n; ++i) indices[i] = i;
    } else if (n) {
        std::mt19937_64 rng(o.seed ^ uint64_t(resolution) ^ std::hash<std::string>{}(chr + kind + norm));
        indices = {0, n - 1, n / 2};
        for (size_t k = 0; k < o.samples * o.windowHigh(); ++k) indices.push_back(rng() % n);
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

const char *usage() {
    return "usage: straw compare <first.hic> <second.hic> [options]\n"
           "  coverage:\n"
           "    --all                     compare every region exhaustively\n"
           "    --exhaustive-at BP        resolutions at or above BP are exhaustive (default 100000)\n"
           "    --samples N               sampled regions per chromosome pair (default 4)\n"
           "    --max-errors N            differences printed before suppression (default 20)\n"
           "    --norm NAME               normalization to check (repeatable; default VC VC_SQRT KR)\n"
           "  region sampling:\n"
           "    --window-bins N           fixed sampled region width in bins (default 256)\n"
           "    --min-window-bins N       sample region widths log-uniformly from N ...\n"
           "    --max-window-bins N       ... up to N bins\n"
           "    --vary-distance           sample the diagonal offset of intra regions log-uniformly\n"
           "    --min-distance BP         lower bound for the sampled diagonal offset\n"
           "    --max-distance BP         upper bound for the sampled diagonal offset\n"
           "    --resolution BP           restrict to this resolution (repeatable)\n"
           "    --sample-resolutions N    randomly keep N of the shared resolutions\n"
           "    --random-regions N        draw N random (pair, resolution, size, distance) regions\n"
           "                              instead of walking every chromosome pair\n"
           "    --intra-fraction P        share of random regions that are intra-chromosomal\n"
           "                              (default 0.5; only these carry a diagonal offset)\n"
           "    --stratified N            per resolution, draw N regions in each of four strata:\n"
           "                              near-diagonal, mid-range, far, and inter-chromosomal,\n"
           "                              cycling through distinct chromosomes\n"
           "    --near-max BP             near/mid boundary (default 1000000)\n"
           "    --mid-max BP              mid/far boundary (default 10000000)\n"
           "    --skip-vectors            compare regions only, skipping vector checks\n"
           "    --seed N                  random seed (default 1)\n"
           "  timing:\n"
           "    --repeat N                read each region N times, timing every read (default 1)\n"
           "    --timing-csv PATH         write one CSV row per query\n"
           "    --no-timing               suppress the timing report\n"
           "  tolerance:\n"
           "    --abs-tol X               absolute tolerance (default 0)\n"
           "    --rel-tol X               relative tolerance (default 1e-6)";
}

Options parse(int argc, char **argv) {
    if (argc < 4) throw std::runtime_error(usage());
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
        else if (a == "--min-window-bins") o.minWindowBins = std::stoull(value());
        else if (a == "--max-window-bins") o.maxWindowBins = std::stoull(value());
        else if (a == "--vary-distance") o.varyDistance = true;
        else if (a == "--min-distance") { o.minDistance = std::stoull(value()); o.varyDistance = true; }
        else if (a == "--max-distance") { o.maxDistance = std::stoull(value()); o.varyDistance = true; }
        else if (a == "--resolution") o.onlyResolutions.push_back(std::stoi(value()));
        else if (a == "--sample-resolutions") o.sampleResolutions = std::stoull(value());
        else if (a == "--random-regions") o.randomRegions = std::stoull(value());
        else if (a == "--intra-fraction") o.intraFraction = std::stod(value());
        else if (a == "--stratified") o.stratified = std::stoull(value());
        else if (a == "--near-max") o.nearMax = std::stoull(value());
        else if (a == "--mid-max") o.midMax = std::stoull(value());
        else if (a == "--skip-vectors") o.skipVectors = true;
        else if (a == "--repeat") o.repeat = std::stoull(value());
        else if (a == "--timing-csv") o.timingCsv = value();
        else if (a == "--no-timing") o.timing = false;
        else if (a == "--seed") o.seed = std::stoull(value());
        else if (a == "--abs-tol") o.absTol = std::stod(value());
        else if (a == "--rel-tol") o.relTol = std::stod(value());
        else if (a == "--max-errors") o.maxErrors = std::stoull(value());
        else if (a == "--norm") {
            if (!customNorm) { o.norms.clear(); customNorm = true; }
            o.norms.push_back(value());
        } else throw std::runtime_error("unknown compare option: " + a + "\n" + usage());
    }
    if (o.exhaustiveAt <= 0 || !o.samples || !o.windowBins || !o.repeat || o.absTol < 0 || o.relTol < 0)
        throw std::runtime_error("compare sizes and tolerances must be positive (absolute tolerance may be zero)");
    if (o.minWindowBins && o.maxWindowBins && o.minWindowBins > o.maxWindowBins)
        throw std::runtime_error("--min-window-bins must not exceed --max-window-bins");
    if (o.maxDistance && o.minDistance > o.maxDistance)
        throw std::runtime_error("--min-distance must not exceed --max-distance");
    if (o.intraFraction < 0 || o.intraFraction > 1)
        throw std::runtime_error("--intra-fraction must be between 0 and 1");
    if (o.nearMax >= o.midMax)
        throw std::runtime_error("--near-max must be smaller than --mid-max");
    return o;
}

} // namespace

int compareMain(int argc, char *argv[]) {
    try {
        auto started = Clock::now();
        Options o = parse(argc, argv);
        Summary s;
        Bench bench;
        bench.open(o);
        bool firstV10 = straw_v10::isV10(o.first), secondV10 = straw_v10::isV10(o.second);
        // Only V10 can synthesize a resolution at query time; V6-V9 stores
        // every advertised resolution.
        if (firstV10)
            for (int32_t r : straw_v10::File(o.first).derivedResolutions()) bench.derived[0].insert(r);
        if (secondV10)
            for (int32_t r : straw_v10::File(o.second).derivedResolutions()) bench.derived[1].insert(r);
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
        if (!o.onlyResolutions.empty()) {
            std::set<int32_t> wanted(o.onlyResolutions.begin(), o.onlyResolutions.end());
            std::vector<int32_t> kept;
            for (int32_t r : resolutions) if (wanted.count(r)) kept.push_back(r);
            for (int32_t r : wanted)
                if (!std::count(resolutions.begin(), resolutions.end(), r))
                    throw std::runtime_error("resolution " + std::to_string(r) + " is not shared by both files");
            resolutions = kept;
        }
        if (o.sampleResolutions && o.sampleResolutions < resolutions.size()) {
            std::mt19937_64 rng(o.seed ^ 0x9e3779b97f4a7c15ull);
            std::shuffle(resolutions.begin(), resolutions.end(), rng);
            resolutions.resize(o.sampleResolutions);
            std::sort(resolutions.begin(), resolutions.end());
        }
        int32_t coarsest = *std::max_element(resolutions.begin(), resolutions.end());
        const auto &v10Chroms = firstV10 ? ca : cb;

        auto chromLength = [&](const std::string &name) {
            return std::min<uint64_t>(ca[name].length, cb[name].length);
        };

        // Reads both files for one region (repeated if requested), records the
        // timings, and compares the first pair of reads.
        auto check = [&](const Region &r) {
            ++s.windows;
            std::string context = r.chr1 + "/" + r.chr2 + " @" + std::to_string(r.resolution) +
                " [" + std::to_string(r.x0) + "," + std::to_string(r.x1) + ")x[" +
                std::to_string(r.y0) + "," + std::to_string(r.y1) + ")";
            Cells first, second;
            for (size_t iteration = 0; iteration < o.repeat; ++iteration) {
                // Whichever file is read second benefits from a warm page
                // cache, so alternate the order: a systematic bias would
                // otherwise show up as a spurious speed difference.
                bool secondLeads = (s.windows + iteration) % 2 == 0;
                Read lead = rawCells(secondLeads ? o.second : o.first, r);
                Read trail = rawCells(secondLeads ? o.first : o.second, r);
                Read &a = secondLeads ? trail : lead;
                Read &b = secondLeads ? lead : trail;
                bench.region(0, r, a.seconds, a.records, iteration);
                bench.region(1, r, b.seconds, b.records, iteration);
                if (iteration == 0) { first = std::move(a.cells); second = std::move(b.cells); }
            }
            compareCells(first, second, o, s, context);
        };

        std::cout << "Comparing " << commonChroms.size() << " chromosomes at " << resolutions.size()
                  << " shared BP resolutions (seed " << o.seed << ")\n";
        if (o.varySize())
            std::cout << "  Region widths sampled log-uniformly over " << o.windowLow() << '-'
                      << o.windowHigh() << " bins\n";
        if (o.varyDistance)
            std::cout << "  Intra-chromosomal diagonal offsets sampled log-uniformly over "
                      << o.minDistance << " bp - "
                      << (o.maxDistance ? std::to_string(o.maxDistance) + " bp"
                                        : std::string("the chromosome length"))
                      << '\n';

        if (o.stratified) {
            std::cout << "  Stratified sweep: " << o.stratified << " regions per stratum per "
                      << "resolution (near <" << o.nearMax << " bp, mid " << o.nearMax << '-'
                      << o.midMax << " bp, far >" << o.midMax << " bp, inter)\n";
            for (int32_t resolution : resolutions) {
                for (int stratum = 0; stratum < STRATUM_COUNT; ++stratum) {
                    // Independent stream per cell keeps a stratum reproducible
                    // regardless of which other strata ran.
                    std::mt19937_64 rng(o.seed ^ uint64_t(resolution) * 0x9e3779b97f4a7c15ull ^
                                        uint64_t(stratum) * 0xbf58476d1ce4e5b9ull);
                    Band band = stratumBand(stratum, o);
                    std::vector<std::string> eligible;
                    for (const auto &chr : commonChroms) {
                        if (stratum == STRATUM_INTER) eligible.push_back(chr);
                        else if (bandFits(chromLength(chr), resolution, band)) eligible.push_back(chr);
                    }
                    if (stratum == STRATUM_INTER && eligible.size() < 2) eligible.clear();
                    if (eligible.empty()) {
                        std::cout << "    " << resolution << " BP " << stratumName(stratum)
                                  << ": skipped, no chromosome is long enough\n";
                        continue;
                    }
                    // Walking a shuffled list rather than sampling with
                    // replacement spreads the draws over distinct chromosomes.
                    std::shuffle(eligible.begin(), eligible.end(), rng);
                    for (size_t k = 0; k < o.stratified; ++k) {
                        const std::string &x = eligible[k % eligible.size()];
                        ++s.matrices;
                        if (stratum == STRATUM_INTER) {
                            size_t offset = 1 + rng() % (eligible.size() - 1);
                            const std::string &y = eligible[(k + offset) % eligible.size()];
                            Region r = sampleRegion(rng, o, std::min(x, y), std::max(x, y),
                                                    resolution, chromLength(std::min(x, y)),
                                                    chromLength(std::max(x, y)));
                            r.stratum = stratum;
                            check(r);
                        } else {
                            Region r = sampleRegion(rng, o, x, x, resolution, chromLength(x),
                                                    chromLength(x), band);
                            r.stratum = stratum;
                            check(r);
                        }
                    }
                }
            }
        }

        if (o.randomRegions) {
            std::cout << "  Drawing " << o.randomRegions << " random regions across "
                      << resolutions.size() << " resolutions\n";
            std::mt19937_64 rng(o.seed ^ 0xd1b54a32d192ed03ull);
            for (size_t k = 0; k < o.randomRegions; ++k) {
                int32_t resolution = resolutions[rng() % resolutions.size()];
                size_t i = rng() % commonChroms.size(), j;
                // Diagonal-offset sampling only means something within a
                // chromosome, and uniform pairs are almost never intra, so the
                // intra share is drawn explicitly.
                if (unit01(rng) < o.intraFraction) j = i;
                else {
                    j = rng() % commonChroms.size();
                    if (i > j) std::swap(i, j);
                }
                const std::string &x = commonChroms[i], &y = commonChroms[j];
                ++s.matrices;
                check(sampleRegion(rng, o, x, y, resolution, chromLength(x), chromLength(y)));
            }
        }

        const bool sweeping = o.randomRegions || o.stratified;
        for (int32_t resolution : resolutions) {
            bool exhaustive = o.exhaustiveAll || resolution >= o.exhaustiveAt || resolution == coarsest;
            if (o.skipVectors && sweeping) continue;
            std::cout << "  " << bench.resolutionLabel(resolution) << ": "
                      << (sweeping ? "vectors only" : (exhaustive ? "exhaustive" : "sampled"))
                      << '\n';
            if (!sweeping) {
                for (size_t i = 0; i < commonChroms.size(); ++i) {
                    for (size_t j = i; j < commonChroms.size(); ++j) {
                        const std::string &x = commonChroms[i], &y = commonChroms[j];
                        uint64_t xlen = chromLength(x), ylen = chromLength(y);
                        ++s.matrices;
                        if (exhaustive) {
                            Region r;
                            r.chr1 = x; r.chr2 = y; r.resolution = resolution;
                            r.x0 = 0; r.x1 = xlen; r.y0 = 0; r.y1 = ylen;
                            r.widthBins = (xlen + uint64_t(resolution) - 1) / uint64_t(resolution);
                            r.distance = x == y ? 0 : -1;
                            check(r);
                        } else {
                            std::mt19937_64 rng(o.seed ^ uint64_t(resolution) ^
                                std::hash<std::string>{}(x + "\0" + y));
                            for (size_t k = 0; k < o.samples; ++k)
                                check(sampleRegion(rng, o, x, y, resolution, xlen, ylen));
                        }
                    }
                }
            }
            if (o.skipVectors) continue;
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
                auto fetch = [&](int which, const std::string &kind, const std::string &norm,
                                 std::vector<double> &out) {
                    const std::string &path = which ? o.second : o.first;
                    auto start = Clock::now();
                    bool ok = kind == "normalization"
                        ? getNormalizationVectorForFile(path, chr, resolution, norm, out)
                        : getExpectedVectorForFile(path, chr, resolution, norm, out);
                    bench.vector(which, kind, chr, resolution, norm, since(start),
                                 ok ? out.size() : 0);
                    return ok;
                };
                bool aa = fetch(0, "expected", "NONE", a);
                bool bb = fetch(1, "expected", "NONE", b);
                if (aa != bb) difference(s, o, "expected NONE availability differs for " + chr + " @" + std::to_string(resolution));
                else if (aa) compareVector("expected", chr, resolution, "NONE", a, b, exhaustive,
                                           firstV10, secondV10, expectedLength, o, s);
                for (const auto &norm : o.norms) {
                    aa = fetch(0, "normalization", norm, a);
                    bb = fetch(1, "normalization", norm, b);
                    if (aa != bb) difference(s, o, "normalization " + norm + " availability differs for " + chr + " @" + std::to_string(resolution));
                    else if (aa) compareVector("normalization", chr, resolution, norm, a, b, exhaustive,
                                               firstV10, secondV10, normalizationLength, o, s);
                    aa = fetch(0, "expected", norm, a);
                    bb = fetch(1, "expected", norm, b);
                    if (aa != bb) difference(s, o, "normalized expected " + norm + " availability differs for " + chr + " @" + std::to_string(resolution));
                    else if (aa) compareVector("expected", chr, resolution, norm, a, b, exhaustive,
                                               firstV10, secondV10, expectedLength, o, s);
                }
            }
        }
        std::cout << "Checked " << s.matrices << " chromosome-pair/resolution matrices, " << s.windows
                  << " regions, " << s.cells << " nonzero-cell union entries, " << s.vectors
                  << " vectors, and " << s.vectorValues << " vector values.\n";
        bench.wall = since(started);
        if (o.timing) report(bench, o);
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
