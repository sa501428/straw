#include "subsample.h"
#include "straw.h"
#include "straw_v10.h"
#include "hbs.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>

namespace {
const char* usage =
    "Usage: straw subsample <hicFile> <--fraction P|--contacts N> "
    "[--resolution BP] [--seed N] [--output output.hbs.gz]\n"
    "Prints chr1 pos1 chr2 pos2 count to stdout (tab-separated).\n"
    "Default: finest nonempty BP resolution, seed 1, all real cis/trans pairs.\n"
    "--contacts sets probability N / total at the coarsest BP resolution;\n"
    "the retained total is random, not an exact target.\n";

uint64_t unsignedValue(const std::string& value, const std::string& option) {
    if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
        throw std::runtime_error(option + " requires an unsigned integer");
    size_t used = 0;
    auto n = std::stoull(value, &used);
    if (used != value.size()) throw std::runtime_error("Invalid " + option);
    return n;
}

uint64_t integerCount(float value) {
    // 2^64 is exactly representable, unlike UINT64_MAX in floating point.
    if (!std::isfinite(value) || value < 0 || std::floor(value) != value ||
        static_cast<double>(value) >= std::ldexp(1.0, 64))
        throw std::runtime_error("Subsampling requires finite, nonnegative integer raw counts");
    return static_cast<uint64_t>(value);
}

uint64_t binomialChunk(uint64_t n, double p, std::mt19937_64& rng) {
    uint64_t kept = 0;
    // Split at the a-th order statistic of n independent uniforms:
    // B ~ Beta(a,n+1-a). Conditional on B, recurse into the side containing p.
    // This is a binomial draw, not a normal/Poisson approximation. It avoids
    // libc++'s modal binomial sampler walking O(sqrt(n)) cells for huge counts.
    while (n > 64 && p > 0 && p < 1) {
        const uint64_t a = n / 2 + 1;
        const uint64_t b = n - a + 1;
        const double x = std::gamma_distribution<double>(static_cast<double>(a), 1)(rng);
        const double y = std::gamma_distribution<double>(static_cast<double>(b), 1)(rng);
        const double split = x / (x + y);
        if (p < split) {
            n = a - 1;
            p /= split;
        } else {
            kept += a;
            n -= a;
            p = (p - split) / (1 - split);
        }
    }
    if (p <= 0) return kept;
    if (p >= 1) return kept + n;
    return kept + std::binomial_distribution<int64_t>(static_cast<int64_t>(n), p)(rng);
}

uint64_t sample(uint64_t n, double p, std::mt19937_64& rng) {
    if (p == 0) return 0;
    if (p == 1) return n;
    // Independent binomials sum to Binomial(n,p). Keep each trial count exactly
    // representable in double and safely inside the distribution's signed type.
    // This also handles V10 uint64 counts without a per-read loop.
    const uint64_t chunk = uint64_t{1} << 52;
    uint64_t kept = 0;
    while (n) {
        auto trials = std::min(n, chunk);
        kept += binomialChunk(trials, p, rng);
        n -= trials;
    }
    return kept;
}
} // namespace

static int exportCounts(const std::string& path, double probability, bool contactsSet,
                        uint64_t target, uint64_t seed, int32_t resolution, bool resolutionSet,
                        const std::string& output, ContactFilter filter);

int subsampleMain(int argc, char* argv[]) {
    if (argc == 3 && std::string(argv[2]) == "--help") {
        std::cout << usage;
        return 0;
    }
    if (argc < 5) throw std::runtime_error(usage);
    const std::string path = argv[2];
    bool fractionSet = false, contactsSet = false, resolutionSet = false, seedSet = false;
    double probability = 1;
    uint64_t target = 0, seed = 1;
    int32_t resolution = 0;
    std::string output;
    for (int i = 3; i < argc; i += 2) {
        const std::string option = argv[i];
        if (i + 1 == argc) throw std::runtime_error("Missing value for " + option);
        const std::string value = argv[i + 1];
        if (option == "--fraction" && !fractionSet) {
            size_t used = 0;
            probability = std::stod(value, &used);
            if (used != value.size() || !std::isfinite(probability) || probability < 0 || probability > 1)
                throw std::runtime_error("--fraction must be between 0 and 1");
            fractionSet = true;
        } else if (option == "--contacts" && !contactsSet) {
            target = unsignedValue(value, option);
            contactsSet = true;
        } else if (option == "--seed" && !seedSet) {
            seed = unsignedValue(value, option);
            seedSet = true;
        } else if (option == "--resolution" && !resolutionSet) {
            auto n = unsignedValue(value, option);
            if (!n || n > INT32_MAX) throw std::runtime_error("Invalid BP resolution");
            resolution = static_cast<int32_t>(n);
            resolutionSet = true;
        } else if ((option == "--output" || option == "-o") && output.empty()) {
            if (!isHbsPath(value)) throw std::runtime_error("Binary output must end in .hbs.gz");
            output = value;
        } else throw std::runtime_error("Unknown or repeated option: " + option);
    }
    if (fractionSet == contactsSet) throw std::runtime_error("Specify exactly one of --fraction or --contacts");
    return exportCounts(path, probability, contactsSet, target, seed, resolution, resolutionSet, output, ContactFilter::ALL);
}

static int exportCounts(const std::string& path, double probability, bool contactsSet,
                        uint64_t target, uint64_t seed, int32_t resolution, bool resolutionSet,
                        const std::string& output, ContactFilter filter) {
    std::unique_ptr<straw_v10::File> v10;
    std::unique_ptr<StrawRawReader> legacy;
    if (straw_v10::isV10(path)) v10.reset(new straw_v10::File(path));
    else legacy.reset(new StrawRawReader(path));
    auto chromosomes = v10 ? v10->chromosomes() : getChromosomesForFile(path);
    auto resolutions = v10 ? v10->resolutions() : getResolutionsForFile(path);
    if (resolutions.empty() || *std::min_element(resolutions.begin(), resolutions.end()) <= 0)
        throw std::runtime_error("Input has no valid BP resolutions");
    if (!resolutionSet) resolution = *std::min_element(resolutions.begin(), resolutions.end());
    if (std::find(resolutions.begin(), resolutions.end(), resolution) == resolutions.end())
        throw std::runtime_error("Requested BP resolution is not available");
    chromosomes.erase(std::remove_if(chromosomes.begin(), chromosomes.end(), [](const chromosome& c) {
        std::string name = c.name;
        std::transform(name.begin(), name.end(), name.begin(), [](unsigned char ch) { return std::tolower(ch); });
        return name == "all";
    }), chromosomes.end());
    std::sort(chromosomes.begin(), chromosomes.end(), [](const chromosome& a, const chromosome& b) {
        return a.index < b.index;
    });

    // Pair order is contiguous for hic_pre / hic_v10 pre. Raw block APIs emit
    // stored cis cells once, without expanding the reflected matrix triangle.
    using Visitor = std::function<void(const chromosome&, const chromosome&, uint64_t, uint64_t, uint64_t)>;
    auto visit = [&](int32_t bp, const Visitor& consume) {
        for (size_t i = 0; i < chromosomes.size(); ++i) {
            for (size_t j = i; j < chromosomes.size(); ++j) {
                const auto& a = chromosomes[i];
                const auto& b = chromosomes[j];
                auto emit = [&](uint64_t x, uint64_t y, uint64_t n) {
                    if (i == j && x > y) std::swap(x, y);
                    consume(a, b, x * static_cast<uint64_t>(bp), y * static_cast<uint64_t>(bp), n);
                };
                if (v10) {
                    v10->raw(a.name, b.name, "BP", bp, 0, UINT32_MAX, 0, UINT32_MAX,
                        [&](const straw_v10::Record& r) {
                            emit(r.binX, r.binY, r.isScore ? integerCount(r.score) : r.count);
                        });
                } else {
                    // One reader for the whole sweep: forEachRawObservedBlock
                    // re-parsed the header and footer for each of the ~N^2/2 pairs.
                    legacy->forEachBlock(a.name, b.name, bp,
                        [&](const std::vector<contactRecord>& records) {
                            for (const auto& r : records) {
                                if (r.binX < 0 || r.binY < 0) throw std::runtime_error("Negative bin coordinate");
                                emit(r.binX, r.binY, integerCount(r.counts));
                            }
                        });
                }
            }
        }
    };

    if (contactsSet) {
        const auto coarse = *std::max_element(resolutions.begin(), resolutions.end());
        uint64_t total = 0;
        visit(coarse, [&](const chromosome&, const chromosome&, uint64_t, uint64_t, uint64_t n) {
            if (n > std::numeric_limits<uint64_t>::max() - total)
                throw std::runtime_error("Total raw contacts exceed uint64 capacity");
            total += n;
        });
        if (target > total) throw std::runtime_error("--contacts exceeds the total raw contacts (" + std::to_string(total) + ")");
        probability = total ? static_cast<double>(static_cast<long double>(target) / total) : 0;
        std::cerr << "Total raw contacts at " << coarse << " BP: " << total << '\n';
    }
    std::mt19937_64 rng(seed);
    std::unique_ptr<HbsWriter> binary;
    std::sort(resolutions.begin(), resolutions.end());
    for (;;) {
        std::cerr << "Sampling at " << resolution << " BP; probability=" << std::setprecision(17)
                  << probability << "; seed=" << seed << '\n';
        bool hasCounts = false;
        visit(resolution, [&](const chromosome& a, const chromosome& b, uint64_t x, uint64_t y, uint64_t n) {
            hasCounts = hasCounts || n != 0;
            if (filter == ContactFilter::INTER && a.index == b.index) return;
            if (filter != ContactFilter::ALL && filter != ContactFilter::INTER) {
                if (a.index != b.index) return;
                const auto distance = (x > y ? x - y : y - x) / resolution;
                if (filter == ContactFilter::INTRA_SHORT && distance >= uint64_t(5000000 / resolution)) return;
                if (filter == ContactFilter::INTRA_LONG && distance <= uint64_t(5000000 / resolution)) return;
            }
            const auto kept = sample(n, probability, rng);
            if (kept && !output.empty()) {
                if (!binary) binary.reset(new HbsWriter(path, output, resolution, chromosomes));
                binary->record(a, x, b, y, kept);
            } else if (kept) std::cout << a.name << '\t' << x << '\t' << b.name << '\t' << y << '\t' << kept << '\n';
            if (!std::cout) throw std::runtime_error("Failed writing short-format output");
        });
        auto next = std::upper_bound(resolutions.begin(), resolutions.end(), resolution);
        if (resolutionSet || hasCounts || next == resolutions.end()) break;
        // Converted V9 files can advertise an ALL-only resolution. No real
        // contacts were emitted (or sampled), so advance without consuming RNG.
        std::cerr << "No real-chromosome counts at " << resolution << " BP; trying " << *next << " BP\n";
        resolution = *next;
    }
    if (!output.empty()) {
        if (!binary) binary.reset(new HbsWriter(path, output, resolution, chromosomes));
        binary->finish();
    }
    std::cout.flush();
    if (!std::cout) throw std::runtime_error("Failed writing short-format output");
    return 0;
}

void dumpHbs(const std::string& input, const std::string& output, int32_t resolution, ContactFilter filter) {
    exportCounts(input, 1, false, 0, 1, resolution, true, output, filter);
}
