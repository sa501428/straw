#include "straw_v10.h"
#include <cstring>
#include <iostream>
#include <sstream>
int main(int argc, char **argv) {
    try {
        if (argc < 3)
            return 2;
        std::string path = argv[1], mode = argv[2];
        straw_v10::File file(path);
        if (mode == "raw") {
            auto cb = [](const straw_v10::Record &r) {
                std::cout << r.binX << ' ' << r.binY << ' ';
                if (r.isScore) {
                    uint32_t bits;
                    std::memcpy(&bits, &r.score, 4);
                    std::cout << "s " << bits;
                } else
                    std::cout << "c " << r.count;
                std::cout << '\n';
            };
            file.raw(argv[3], argv[4], argc > 6 ? argv[6] : "BP", std::stoi(argv[5]), 0, UINT32_MAX,
                     0, UINT32_MAX, cb);
        } else if (mode == "norm") {
            for (double v : file.normalization(argv[3], "BP", std::stoi(argv[4]), "VC"))
                std::cout << v << '\n';
        } else if (mode == "meta") {
            for (auto c : getChromosomesForFile(path))
                std::cout << c.index << ' ' << c.name << ' ' << c.length << '\n';
            for (auto r : getResolutionsForFile(path))
                std::cout << r << '\n';
            std::ifstream in(path, std::ios::binary);
            int64_t master = 0, nvi = 0, nvil = 0;
            int32_t n = 0, version = 0;
            std::string genome;
            auto cs = readHeader(in, master, genome, n, version, nvi, nvil);
            if (version != 10 || cs.size() != file.chromosomes().size())
                return 3;
        } else if (mode == "callbacks") {
            std::vector<double> norm;
            bool ok = forEachRawObservedBlockWithNorm(
                path, "chrA", 10, "VC", norm, [](const std::vector<contactRecord> &b) {
                    for (auto r : b)
                        std::cout << r.binX << ' ' << r.binY << ' ' << r.counts << '\n';
                });
            if (!ok)
                return 3;
            std::vector<StrawRegion> regions{{0, 40, 0, 40}, {0, 20, 20, 40}};
            if (!strawStreamRegions(
                    path, "chrA", 10, "VC", regions, [](size_t i, const contactRecord &r) {
                        std::cout << i << ' ' << r.binX << ' ' << r.binY << ' ' << r.counts << '\n';
                    }))
                return 4;
            std::cout << getNumRecordsForFile(path, 10, false) << '\n';
        }
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
