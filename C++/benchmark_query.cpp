#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include <cstring>

#include "straw.h"
#include "norm_vectors.h"
#include "straw_v10.h"
#include "file_source.h"
#include "legacy_norm_index.h"

namespace {
std::string location(const char *chromosome, const char *start, const char *end) {
    return std::string(chromosome) + ":" + start + ":" + end;
}
}

int main(int argc, char **argv) {
    if (argc == 3 && std::string(argv[1]) == "--inspect") {
        try {
            const auto info = getFileInfo(argv[2]);
            std::cout << "version\t" << info.version << '\n';
            for (const auto &c : info.chromosomes)
                std::cout << "chromosome\t" << c.name << '\t' << c.length << '\n';
            for (auto r : info.bpResolutions) std::cout << "resolution\t" << r << '\n';
            for (const auto &n : info.normalizations) std::cout << "normalization\t" << n << '\n';
            if (info.version == 10) {
                straw_v10::File file(argv[2]);
                for (const auto &e : file.normalizationEntries())
                    std::cout << "vector\t" << e.normalization << '\t' << e.chromosome
                              << '\t' << e.unit << '\t' << e.resolution << '\n';
            } else {
                straw_file::Source source(argv[2]);
                auto header = source.read(0, 4096);
                straw_v10::Cursor cursor(header);
                cursor.at = 16;
                cursor.str();
                const auto position = cursor.wide(), length = cursor.wide();
                if (!position || !length) throw std::runtime_error("inspection requires direct v9 NVI locator");
                auto bytes = source.read(position, length);
                for (const auto &e : straw_v9::normalizationIndex(bytes.data(), bytes.size())) {
                    const auto c = std::find_if(info.chromosomes.begin(), info.chromosomes.end(),
                                               [&](const chromosome &c) { return c.index == e.chromosome; });
                    if (c == info.chromosomes.end()) throw std::runtime_error("invalid NVI chromosome");
                    std::cout << "vector\t" << e.normalization << '\t' << c->name
                              << '\t' << e.unit << '\t' << e.resolution << '\n';
                }
            }
            return 0;
        } catch (const std::exception &error) {
            std::cerr << error.what() << '\n';
            return 1;
        }
    }
    if (argc != 9 && argc != 10) {
        std::cerr << "usage: straw_benchmark_query URL NORM CHR X_START X_END Y_START Y_END RESOLUTION [RECORDS.bin]\n";
        return 2;
    }
    try {
        const auto resolution = static_cast<int32_t>(std::stol(argv[8]));
        const std::string first = location(argv[3], argv[4], argv[5]);
        const std::string second = location(argv[3], argv[6], argv[7]);
        const auto begin = std::chrono::steady_clock::now();
        const auto records = straw("observed", argv[2], argv[1], first, second, "BP", resolution);
        const auto end = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed = end - begin;
        // Optional audit payload: little-endian int32 binX, int32 binY, float64
        // value. Export/validation are deliberately outside the query timer.
        if (argc == 10) {
            std::ofstream output(argv[9], std::ios::binary);
            if (!output) throw std::runtime_error("cannot create record audit payload");
            auto little = [&](uint64_t value, unsigned bytes) {
                for (unsigned i = 0; i < bytes; ++i) output.put(char((value >> (8 * i)) & 255));
            };
            for (const auto &record : records) {
                little(uint32_t(record.binX / resolution), 4);
                little(uint32_t(record.binY / resolution), 4);
                double value = record.counts;
                uint64_t bits;
                std::memcpy(&bits, &value, sizeof(bits));
                little(bits, 8);
            }
            output.close();
            if (!output) throw std::runtime_error("failed writing record audit payload");
        }
        std::cout << std::setprecision(17) << elapsed.count() << '\t' << records.size() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
