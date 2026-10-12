#include "straw_v10.h"
#include <cassert>
#include <iostream>
int main(int argc, char **argv) {
    straw_v10::File file(argv[1]);
    straw_v10::AllOverview all;
    assert(file.allOverview(all));
    assert(all.chromosomeIndex == 2 && all.finestResolution == 1 && all.binCount == 18);
    assert(all.effectiveResolution == 10 && all.sourceOffsets.at(1) == 10);
    uint64_t sum = 0; unsigned cells = 0;
    file.raw("All", "All", "BP", 1, 0, 18, 0, 18, [&](const straw_v10::Record &r) {
        if (r.binX == 1 && r.binY == 2) assert(r.count == (uint64_t(1) << 53) + 17);
        sum += r.count; ++cells;
    });
    assert(cells == 3 && sum == (uint64_t(1) << 53) + 49);
    unsigned trans = 0;
    file.raw("All", "chrA", "BP", 1, 0, 18, 0, 95, [&](const straw_v10::Record &) { ++trans; });
    assert(trans == 0);
    sum = 0; cells = 0;
    file.raw("All", "All", "BP", 100, 0, 1, 0, 1, [&](const straw_v10::Record &r) { sum += r.count; ++cells; });
    assert(cells == 1 && sum == (uint64_t(1) << 53) + 49);
    std::cout << "Scaled All: layout, finest query, exact counts, collapsed level and omitted trans passed\n";
}
