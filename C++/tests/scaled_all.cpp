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
    if (argc > 2) {
        straw_v10::File kb(argv[2]);
        assert(kb.allOverview(all));
        assert(all.version == 2 && all.syntheticLength == 170 && all.displayResolution == 17);
        assert(all.storageStride == 17 && all.coordinateScale == 1000 && all.binCount == 10);
        assert(all.sourceOffsets.at(1) == 95000 && all.chromosomeLengths.at(0) == 95000);
        sum = 0; cells = 0;
        kb.raw("All", "All", "BP", 1, 0, 170, 0, 170, [&](const straw_v10::Record &r) {
            assert(r.binX % 17 == 0 && r.binY % 17 == 0);
            sum += r.count; ++cells;
        });
        assert(cells == 2 && sum == (uint64_t(1) << 53) + 49);
    }
    std::cout << "Scaled All: layout, finest query, exact counts, collapsed level and omitted trans passed\n";
}
