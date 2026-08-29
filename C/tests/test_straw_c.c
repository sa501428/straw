#include "straw/straw.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>

static void check(straw_status_t status, straw_error_t *error) {
    if (status != STRAW_STATUS_OK) {
        fprintf(stderr, "%s: %s\n", straw_status_name(status), straw_error_message(error));
        straw_error_free(error);
        assert(status == STRAW_STATUS_OK);
    }
}

int main(int argc, char **argv) {
    assert(argc == 2);
    assert(straw_abi_version() == STRAW_ABI_VERSION);

    straw_file_t *file = NULL;
    straw_error_t *error = NULL;
    check(straw_file_open(argv[1], &file, &error), error);
    assert(file != NULL);
    assert(straw_file_format_version(file) >= 6);
    assert(straw_file_chromosome_count(file) > 2);
    assert(straw_file_resolution_count(file, "BP") > 0);

    straw_query_options_t options = {
        sizeof(options), "observed", "NONE", "1", "2", "BP", 2500000
    };
    straw_records_t *forward = NULL;
    check(straw_query_records(file, &options, &forward, &error), error);
    assert(straw_records_size(forward) > 0);

    options.first_location = "2";
    options.second_location = "1";
    straw_records_t *reverse = NULL;
    check(straw_query_records(file, &options, &reverse, &error), error);
    assert(straw_records_size(forward) == straw_records_size(reverse));
    assert(straw_records_x(forward)[0] == straw_records_y(reverse)[0]);
    assert(straw_records_y(forward)[0] == straw_records_x(reverse)[0]);
    assert(straw_records_values(forward)[0] == straw_records_values(reverse)[0]);

    straw_query_t *query = NULL;
    check(straw_query_prepare(file, "observed", "NONE", "1", "1", "BP", 2500000,
                              &query, &error), error);
    straw_region_t regions[2] = {{0, 5000000, 10000000, 15000000},
                                 {10000000, 15000000, 0, 5000000}};
    straw_batch_t *batch = NULL;
    check(straw_query_regions(query, regions, 2, &batch, &error), error);
    assert(straw_batch_region_count(batch) == 2);
    assert(straw_batch_offsets(batch)[1] > 0);
    assert(straw_batch_offsets(batch)[1] * 2 == straw_batch_offsets(batch)[2]);

    straw_chromosome_counts_t *counts = NULL;
    check(straw_file_chromosome_record_counts(file, 2500000, &counts, &error), error);
    assert(straw_chromosome_counts_size(counts) > 0);
    assert(strlen(straw_chromosome_counts_name(counts, 0)) > 0);
    assert(straw_chromosome_counts_values(counts) != NULL);
    straw_chromosome_counts_free(counts);

    straw_batch_free(batch);
    straw_query_close(query);
    straw_records_free(reverse);
    straw_records_free(forward);
    straw_file_close(file);

    file = NULL;
    error = NULL;
    assert(straw_file_open("this-file-does-not-exist.hic", &file, &error) != STRAW_STATUS_OK);
    assert(file == NULL);
    assert(error != NULL);
    straw_error_free(error);
    return 0;
}
