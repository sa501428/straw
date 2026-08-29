#include "straw/straw.h"

#include <stdio.h>

int main(int argc, char **argv) {
    if (argc != 6) {
        fprintf(stderr, "usage: %s file.hic chr1 chr2 unit resolution\n", argv[0]);
        return 2;
    }
    straw_file_t *file = NULL;
    straw_error_t *error = NULL;
    straw_status_t status = straw_file_open(argv[1], &file, &error);
    if (status != STRAW_STATUS_OK) {
        fprintf(stderr, "%s\n", straw_error_message(error));
        straw_error_free(error);
        return 1;
    }
    straw_query_options_t options = {
        sizeof(options), "observed", "NONE", argv[2], argv[3], argv[4], 0
    };
    if (sscanf(argv[5], "%d", &options.resolution) != 1) return 2;
    straw_records_t *records = NULL;
    status = straw_query_records(file, &options, &records, &error);
    if (status != STRAW_STATUS_OK) {
        fprintf(stderr, "%s\n", straw_error_message(error));
        straw_error_free(error);
        straw_file_close(file);
        return 1;
    }
    for (size_t i = 0; i < straw_records_size(records); ++i)
        printf("%lld\t%lld\t%.9g\n", (long long)straw_records_x(records)[i],
               (long long)straw_records_y(records)[i], straw_records_values(records)[i]);
    straw_records_free(records);
    straw_file_close(file);
    return 0;
}
