#include "straw/straw.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void check(straw_status_t status, straw_error_t *error) {
    if (status != STRAW_STATUS_OK) {
        fprintf(stderr, "%s\n", straw_error_message(error));
        exit(1);
    }
}

static void close_to(double expected, double actual) {
    if (!isfinite(actual) || fabs(expected - actual) > 1e-5) {
        fprintf(stderr, "expected %.9g, got %.9g\n", expected, actual);
        exit(1);
    }
}

int main(int argc, char **argv) {
    if (argc != 2) return 2;
    straw_file_t *file = NULL;
    straw_error_t *error = NULL;
    straw_status_t status = straw_file_open(argv[1], &file, &error);
    check(status, error);
    for (int resolution = 10; resolution <= 20; resolution += 10) {
        double mean = resolution == 10 ? 7.0 / 56 : 7.0 / 16;
        double product = resolution == 10 ? 6 : 24;
        for (int normalized = 0; normalized < 2; ++normalized) {
            const char *norm = normalized ? "VC" : "NONE";
            for (int expected = 0; expected < 2; ++expected) {
                straw_records_t *records = NULL;
                status = straw_query_records_simple(file, expected ? "expected" : "oe", norm,
                    "chrB", "chrA", "BP", resolution, &records, &error);
                check(status, error);
                if (straw_records_size(records) != 3) return 1;
                const float *values = straw_records_values(records);
                for (int i = 0; i < 3; ++i) {
                    double count = i == 2 ? 5 : 1;
                    close_to(expected ? mean : count / mean / (normalized ? product : 1), values[i]);
                }
                straw_records_free(records);
            }
            straw_query_t *query = NULL;
            status = straw_query_prepare(file, "oe", norm, "chrB", "chrA", "BP", resolution,
                &query, &error);
            check(status, error);
            straw_region_t region = {0, 10, 0, 10};
            straw_records_t *window = NULL;
            status = straw_query_window(query, &region, &window, &error);
            check(status, error);
            if (straw_records_size(window) != 1) return 1;
            close_to(1 / mean / (normalized ? product : 1), straw_records_values(window)[0]);
            straw_records_free(window);
            straw_query_close(query);
        }
    }
    straw_file_close(file);
    return 0;
}
