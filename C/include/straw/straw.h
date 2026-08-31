#ifndef LIBSTRAW_STRAW_H
#define LIBSTRAW_STRAW_H

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32) && defined(STRAW_BUILDING_LIBRARY)
#  define STRAW_API __declspec(dllexport)
#elif defined(_WIN32)
#  define STRAW_API __declspec(dllimport)
#elif defined(__GNUC__) || defined(__clang__)
#  define STRAW_API __attribute__((visibility("default")))
#else
#  define STRAW_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define STRAW_ABI_VERSION_MAJOR 1u
#define STRAW_ABI_VERSION_MINOR 0u
#define STRAW_ABI_VERSION ((STRAW_ABI_VERSION_MAJOR << 16u) | STRAW_ABI_VERSION_MINOR)

typedef enum straw_status {
    STRAW_STATUS_OK = 0,
    STRAW_STATUS_INVALID_ARGUMENT = 1,
    STRAW_STATUS_IO_ERROR = 2,
    STRAW_STATUS_UNSUPPORTED_VERSION = 3,
    STRAW_STATUS_CORRUPT_FILE = 4,
    STRAW_STATUS_NOT_FOUND = 5,
    STRAW_STATUS_UNAVAILABLE = 6,
    STRAW_STATUS_UNSUPPORTED_OPERATION = 7,
    STRAW_STATUS_ALLOCATION_FAILED = 8,
    STRAW_STATUS_INTERNAL_ERROR = 9
} straw_status_t;

typedef enum straw_raw_value_kind {
    STRAW_RAW_COUNT = 0,
    STRAW_RAW_SCORE = 1
} straw_raw_value_kind_t;

typedef struct straw_file straw_file_t;
typedef struct straw_query straw_query_t;
typedef struct straw_error straw_error_t;
typedef struct straw_records straw_records_t;
typedef struct straw_raw_records straw_raw_records_t;
typedef struct straw_dense_matrix straw_dense_matrix_t;
typedef struct straw_vector straw_vector_t;
typedef struct straw_batch straw_batch_t;
typedef struct straw_chromosome_counts straw_chromosome_counts_t;

typedef struct straw_query_options {
    size_t struct_size;
    const char *matrix_type;
    const char *normalization;
    const char *first_location;
    const char *second_location;
    const char *unit;
    int32_t resolution;
} straw_query_options_t;

typedef struct straw_region {
    int64_t x_start;
    int64_t x_end;
    int64_t y_start;
    int64_t y_end;
} straw_region_t;

STRAW_API uint32_t straw_abi_version(void);
STRAW_API const char *straw_library_version(void);
STRAW_API const char *straw_status_name(straw_status_t status);

STRAW_API straw_status_t straw_error_status(const straw_error_t *error);
STRAW_API const char *straw_error_message(const straw_error_t *error);
STRAW_API void straw_error_free(straw_error_t *error);

STRAW_API straw_status_t straw_file_open(const char *path_or_url,
                                         straw_file_t **out_file,
                                         straw_error_t **out_error);
STRAW_API void straw_file_close(straw_file_t *file);
STRAW_API const char *straw_file_path(const straw_file_t *file);
STRAW_API int32_t straw_file_format_version(const straw_file_t *file);
STRAW_API const char *straw_file_genome(const straw_file_t *file);

STRAW_API size_t straw_file_chromosome_count(const straw_file_t *file);
STRAW_API const char *straw_file_chromosome_name(const straw_file_t *file, size_t index);
STRAW_API int32_t straw_file_chromosome_index(const straw_file_t *file, size_t index);
STRAW_API int64_t straw_file_chromosome_length(const straw_file_t *file, size_t index);

STRAW_API size_t straw_file_resolution_count(const straw_file_t *file, const char *unit);
STRAW_API int32_t straw_file_resolution(const straw_file_t *file, const char *unit, size_t index);
STRAW_API size_t straw_file_normalization_count(const straw_file_t *file);
STRAW_API const char *straw_file_normalization(const straw_file_t *file, size_t index);
STRAW_API size_t straw_file_attribute_count(const straw_file_t *file);
STRAW_API const char *straw_file_attribute_key(const straw_file_t *file, size_t index);
STRAW_API const char *straw_file_attribute_value(const straw_file_t *file, size_t index);

STRAW_API straw_status_t straw_query_records(const straw_file_t *file,
                                             const straw_query_options_t *options,
                                             straw_records_t **out_records,
                                             straw_error_t **out_error);
STRAW_API straw_status_t straw_query_records_simple(const straw_file_t *file,
                                                    const char *matrix_type,
                                                    const char *normalization,
                                                    const char *first_location,
                                                    const char *second_location,
                                                    const char *unit,
                                                    int32_t resolution,
                                                    straw_records_t **out_records,
                                                    straw_error_t **out_error);
STRAW_API size_t straw_records_size(const straw_records_t *records);
STRAW_API const int64_t *straw_records_x(const straw_records_t *records);
STRAW_API const int64_t *straw_records_y(const straw_records_t *records);
STRAW_API const float *straw_records_values(const straw_records_t *records);
STRAW_API void straw_records_free(straw_records_t *records);

STRAW_API straw_status_t straw_query_dense(const straw_file_t *file,
                                           const straw_query_options_t *options,
                                           straw_dense_matrix_t **out_matrix,
                                           straw_error_t **out_error);
STRAW_API size_t straw_dense_rows(const straw_dense_matrix_t *matrix);
STRAW_API size_t straw_dense_columns(const straw_dense_matrix_t *matrix);
STRAW_API const float *straw_dense_values_row_major(const straw_dense_matrix_t *matrix);
STRAW_API void straw_dense_matrix_free(straw_dense_matrix_t *matrix);

STRAW_API straw_status_t straw_file_normalization_vector(const straw_file_t *file,
                                                         const char *chromosome,
                                                         const char *unit,
                                                         int32_t resolution,
                                                         const char *normalization,
                                                         straw_vector_t **out_vector,
                                                         straw_error_t **out_error);
STRAW_API straw_status_t straw_file_expected_vector(const straw_file_t *file,
                                                    const char *chromosome,
                                                    const char *unit,
                                                    int32_t resolution,
                                                    const char *normalization,
                                                    straw_vector_t **out_vector,
                                                    straw_error_t **out_error);
STRAW_API size_t straw_vector_size(const straw_vector_t *vector);
STRAW_API const double *straw_vector_values(const straw_vector_t *vector);
STRAW_API void straw_vector_free(straw_vector_t *vector);

STRAW_API straw_status_t straw_query_raw(const straw_file_t *file,
                                         const char *first_location,
                                         const char *second_location,
                                         const char *unit,
                                         int32_t resolution,
                                         straw_raw_records_t **out_records,
                                         straw_error_t **out_error);
STRAW_API size_t straw_raw_records_size(const straw_raw_records_t *records);
STRAW_API const uint64_t *straw_raw_records_x(const straw_raw_records_t *records);
STRAW_API const uint64_t *straw_raw_records_y(const straw_raw_records_t *records);
STRAW_API const uint64_t *straw_raw_records_counts(const straw_raw_records_t *records);
STRAW_API const float *straw_raw_records_scores(const straw_raw_records_t *records);
STRAW_API const uint8_t *straw_raw_records_kinds(const straw_raw_records_t *records);
STRAW_API void straw_raw_records_free(straw_raw_records_t *records);

STRAW_API straw_status_t straw_query_prepare(const straw_file_t *file,
                                             const char *matrix_type,
                                             const char *normalization,
                                             const char *first_chromosome,
                                             const char *second_chromosome,
                                             const char *unit,
                                             int32_t resolution,
                                             straw_query_t **out_query,
                                             straw_error_t **out_error);
STRAW_API void straw_query_close(straw_query_t *query);
STRAW_API straw_status_t straw_query_window(const straw_query_t *query,
                                            const straw_region_t *region,
                                            straw_records_t **out_records,
                                            straw_error_t **out_error);
STRAW_API straw_status_t straw_query_regions(const straw_query_t *query,
                                             const straw_region_t *regions,
                                             size_t region_count,
                                             straw_batch_t **out_batch,
                                             straw_error_t **out_error);
STRAW_API size_t straw_batch_region_count(const straw_batch_t *batch);
STRAW_API size_t straw_batch_record_count(const straw_batch_t *batch);
STRAW_API const size_t *straw_batch_offsets(const straw_batch_t *batch);
STRAW_API const int64_t *straw_batch_x(const straw_batch_t *batch);
STRAW_API const int64_t *straw_batch_y(const straw_batch_t *batch);
STRAW_API const float *straw_batch_values(const straw_batch_t *batch);
STRAW_API void straw_batch_free(straw_batch_t *batch);

STRAW_API straw_status_t straw_file_record_count(const straw_file_t *file,
                                                 int32_t resolution,
                                                 int inter_only,
                                                 uint64_t *out_count,
                                                 straw_error_t **out_error);

STRAW_API straw_status_t straw_file_chromosome_record_counts(const straw_file_t *file,
                                                              int32_t resolution,
                                                              straw_chromosome_counts_t **out_counts,
                                                              straw_error_t **out_error);
STRAW_API size_t straw_chromosome_counts_size(const straw_chromosome_counts_t *counts);
STRAW_API const char *straw_chromosome_counts_name(const straw_chromosome_counts_t *counts, size_t index);
STRAW_API const uint64_t *straw_chromosome_counts_values(const straw_chromosome_counts_t *counts);
STRAW_API void straw_chromosome_counts_free(straw_chromosome_counts_t *counts);

#ifdef __cplusplus
}
#endif

#endif
