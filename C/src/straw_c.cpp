#include "straw/straw.h"

#include "straw.h"
#include "straw_v10.h"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct straw_error {
    straw_status_t status;
    std::string message;
};

struct straw_file {
    std::string path;
    int32_t version = 0;
    std::string genome;
    std::vector<chromosome> chromosomes;
    std::vector<int32_t> bp_resolutions;
    std::vector<int32_t> frag_resolutions;
    std::vector<std::string> normalizations;
    std::vector<std::pair<std::string, std::string>> attributes;
};

struct straw_query {
    std::unique_ptr<StrawPreparedQuery> native;
};

struct straw_records {
    std::vector<int64_t> x, y;
    std::vector<float> values;
};

struct straw_raw_records {
    std::vector<uint64_t> x, y, counts;
    std::vector<float> scores;
    std::vector<uint8_t> kinds;
};

struct straw_dense_matrix {
    size_t rows = 0, columns = 0;
    std::vector<float> values;
};

struct straw_vector {
    std::vector<double> values;
};

struct straw_batch {
    size_t regions = 0;
    std::vector<size_t> offsets;
    std::vector<int64_t> x, y;
    std::vector<float> values;
};

struct straw_chromosome_counts {
    std::vector<std::string> names;
    std::vector<uint64_t> values;
};

namespace {

void clear_error(straw_error_t **out_error) {
    if (out_error) *out_error = nullptr;
}

straw_status_t fail(straw_status_t status, const std::string &message,
                    straw_error_t **out_error) noexcept {
    if (out_error) {
        try {
            *out_error = new straw_error{status, message};
        } catch (...) {
            *out_error = nullptr;
        }
    }
    return status;
}

straw_status_t classify(const std::exception &error, straw_error_t **out_error) noexcept {
    straw_status_t status = STRAW_STATUS_INTERNAL_ERROR;
    if (const auto *typed = dynamic_cast<const StrawException *>(&error)) {
        switch (typed->code()) {
            case StrawErrorCode::InvalidArgument: status = STRAW_STATUS_INVALID_ARGUMENT; break;
            case StrawErrorCode::Io: status = STRAW_STATUS_IO_ERROR; break;
            case StrawErrorCode::UnsupportedVersion: status = STRAW_STATUS_UNSUPPORTED_VERSION; break;
            case StrawErrorCode::CorruptFile: status = STRAW_STATUS_CORRUPT_FILE; break;
            case StrawErrorCode::NotFound: status = STRAW_STATUS_NOT_FOUND; break;
            case StrawErrorCode::Unavailable: status = STRAW_STATUS_UNAVAILABLE; break;
            case StrawErrorCode::UnsupportedOperation: status = STRAW_STATUS_UNSUPPORTED_OPERATION; break;
            case StrawErrorCode::Internal: status = STRAW_STATUS_INTERNAL_ERROR; break;
        }
    } else if (dynamic_cast<const std::invalid_argument *>(&error)) {
        status = STRAW_STATUS_INVALID_ARGUMENT;
    }
    return fail(status, error.what(), out_error);
}

straw_status_t allocation_failure(straw_error_t **out_error) noexcept {
    return fail(STRAW_STATUS_ALLOCATION_FAILED, "native allocation failed", out_error);
}

bool valid_text(const char *value) { return value && *value; }

std::string chromosome_from_location(const std::string &value) {
    return value.substr(0, value.find(':'));
}

bool has_chromosome(const straw_file &file, const std::string &name) {
    return std::any_of(file.chromosomes.begin(), file.chromosomes.end(),
                       [&](const chromosome &value) { return value.name == name; });
}

void validate_query_capabilities(const straw_file &file, const std::string &matrix_type,
                                 const std::string &normalization, const std::string &first,
                                 const std::string &second, const std::string &unit,
                                 int32_t resolution) {
    if (matrix_type != "observed" && matrix_type != "oe" && matrix_type != "expected")
        throw StrawException(StrawErrorCode::InvalidArgument,
                             "matrix type must be observed, oe, or expected");
    if (unit != "BP" && unit != "FRAG")
        throw StrawException(StrawErrorCode::InvalidArgument, "unit must be BP or FRAG");
    if (!has_chromosome(file, chromosome_from_location(first)) ||
        !has_chromosome(file, chromosome_from_location(second)))
        throw StrawException(StrawErrorCode::NotFound,
                             "requested chromosome is not present in the file");
    const auto &resolutions = unit == "BP" ? file.bp_resolutions : file.frag_resolutions;
    if (std::find(resolutions.begin(), resolutions.end(), resolution) == resolutions.end())
        throw StrawException(StrawErrorCode::Unavailable,
                             "requested resolution is not available");
    if (file.version == 10 && normalization != "NONE" &&
        std::find(file.normalizations.begin(), file.normalizations.end(), normalization) ==
            file.normalizations.end())
        throw StrawException(StrawErrorCode::Unavailable,
                             "requested normalization is not available");
}

template <typename Function>
auto native_file_call(const straw_file &file, Function function) -> decltype(function()) {
    try {
        return function();
    } catch (const StrawException &) {
        throw;
    } catch (const std::exception &error) {
        const bool remote = file.path.compare(0, 4, "http") == 0;
        throw StrawException(remote ? StrawErrorCode::Io : StrawErrorCode::CorruptFile,
                             error.what());
    }
}

straw_status_t validate_options(const straw_query_options_t *options,
                                straw_error_t **out_error) {
    const size_t minimum = offsetof(straw_query_options_t, resolution) + sizeof(options->resolution);
    if (!options || options->struct_size < minimum || !valid_text(options->matrix_type) ||
        !valid_text(options->normalization) || !valid_text(options->first_location) ||
        !valid_text(options->second_location) || !valid_text(options->unit) ||
        options->resolution <= 0) {
        return fail(STRAW_STATUS_INVALID_ARGUMENT, "invalid query options", out_error);
    }
    return STRAW_STATUS_OK;
}

std::string location(const std::string &chromosome, int64_t start, int64_t end) {
    return chromosome + ":" + std::to_string(start) + ":" + std::to_string(end);
}

straw_records *read_records(const std::string &path, const std::string &matrix_type,
                            const std::string &normalization, const std::string &first,
                            const std::string &second, const std::string &unit,
                            int32_t resolution) {
    std::vector<contactRecord> source = straw(matrix_type, normalization, path, first, second,
                                               unit, resolution);
    straw_records *result = new straw_records;
    result->x.reserve(source.size());
    result->y.reserve(source.size());
    result->values.reserve(source.size());
    for (const auto &record : source) {
        result->x.push_back(record.binX);
        result->y.push_back(record.binY);
        result->values.push_back(record.counts);
    }
    return result;
}

straw_records *read_records(StrawPreparedQuery &query, int64_t x_start, int64_t x_end,
                            int64_t y_start, int64_t y_end) {
    std::unique_ptr<straw_records> result(new straw_records);
    query.streamWindow(x_start, x_end, y_start, y_end, [&](const contactRecord &record) {
        result->x.push_back(record.binX);
        result->y.push_back(record.binY);
        result->values.push_back(record.counts);
    });
    return result.release();
}

template <typename Function>
straw_status_t guarded(straw_error_t **out_error, Function function) noexcept {
    clear_error(out_error);
    try {
        function();
        return STRAW_STATUS_OK;
    } catch (const std::bad_alloc &) {
        return allocation_failure(out_error);
    } catch (const std::exception &error) {
        return classify(error, out_error);
    } catch (...) {
        return fail(STRAW_STATUS_INTERNAL_ERROR, "unknown native error", out_error);
    }
}

} // namespace

extern "C" {

uint32_t straw_abi_version(void) { return STRAW_ABI_VERSION; }
const char *straw_library_version(void) { return "1.0.0"; }

const char *straw_status_name(straw_status_t status) {
    switch (status) {
        case STRAW_STATUS_OK: return "ok";
        case STRAW_STATUS_INVALID_ARGUMENT: return "invalid_argument";
        case STRAW_STATUS_IO_ERROR: return "io_error";
        case STRAW_STATUS_UNSUPPORTED_VERSION: return "unsupported_version";
        case STRAW_STATUS_CORRUPT_FILE: return "corrupt_file";
        case STRAW_STATUS_NOT_FOUND: return "not_found";
        case STRAW_STATUS_UNAVAILABLE: return "unavailable";
        case STRAW_STATUS_UNSUPPORTED_OPERATION: return "unsupported_operation";
        case STRAW_STATUS_ALLOCATION_FAILED: return "allocation_failed";
        case STRAW_STATUS_INTERNAL_ERROR: return "internal_error";
    }
    return "unknown";
}

straw_status_t straw_error_status(const straw_error_t *error) {
    return error ? error->status : STRAW_STATUS_OK;
}
const char *straw_error_message(const straw_error_t *error) {
    return error ? error->message.c_str() : "";
}
void straw_error_free(straw_error_t *error) { delete error; }

straw_status_t straw_file_open(const char *path_or_url, straw_file_t **out_file,
                               straw_error_t **out_error) {
    if (out_file) *out_file = nullptr;
    if (!out_file || !valid_text(path_or_url)) {
        clear_error(out_error);
        return fail(STRAW_STATUS_INVALID_ARGUMENT, "path and output file are required", out_error);
    }
    return guarded(out_error, [&] {
        straw_file *file = new straw_file;
        try {
            file->path = path_or_url;
            const bool remote = file->path.compare(0, 4, "http") == 0;
            if (!remote) {
                std::ifstream input(file->path, std::ios::binary);
                if (!input)
                    throw StrawException(StrawErrorCode::Io,
                                         "file cannot be opened for reading");
            }
            file->version = getVersionForFile(file->path);
            if (file->version < 6 || file->version > 10) {
                throw StrawException(StrawErrorCode::UnsupportedVersion,
                                     "unsupported .hic version");
            }
            file->genome = getGenomeForFile(file->path);
            file->chromosomes = getChromosomesForFile(file->path);
            file->bp_resolutions = getResolutionsForFile(file->path, "BP");
            file->frag_resolutions = getResolutionsForFile(file->path, "FRAG");
            file->normalizations = getNormalizationsForFile(file->path);
            file->attributes = getAttributesForFile(file->path);
            if (file->chromosomes.empty()) {
                throw StrawException(StrawErrorCode::CorruptFile,
                                     "invalid or empty .hic chromosome table");
            }
        } catch (const StrawException &) {
            delete file;
            throw;
        } catch (const std::exception &error) {
            const bool remote = file->path.compare(0, 4, "http") == 0;
            delete file;
            throw StrawException(remote ? StrawErrorCode::Io : StrawErrorCode::CorruptFile,
                                 error.what());
        } catch (...) {
            delete file;
            throw;
        }
        *out_file = file;
    });
}

void straw_file_close(straw_file_t *file) { delete file; }
const char *straw_file_path(const straw_file_t *file) { return file ? file->path.c_str() : ""; }
int32_t straw_file_format_version(const straw_file_t *file) { return file ? file->version : 0; }
const char *straw_file_genome(const straw_file_t *file) { return file ? file->genome.c_str() : ""; }

size_t straw_file_chromosome_count(const straw_file_t *file) {
    return file ? file->chromosomes.size() : 0;
}
const char *straw_file_chromosome_name(const straw_file_t *file, size_t index) {
    return file && index < file->chromosomes.size() ? file->chromosomes[index].name.c_str() : "";
}
int32_t straw_file_chromosome_index(const straw_file_t *file, size_t index) {
    return file && index < file->chromosomes.size() ? file->chromosomes[index].index : -1;
}
int64_t straw_file_chromosome_length(const straw_file_t *file, size_t index) {
    return file && index < file->chromosomes.size() ? file->chromosomes[index].length : -1;
}

size_t straw_file_resolution_count(const straw_file_t *file, const char *unit) {
    if (!file || !unit) return 0;
    if (std::string(unit) == "BP") return file->bp_resolutions.size();
    if (std::string(unit) == "FRAG") return file->frag_resolutions.size();
    return 0;
}
int32_t straw_file_resolution(const straw_file_t *file, const char *unit, size_t index) {
    if (!file || !unit) return 0;
    const auto &values = std::string(unit) == "FRAG" ? file->frag_resolutions : file->bp_resolutions;
    return index < values.size() ? values[index] : 0;
}
size_t straw_file_normalization_count(const straw_file_t *file) {
    return file ? file->normalizations.size() : 0;
}
const char *straw_file_normalization(const straw_file_t *file, size_t index) {
    return file && index < file->normalizations.size() ? file->normalizations[index].c_str() : "";
}
size_t straw_file_attribute_count(const straw_file_t *file) {
    return file ? file->attributes.size() : 0;
}
const char *straw_file_attribute_key(const straw_file_t *file, size_t index) {
    return file && index < file->attributes.size() ? file->attributes[index].first.c_str() : "";
}
const char *straw_file_attribute_value(const straw_file_t *file, size_t index) {
    return file && index < file->attributes.size() ? file->attributes[index].second.c_str() : "";
}

straw_status_t straw_query_records(const straw_file_t *file, const straw_query_options_t *options,
                                   straw_records_t **out_records, straw_error_t **out_error) {
    if (out_records) *out_records = nullptr;
    clear_error(out_error);
    if (!file || !out_records) return fail(STRAW_STATUS_INVALID_ARGUMENT, "file and output are required", out_error);
    straw_status_t status = validate_options(options, out_error);
    if (status != STRAW_STATUS_OK) return status;
    return guarded(out_error, [&] {
        validate_query_capabilities(*file, options->matrix_type, options->normalization,
                                    options->first_location, options->second_location,
                                    options->unit, options->resolution);
        *out_records = native_file_call(*file, [&] {
            return read_records(file->path, options->matrix_type, options->normalization,
                                options->first_location, options->second_location, options->unit,
                                options->resolution);
        });
    });
}

straw_status_t straw_query_records_simple(const straw_file_t *file, const char *matrix_type,
                                          const char *normalization, const char *first_location,
                                          const char *second_location, const char *unit,
                                          int32_t resolution, straw_records_t **out_records,
                                          straw_error_t **out_error) {
    straw_query_options_t options = {sizeof(straw_query_options_t), matrix_type, normalization,
                                     first_location, second_location, unit, resolution};
    return straw_query_records(file, &options, out_records, out_error);
}

size_t straw_records_size(const straw_records_t *records) { return records ? records->x.size() : 0; }
const int64_t *straw_records_x(const straw_records_t *records) { return records && !records->x.empty() ? records->x.data() : nullptr; }
const int64_t *straw_records_y(const straw_records_t *records) { return records && !records->y.empty() ? records->y.data() : nullptr; }
const float *straw_records_values(const straw_records_t *records) { return records && !records->values.empty() ? records->values.data() : nullptr; }
void straw_records_free(straw_records_t *records) { delete records; }

straw_status_t straw_query_dense(const straw_file_t *file, const straw_query_options_t *options,
                                 straw_dense_matrix_t **out_matrix, straw_error_t **out_error) {
    if (out_matrix) *out_matrix = nullptr;
    clear_error(out_error);
    if (!file || !out_matrix) return fail(STRAW_STATUS_INVALID_ARGUMENT, "file and output are required", out_error);
    straw_status_t status = validate_options(options, out_error);
    if (status != STRAW_STATUS_OK) return status;
    return guarded(out_error, [&] {
        validate_query_capabilities(*file, options->matrix_type, options->normalization,
                                    options->first_location, options->second_location,
                                    options->unit, options->resolution);
        auto source = native_file_call(*file, [&] {
            return strawAsMatrix(options->matrix_type, options->normalization, file->path,
                                 options->first_location, options->second_location, options->unit,
                                 options->resolution);
        });
        straw_dense_matrix *matrix = new straw_dense_matrix;
        matrix->rows = source.size();
        matrix->columns = source.empty() ? 0 : source.front().size();
        if (matrix->columns && matrix->rows > std::numeric_limits<size_t>::max() / matrix->columns) {
            delete matrix;
            throw std::bad_alloc();
        }
        matrix->values.reserve(matrix->rows * matrix->columns);
        for (const auto &row : source) {
            if (row.size() != matrix->columns) {
                delete matrix;
                throw StrawException(StrawErrorCode::CorruptFile,
                                     "non-rectangular native dense result");
            }
            matrix->values.insert(matrix->values.end(), row.begin(), row.end());
        }
        *out_matrix = matrix;
    });
}

size_t straw_dense_rows(const straw_dense_matrix_t *matrix) { return matrix ? matrix->rows : 0; }
size_t straw_dense_columns(const straw_dense_matrix_t *matrix) { return matrix ? matrix->columns : 0; }
const float *straw_dense_values_row_major(const straw_dense_matrix_t *matrix) { return matrix && !matrix->values.empty() ? matrix->values.data() : nullptr; }
void straw_dense_matrix_free(straw_dense_matrix_t *matrix) { delete matrix; }

straw_status_t straw_file_normalization_vector(const straw_file_t *file, const char *chromosome,
                                               const char *unit, int32_t resolution,
                                               const char *normalization, straw_vector_t **out_vector,
                                               straw_error_t **out_error) {
    if (out_vector) *out_vector = nullptr;
    clear_error(out_error);
    if (!file || !valid_text(chromosome) || !valid_text(unit) || !valid_text(normalization) ||
        resolution <= 0 || !out_vector)
        return fail(STRAW_STATUS_INVALID_ARGUMENT, "invalid vector arguments", out_error);
    return guarded(out_error, [&] {
        validate_query_capabilities(*file, "observed", normalization, chromosome, chromosome,
                                    unit, resolution);
        std::unique_ptr<straw_vector> result(new straw_vector);
        bool available = native_file_call(*file, [&] {
            return getNormalizationVectorForFile(file->path, chromosome, resolution,
                                                 normalization, result->values, unit);
        });
        if (!available) {
            throw StrawException(StrawErrorCode::Unavailable,
                                 "normalization vector is not available");
        }
        *out_vector = result.release();
    });
}

straw_status_t straw_file_expected_vector(const straw_file_t *file, const char *chromosome,
                                          const char *unit, int32_t resolution,
                                          const char *normalization, straw_vector_t **out_vector,
                                          straw_error_t **out_error) {
    if (out_vector) *out_vector = nullptr;
    clear_error(out_error);
    if (!file || !valid_text(chromosome) || !valid_text(unit) || !valid_text(normalization) ||
        resolution <= 0 || !out_vector)
        return fail(STRAW_STATUS_INVALID_ARGUMENT, "invalid vector arguments", out_error);
    return guarded(out_error, [&] {
        validate_query_capabilities(*file, "expected", normalization, chromosome, chromosome,
                                    unit, resolution);
        std::unique_ptr<straw_vector> result(new straw_vector);
        bool available = native_file_call(*file, [&] {
            return getExpectedVectorForFile(file->path, chromosome, resolution,
                                            normalization, result->values, unit);
        });
        if (!available) {
            throw StrawException(StrawErrorCode::Unavailable,
                                 "expected vector is not available");
        }
        *out_vector = result.release();
    });
}

size_t straw_vector_size(const straw_vector_t *vector) { return vector ? vector->values.size() : 0; }
const double *straw_vector_values(const straw_vector_t *vector) { return vector && !vector->values.empty() ? vector->values.data() : nullptr; }
void straw_vector_free(straw_vector_t *vector) { delete vector; }

straw_status_t straw_query_raw(const straw_file_t *file, const char *first_location,
                               const char *second_location, const char *unit, int32_t resolution,
                               straw_raw_records_t **out_records, straw_error_t **out_error) {
    if (out_records) *out_records = nullptr;
    clear_error(out_error);
    if (!file || !valid_text(first_location) || !valid_text(second_location) || !valid_text(unit) ||
        resolution <= 0 || !out_records)
        return fail(STRAW_STATUS_INVALID_ARGUMENT, "invalid raw query arguments", out_error);
    if (file->version != 10)
        return fail(STRAW_STATUS_UNSUPPORTED_OPERATION,
                    "exact raw records are available only for V10 files", out_error);
    return guarded(out_error, [&] {
        validate_query_capabilities(*file, "observed", "NONE", first_location,
                                    second_location, unit, resolution);
        straw_raw_records *result = new straw_raw_records;
        try {
            native_file_call(*file, [&] {
                straw_v10::File native(file->path);
                native.streamRaw(first_location, second_location, unit, resolution,
                             [&](const straw_v10::Record &record) {
                result->x.push_back(record.binX);
                result->y.push_back(record.binY);
                result->counts.push_back(record.isScore ? 0 : record.count);
                result->scores.push_back(record.isScore ? record.score : 0.0f);
                result->kinds.push_back(record.isScore ? STRAW_RAW_SCORE : STRAW_RAW_COUNT);
            });
                return true;
            });
        } catch (...) {
            delete result;
            throw;
        }
        *out_records = result;
    });
}

size_t straw_raw_records_size(const straw_raw_records_t *records) { return records ? records->x.size() : 0; }
const uint64_t *straw_raw_records_x(const straw_raw_records_t *records) { return records && !records->x.empty() ? records->x.data() : nullptr; }
const uint64_t *straw_raw_records_y(const straw_raw_records_t *records) { return records && !records->y.empty() ? records->y.data() : nullptr; }
const uint64_t *straw_raw_records_counts(const straw_raw_records_t *records) { return records && !records->counts.empty() ? records->counts.data() : nullptr; }
const float *straw_raw_records_scores(const straw_raw_records_t *records) { return records && !records->scores.empty() ? records->scores.data() : nullptr; }
const uint8_t *straw_raw_records_kinds(const straw_raw_records_t *records) { return records && !records->kinds.empty() ? records->kinds.data() : nullptr; }
void straw_raw_records_free(straw_raw_records_t *records) { delete records; }

straw_status_t straw_query_prepare(const straw_file_t *file, const char *matrix_type,
                                   const char *normalization, const char *first_chromosome,
                                   const char *second_chromosome, const char *unit,
                                   int32_t resolution, straw_query_t **out_query,
                                   straw_error_t **out_error) {
    if (out_query) *out_query = nullptr;
    clear_error(out_error);
    if (!file || !valid_text(matrix_type) || !valid_text(normalization) ||
        !valid_text(first_chromosome) || !valid_text(second_chromosome) || !valid_text(unit) ||
        resolution <= 0 || !out_query)
        return fail(STRAW_STATUS_INVALID_ARGUMENT, "invalid prepared query arguments", out_error);
    return guarded(out_error, [&] {
        std::unique_ptr<straw_query> query(new straw_query);
        query->native.reset(new StrawPreparedQuery(file->path, matrix_type, normalization,
                                                   first_chromosome, second_chromosome,
                                                   unit, resolution));
        *out_query = query.release();
    });
}

void straw_query_close(straw_query_t *query) { delete query; }

straw_status_t straw_query_window(const straw_query_t *query, const straw_region_t *region,
                                  straw_records_t **out_records, straw_error_t **out_error) {
    if (out_records) *out_records = nullptr;
    clear_error(out_error);
    if (!query || !region || !out_records || region->x_start < 0 || region->y_start < 0 ||
        region->x_end < region->x_start || region->y_end < region->y_start)
        return fail(STRAW_STATUS_INVALID_ARGUMENT, "invalid query window", out_error);
    return guarded(out_error, [&] {
        *out_records = read_records(*query->native, region->x_start, region->x_end,
                                    region->y_start, region->y_end);
    });
}

straw_status_t straw_query_regions(const straw_query_t *query, const straw_region_t *regions,
                                   size_t region_count, straw_batch_t **out_batch,
                                   straw_error_t **out_error) {
    if (out_batch) *out_batch = nullptr;
    clear_error(out_error);
    if (!query || (!regions && region_count) || !out_batch)
        return fail(STRAW_STATUS_INVALID_ARGUMENT, "invalid batch query arguments", out_error);
    return guarded(out_error, [&] {
        straw_batch *batch = new straw_batch;
        try {
            batch->regions = region_count;
            batch->offsets.reserve(region_count + 1);
            batch->offsets.push_back(0);
            for (size_t i = 0; i < region_count; ++i) {
                const auto &region = regions[i];
                if (region.x_start < 0 || region.y_start < 0 || region.x_end < region.x_start ||
                    region.y_end < region.y_start)
                    throw std::invalid_argument("invalid region in batch query");
                straw_records *records = read_records(*query->native,
                    region.x_start, region.x_end, region.y_start, region.y_end);
                batch->x.insert(batch->x.end(), records->x.begin(), records->x.end());
                batch->y.insert(batch->y.end(), records->y.begin(), records->y.end());
                batch->values.insert(batch->values.end(), records->values.begin(), records->values.end());
                delete records;
                batch->offsets.push_back(batch->x.size());
            }
        } catch (...) {
            delete batch;
            throw;
        }
        *out_batch = batch;
    });
}

size_t straw_batch_region_count(const straw_batch_t *batch) { return batch ? batch->regions : 0; }
size_t straw_batch_record_count(const straw_batch_t *batch) { return batch ? batch->x.size() : 0; }
const size_t *straw_batch_offsets(const straw_batch_t *batch) { return batch && !batch->offsets.empty() ? batch->offsets.data() : nullptr; }
const int64_t *straw_batch_x(const straw_batch_t *batch) { return batch && !batch->x.empty() ? batch->x.data() : nullptr; }
const int64_t *straw_batch_y(const straw_batch_t *batch) { return batch && !batch->y.empty() ? batch->y.data() : nullptr; }
const float *straw_batch_values(const straw_batch_t *batch) { return batch && !batch->values.empty() ? batch->values.data() : nullptr; }
void straw_batch_free(straw_batch_t *batch) { delete batch; }

straw_status_t straw_file_record_count(const straw_file_t *file, int32_t resolution, int inter_only,
                                       uint64_t *out_count, straw_error_t **out_error) {
    if (out_count) *out_count = 0;
    clear_error(out_error);
    if (!file || resolution <= 0 || !out_count)
        return fail(STRAW_STATUS_INVALID_ARGUMENT, "invalid record count arguments", out_error);
    return guarded(out_error, [&] {
        int64_t count = getNumRecordsForFile(file->path, resolution, inter_only != 0);
        if (count < 0) throw StrawException(StrawErrorCode::CorruptFile,
                                            "negative native record count");
        *out_count = static_cast<uint64_t>(count);
    });
}

straw_status_t straw_file_chromosome_record_counts(const straw_file_t *file, int32_t resolution,
                                                   straw_chromosome_counts_t **out_counts,
                                                   straw_error_t **out_error) {
    if (out_counts) *out_counts = nullptr;
    clear_error(out_error);
    if (!file || resolution <= 0 || !out_counts)
        return fail(STRAW_STATUS_INVALID_ARGUMENT, "invalid record count arguments", out_error);
    return guarded(out_error, [&] {
        std::vector<std::pair<std::string, int64_t>> source =
            getRecordCountsByChromosome(file->path, resolution);
        straw_chromosome_counts *result = new straw_chromosome_counts;
        try {
            result->names.reserve(source.size());
            result->values.reserve(source.size());
            for (const auto &entry : source) {
                if (entry.second < 0) throw std::runtime_error("negative native record count");
                result->names.push_back(entry.first);
                result->values.push_back(static_cast<uint64_t>(entry.second));
            }
        } catch (...) {
            delete result;
            throw;
        }
        *out_counts = result;
    });
}

size_t straw_chromosome_counts_size(const straw_chromosome_counts_t *counts) {
    return counts ? counts->names.size() : 0;
}
const char *straw_chromosome_counts_name(const straw_chromosome_counts_t *counts, size_t index) {
    return counts && index < counts->names.size() ? counts->names[index].c_str() : "";
}
const uint64_t *straw_chromosome_counts_values(const straw_chromosome_counts_t *counts) {
    return counts && !counts->values.empty() ? counts->values.data() : nullptr;
}
void straw_chromosome_counts_free(straw_chromosome_counts_t *counts) { delete counts; }

} // extern "C"
