#include "mex.h"
#include "straw/straw.h"

#include <cstdint>
#include <cstring>
#include <string>

namespace {
class Text {
public:
    explicit Text(const mxArray *value) : data(mxArrayToUTF8String(value)) {
        if (!data) mexErrMsgIdAndTxt("straw:argument", "Expected a UTF-8 string argument");
    }
    ~Text() { mxFree(data); }
    const char *get() const { return data; }
private:
    char *data;
};

[[noreturn]] void raise(straw_status_t status, straw_error_t *error) {
    std::string message = error ? straw_error_message(error) : straw_status_name(status);
    straw_error_free(error);
    mexErrMsgIdAndTxt("straw:native", "%s", message.c_str());
}

straw_file_t *open_file(const char *path) {
    straw_file_t *file = nullptr; straw_error_t *error = nullptr;
    straw_status_t status = straw_file_open(path, &file, &error);
    if (status != STRAW_STATUS_OK) raise(status, error);
    return file;
}

mxArray *int64_column(const int64_t *values, size_t count) {
    mxArray *result = mxCreateNumericMatrix(count, 1, mxINT64_CLASS, mxREAL);
    if (count) std::memcpy(mxGetData(result), values, count * sizeof(int64_t));
    return result;
}
mxArray *uint64_column(const uint64_t *values, size_t count) {
    mxArray *result = mxCreateNumericMatrix(count, 1, mxUINT64_CLASS, mxREAL);
    if (count) std::memcpy(mxGetData(result), values, count * sizeof(uint64_t));
    return result;
}
mxArray *single_column(const float *values, size_t count) {
    mxArray *result = mxCreateNumericMatrix(count, 1, mxSINGLE_CLASS, mxREAL);
    if (count) std::memcpy(mxGetData(result), values, count * sizeof(float));
    return result;
}

void return_records(int nlhs, mxArray **plhs, straw_records_t *records) {
    const size_t count = straw_records_size(records);
    mxArray *x = int64_column(straw_records_x(records), count);
    mxArray *y = int64_column(straw_records_y(records), count);
    mxArray *values = single_column(straw_records_values(records), count);
    straw_records_free(records);
    if (nlhs <= 1) {
        const char *fields[] = {"x", "y", "counts"};
        plhs[0] = mxCreateStructMatrix(1, 1, 3, fields);
        mxSetField(plhs[0], 0, "x", x); mxSetField(plhs[0], 0, "y", y);
        mxSetField(plhs[0], 0, "counts", values);
    } else {
        plhs[0] = x; plhs[1] = y;
        if (nlhs > 2) plhs[2] = values; else mxDestroyArray(values);
    }
}

void query(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs) {
    if (nrhs < 6 || nrhs > 7)
        mexErrMsgIdAndTxt("straw:usage", "Usage: straw(norm,file,chr1,chr2,unit,resolution[,matrix])");
    Text norm(prhs[0]), path(prhs[1]), first(prhs[2]), second(prhs[3]), unit(prhs[4]);
    std::string matrix = "observed";
    if (nrhs == 7) { Text supplied(prhs[6]); matrix = supplied.get(); }
    straw_file_t *file = open_file(path.get());
    straw_query_options_t options = {sizeof(options), matrix.c_str(), norm.get(), first.get(),
        second.get(), unit.get(), static_cast<int32_t>(mxGetScalar(prhs[5]))};
    straw_records_t *records = nullptr; straw_error_t *error = nullptr;
    straw_status_t status = straw_query_records(file, &options, &records, &error);
    straw_file_close(file);
    if (status != STRAW_STATUS_OK) raise(status, error);
    return_records(nlhs, plhs, records);
}

void raw(mxArray **plhs, int nrhs, const mxArray **prhs) {
    if (nrhs != 6)
        mexErrMsgIdAndTxt("straw:usage", "Usage: straw_mex('raw',file,chr1,chr2,unit,resolution)");
    Text path(prhs[1]), first(prhs[2]), second(prhs[3]), unit(prhs[4]);
    straw_file_t *file = open_file(path.get());
    straw_raw_records_t *records = nullptr; straw_error_t *error = nullptr;
    straw_status_t status = straw_query_raw(file, first.get(), second.get(), unit.get(),
        static_cast<int32_t>(mxGetScalar(prhs[5])), &records, &error);
    straw_file_close(file);
    if (status != STRAW_STATUS_OK) raise(status, error);
    const size_t count = straw_raw_records_size(records);
    const char *fields[] = {"binX", "binY", "count", "score", "isScore"};
    plhs[0] = mxCreateStructMatrix(1, 1, 5, fields);
    mxSetField(plhs[0], 0, "binX", uint64_column(straw_raw_records_x(records), count));
    mxSetField(plhs[0], 0, "binY", uint64_column(straw_raw_records_y(records), count));
    mxSetField(plhs[0], 0, "count", uint64_column(straw_raw_records_counts(records), count));
    mxSetField(plhs[0], 0, "score", single_column(straw_raw_records_scores(records), count));
    mxArray *kinds = mxCreateLogicalMatrix(count, 1);
    for (size_t i = 0; i < count; ++i)
        static_cast<mxLogical *>(mxGetData(kinds))[i] = straw_raw_records_kinds(records)[i] != 0;
    mxSetField(plhs[0], 0, "isScore", kinds);
    straw_raw_records_free(records);
}

void metadata(mxArray **plhs, int nrhs, const mxArray **prhs) {
    if (nrhs != 2) mexErrMsgIdAndTxt("straw:usage", "Usage: straw_mex('metadata',file)");
    Text path(prhs[1]); straw_file_t *file = open_file(path.get());
    const char *fields[] = {"version", "genome", "chromosomes", "bpResolutions",
                            "fragResolutions", "normalizations"};
    plhs[0] = mxCreateStructMatrix(1, 1, 6, fields);
    mxSetField(plhs[0], 0, "version", mxCreateDoubleScalar(straw_file_format_version(file)));
    mxSetField(plhs[0], 0, "genome", mxCreateString(straw_file_genome(file)));
    const char *chrom_fields[] = {"index", "name", "length"};
    size_t count = straw_file_chromosome_count(file);
    mxArray *chromosomes = mxCreateStructMatrix(count, 1, 3, chrom_fields);
    for (size_t i = 0; i < count; ++i) {
        mxSetField(chromosomes, i, "index", mxCreateDoubleScalar(straw_file_chromosome_index(file, i)));
        mxSetField(chromosomes, i, "name", mxCreateString(straw_file_chromosome_name(file, i)));
        mxArray *length = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
        *static_cast<int64_t *>(mxGetData(length)) = straw_file_chromosome_length(file, i);
        mxSetField(chromosomes, i, "length", length);
    }
    mxSetField(plhs[0], 0, "chromosomes", chromosomes);
    const char *units[] = {"BP", "FRAG"};
    for (const char *unit_name : units) {
        count = straw_file_resolution_count(file, unit_name);
        mxArray *values = mxCreateNumericMatrix(count, 1, mxINT32_CLASS, mxREAL);
        for (size_t i = 0; i < count; ++i)
            static_cast<int32_t *>(mxGetData(values))[i] = straw_file_resolution(file, unit_name, i);
        mxSetField(plhs[0], 0, unit_name[0] == 'B' ? "bpResolutions" : "fragResolutions", values);
    }
    count = straw_file_normalization_count(file);
    mxArray *norms = mxCreateCellMatrix(count, 1);
    for (size_t i = 0; i < count; ++i)
        mxSetCell(norms, i, mxCreateString(straw_file_normalization(file, i)));
    mxSetField(plhs[0], 0, "normalizations", norms);
    straw_file_close(file);
}
} // namespace

void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs) {
    if (nrhs == 0) mexErrMsgIdAndTxt("straw:usage", "No arguments supplied");
    if (mxIsChar(prhs[0])) {
        Text command(prhs[0]);
        if (std::string(command.get()) == "metadata") return metadata(plhs, nrhs, prhs);
        if (std::string(command.get()) == "raw") return raw(plhs, nrhs, prhs);
    }
    query(nlhs, plhs, nrhs, prhs);
}
