module Straw

using Libdl

export HicFile, Chromosome, ContactRecords, RawRecords, PreparedQuery,
       chromosomes, resolutions, normalizations, records, dense,
       normalization_vector, expected_vector, raw_records, prepare, query, query_regions

const libstraw = get(ENV, "LIBSTRAW_PATH", let
    found = Libdl.find_library(["straw", "libstraw"])
    isempty(found) ? "libstraw" : found
end)

const Handle = Ptr{Cvoid}

struct CQueryOptions
    struct_size::Csize_t
    matrix_type::Cstring
    normalization::Cstring
    first_location::Cstring
    second_location::Cstring
    unit::Cstring
    resolution::Int32
end

struct CRegion
    x_start::Int64
    x_end::Int64
    y_start::Int64
    y_end::Int64
end

struct Chromosome
    index::Int32
    name::String
    length::Int64
end

struct ContactRecords
    x::Vector{Int64}
    y::Vector{Int64}
    values::Vector{Float32}
end

struct RawRecords
    x::Vector{UInt64}
    y::Vector{UInt64}
    counts::Vector{UInt64}
    scores::Vector{Float32}
    is_score::BitVector
end

mutable struct HicFile
    handle::Handle
end

mutable struct PreparedQuery
    handle::Handle
end

function native_error(status::Integer, error::Handle)
    status == 0 && return
    message = error == C_NULL ? "libstraw error $(status)" : unsafe_string(
        ccall((:straw_error_message, libstraw), Cstring, (Handle,), error))
    error != C_NULL && ccall((:straw_error_free, libstraw), Cvoid, (Handle,), error)
    throw(ErrorException(message))
end

function HicFile(path::AbstractString)
    output = Ref{Handle}(C_NULL)
    error = Ref{Handle}(C_NULL)
    status = ccall((:straw_file_open, libstraw), Cint,
                   (Cstring, Ref{Handle}, Ref{Handle}), path, output, error)
    native_error(status, error[])
    file = HicFile(output[])
    finalizer(close, file)
    file
end

function Base.close(file::HicFile)
    if file.handle != C_NULL
        ccall((:straw_file_close, libstraw), Cvoid, (Handle,), file.handle)
        file.handle = C_NULL
    end
    nothing
end

function Base.close(query::PreparedQuery)
    if query.handle != C_NULL
        ccall((:straw_query_close, libstraw), Cvoid, (Handle,), query.handle)
        query.handle = C_NULL
    end
    nothing
end

function chromosomes(file::HicFile)
    n = ccall((:straw_file_chromosome_count, libstraw), Csize_t, (Handle,), file.handle)
    [Chromosome(
        ccall((:straw_file_chromosome_index, libstraw), Int32, (Handle, Csize_t), file.handle, i),
        unsafe_string(ccall((:straw_file_chromosome_name, libstraw), Cstring,
                            (Handle, Csize_t), file.handle, i)),
        ccall((:straw_file_chromosome_length, libstraw), Int64,
              (Handle, Csize_t), file.handle, i)) for i in 0:n-1]
end

function resolutions(file::HicFile, unit::AbstractString="BP")
    n = ccall((:straw_file_resolution_count, libstraw), Csize_t,
              (Handle, Cstring), file.handle, unit)
    [ccall((:straw_file_resolution, libstraw), Int32,
           (Handle, Cstring, Csize_t), file.handle, unit, i) for i in 0:n-1]
end

function normalizations(file::HicFile)
    n = ccall((:straw_file_normalization_count, libstraw), Csize_t, (Handle,), file.handle)
    [unsafe_string(ccall((:straw_file_normalization, libstraw), Cstring,
                         (Handle, Csize_t), file.handle, i)) for i in 0:n-1]
end

copy_native(::Type{T}, pointer::Ptr{T}, n::Integer) where {T} =
    n == 0 ? T[] : copy(unsafe_wrap(Vector{T}, pointer, n; own=false))

function records(file::HicFile, first::AbstractString, second::AbstractString;
                 matrix_type::AbstractString="observed", normalization::AbstractString="NONE",
                 unit::AbstractString="BP", resolution::Integer)
    output = Ref{Handle}(C_NULL)
    error = Ref{Handle}(C_NULL)
    GC.@preserve matrix_type normalization first second unit begin
        options = CQueryOptions(sizeof(CQueryOptions), Base.unsafe_convert(Cstring, matrix_type),
            Base.unsafe_convert(Cstring, normalization), Base.unsafe_convert(Cstring, first),
            Base.unsafe_convert(Cstring, second), Base.unsafe_convert(Cstring, unit), Int32(resolution))
        status = ccall((:straw_query_records, libstraw), Cint,
                       (Handle, Ref{CQueryOptions}, Ref{Handle}, Ref{Handle}),
                       file.handle, options, output, error)
        native_error(status, error[])
    end
    owner = output[]
    try
        n = ccall((:straw_records_size, libstraw), Csize_t, (Handle,), owner)
        ContactRecords(
            copy_native(Int64, ccall((:straw_records_x, libstraw), Ptr{Int64}, (Handle,), owner), n),
            copy_native(Int64, ccall((:straw_records_y, libstraw), Ptr{Int64}, (Handle,), owner), n),
            copy_native(Float32, ccall((:straw_records_values, libstraw), Ptr{Float32}, (Handle,), owner), n))
    finally
        ccall((:straw_records_free, libstraw), Cvoid, (Handle,), owner)
    end
end

function dense(file::HicFile, first::AbstractString, second::AbstractString;
               matrix_type::AbstractString="observed", normalization::AbstractString="NONE",
               unit::AbstractString="BP", resolution::Integer)
    output = Ref{Handle}(C_NULL); error = Ref{Handle}(C_NULL)
    GC.@preserve matrix_type normalization first second unit begin
        options = CQueryOptions(sizeof(CQueryOptions), Base.unsafe_convert(Cstring, matrix_type),
            Base.unsafe_convert(Cstring, normalization), Base.unsafe_convert(Cstring, first),
            Base.unsafe_convert(Cstring, second), Base.unsafe_convert(Cstring, unit), Int32(resolution))
        status = ccall((:straw_query_dense, libstraw), Cint,
                       (Handle, Ref{CQueryOptions}, Ref{Handle}, Ref{Handle}),
                       file.handle, options, output, error)
        native_error(status, error[])
    end
    owner = output[]
    try
        rows = ccall((:straw_dense_rows, libstraw), Csize_t, (Handle,), owner)
        columns = ccall((:straw_dense_columns, libstraw), Csize_t, (Handle,), owner)
        values = copy_native(Float32, ccall((:straw_dense_values_row_major, libstraw),
                                             Ptr{Float32}, (Handle,), owner), rows * columns)
        permutedims(reshape(values, columns, rows), (2, 1))
    finally
        ccall((:straw_dense_matrix_free, libstraw), Cvoid, (Handle,), owner)
    end
end

function native_vector(symbol::Symbol, file::HicFile, chromosome::AbstractString,
                       unit::AbstractString, resolution::Integer, normalization::AbstractString)
    output = Ref{Handle}(C_NULL); error = Ref{Handle}(C_NULL)
    status = if symbol === :straw_file_normalization_vector
        ccall((:straw_file_normalization_vector, libstraw), Cint,
              (Handle, Cstring, Cstring, Int32, Cstring, Ref{Handle}, Ref{Handle}),
              file.handle, chromosome, unit, Int32(resolution), normalization, output, error)
    else
        ccall((:straw_file_expected_vector, libstraw), Cint,
              (Handle, Cstring, Cstring, Int32, Cstring, Ref{Handle}, Ref{Handle}),
              file.handle, chromosome, unit, Int32(resolution), normalization, output, error)
    end
    native_error(status, error[])
    owner = output[]
    try
        n = ccall((:straw_vector_size, libstraw), Csize_t, (Handle,), owner)
        copy_native(Float64, ccall((:straw_vector_values, libstraw), Ptr{Float64}, (Handle,), owner), n)
    finally
        ccall((:straw_vector_free, libstraw), Cvoid, (Handle,), owner)
    end
end

normalization_vector(file::HicFile, chromosome::AbstractString; unit="BP", resolution, normalization="NONE") =
    native_vector(:straw_file_normalization_vector, file, chromosome, unit, resolution, normalization)
expected_vector(file::HicFile, chromosome::AbstractString; unit="BP", resolution, normalization="NONE") =
    native_vector(:straw_file_expected_vector, file, chromosome, unit, resolution, normalization)

function raw_records(file::HicFile, first::AbstractString, second::AbstractString;
                     unit::AbstractString="BP", resolution::Integer)
    output = Ref{Handle}(C_NULL); error = Ref{Handle}(C_NULL)
    status = ccall((:straw_query_raw, libstraw), Cint,
                   (Handle, Cstring, Cstring, Cstring, Int32, Ref{Handle}, Ref{Handle}),
                   file.handle, first, second, unit, Int32(resolution), output, error)
    native_error(status, error[])
    owner = output[]
    try
        n = ccall((:straw_raw_records_size, libstraw), Csize_t, (Handle,), owner)
        kinds = copy_native(UInt8, ccall((:straw_raw_records_kinds, libstraw), Ptr{UInt8}, (Handle,), owner), n)
        RawRecords(
            copy_native(UInt64, ccall((:straw_raw_records_x, libstraw), Ptr{UInt64}, (Handle,), owner), n),
            copy_native(UInt64, ccall((:straw_raw_records_y, libstraw), Ptr{UInt64}, (Handle,), owner), n),
            copy_native(UInt64, ccall((:straw_raw_records_counts, libstraw), Ptr{UInt64}, (Handle,), owner), n),
            copy_native(Float32, ccall((:straw_raw_records_scores, libstraw), Ptr{Float32}, (Handle,), owner), n),
            BitVector(kinds .!= 0))
    finally
        ccall((:straw_raw_records_free, libstraw), Cvoid, (Handle,), owner)
    end
end

function prepare(file::HicFile, first_chromosome::AbstractString, second_chromosome::AbstractString;
                 matrix_type::AbstractString="observed", normalization::AbstractString="NONE",
                 unit::AbstractString="BP", resolution::Integer)
    output = Ref{Handle}(C_NULL); error = Ref{Handle}(C_NULL)
    status = ccall((:straw_query_prepare, libstraw), Cint,
                   (Handle, Cstring, Cstring, Cstring, Cstring, Cstring, Int32,
                    Ref{Handle}, Ref{Handle}), file.handle, matrix_type, normalization,
                   first_chromosome, second_chromosome, unit, Int32(resolution), output, error)
    native_error(status, error[])
    prepared = PreparedQuery(output[])
    finalizer(close, prepared)
    prepared
end

function query(prepared::PreparedQuery, region::NTuple{4,<:Integer})
    native = CRegion(map(Int64, region)...)
    output = Ref{Handle}(C_NULL); error = Ref{Handle}(C_NULL)
    status = ccall((:straw_query_window, libstraw), Cint,
                   (Handle, Ref{CRegion}, Ref{Handle}, Ref{Handle}),
                   prepared.handle, native, output, error)
    native_error(status, error[])
    owner = output[]
    try
        n = ccall((:straw_records_size, libstraw), Csize_t, (Handle,), owner)
        ContactRecords(copy_native(Int64, ccall((:straw_records_x, libstraw), Ptr{Int64}, (Handle,), owner), n),
                       copy_native(Int64, ccall((:straw_records_y, libstraw), Ptr{Int64}, (Handle,), owner), n),
                       copy_native(Float32, ccall((:straw_records_values, libstraw), Ptr{Float32}, (Handle,), owner), n))
    finally
        ccall((:straw_records_free, libstraw), Cvoid, (Handle,), owner)
    end
end

function query_regions(prepared::PreparedQuery, regions)
    native = CRegion[CRegion(map(Int64, region)...) for region in regions]
    output = Ref{Handle}(C_NULL); error = Ref{Handle}(C_NULL)
    status = ccall((:straw_query_regions, libstraw), Cint,
                   (Handle, Ptr{CRegion}, Csize_t, Ref{Handle}, Ref{Handle}),
                   prepared.handle, native, length(native), output, error)
    native_error(status, error[])
    owner = output[]
    try
        count = ccall((:straw_batch_record_count, libstraw), Csize_t, (Handle,), owner)
        offsets = copy_native(Csize_t, ccall((:straw_batch_offsets, libstraw), Ptr{Csize_t}, (Handle,), owner), length(native) + 1)
        records = ContactRecords(
            copy_native(Int64, ccall((:straw_batch_x, libstraw), Ptr{Int64}, (Handle,), owner), count),
            copy_native(Int64, ccall((:straw_batch_y, libstraw), Ptr{Int64}, (Handle,), owner), count),
            copy_native(Float32, ccall((:straw_batch_values, libstraw), Ptr{Float32}, (Handle,), owner), count))
        (offsets=Int.(offsets .+ 1), records=records)
    finally
        ccall((:straw_batch_free, libstraw), Cvoid, (Handle,), owner)
    end
end

end
