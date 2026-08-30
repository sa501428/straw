require "ffi"
require "rbconfig"

module HicStraw
  module Native
    extend FFI::Library
    platform = if RbConfig::CONFIG["host_os"] =~ /darwin/
                 RbConfig::CONFIG["host_cpu"] =~ /arm|aarch64/ ? "osx-arm64" : "osx-x64"
               elsif RbConfig::CONFIG["host_os"] =~ /mswin|mingw/
                 "win-x64"
               else
                 "linux-x64"
               end
    filename = platform.start_with?("win") ? "straw.dll" :
      (platform.start_with?("osx") ? "libstraw.dylib" : "libstraw.so")
    bundled = File.expand_path(File.join(__dir__, "hic_straw", "native", platform, filename))
    ffi_lib [ENV["LIBSTRAW_PATH"], bundled, "straw"].compact

    attach_function :straw_file_open, [:string, :pointer, :pointer], :int
    attach_function :straw_file_close, [:pointer], :void
    attach_function :straw_file_format_version, [:pointer], :int32
    attach_function :straw_file_genome, [:pointer], :string
    attach_function :straw_file_chromosome_count, [:pointer], :size_t
    attach_function :straw_file_chromosome_name, [:pointer, :size_t], :string
    attach_function :straw_file_chromosome_index, [:pointer, :size_t], :int32
    attach_function :straw_file_chromosome_length, [:pointer, :size_t], :int64
    attach_function :straw_file_resolution_count, [:pointer, :string], :size_t
    attach_function :straw_file_resolution, [:pointer, :string, :size_t], :int32

    attach_function :straw_query_records_simple,
                    [:pointer, :string, :string, :string, :string, :string, :int32,
                     :pointer, :pointer], :int
    attach_function :straw_records_size, [:pointer], :size_t
    attach_function :straw_records_x, [:pointer], :pointer
    attach_function :straw_records_y, [:pointer], :pointer
    attach_function :straw_records_values, [:pointer], :pointer
    attach_function :straw_records_free, [:pointer], :void

    attach_function :straw_query_raw,
                    [:pointer, :string, :string, :string, :int32, :pointer, :pointer], :int
    attach_function :straw_raw_records_size, [:pointer], :size_t
    attach_function :straw_raw_records_x, [:pointer], :pointer
    attach_function :straw_raw_records_y, [:pointer], :pointer
    attach_function :straw_raw_records_counts, [:pointer], :pointer
    attach_function :straw_raw_records_scores, [:pointer], :pointer
    attach_function :straw_raw_records_kinds, [:pointer], :pointer
    attach_function :straw_raw_records_free, [:pointer], :void

    attach_function :straw_error_message, [:pointer], :string
    attach_function :straw_error_free, [:pointer], :void
  end

  class Error < StandardError
    attr_reader :status
    def initialize(status, message)
      @status = status
      super(message)
    end
  end

  class FileHandle < FFI::AutoPointer
    def self.release(pointer)
      Native.straw_file_close(pointer)
    end
  end

  def self.check(status, error_output)
    return if status.zero?
    error = error_output.read_pointer
    message = error.null? ? "libstraw error #{status}" : Native.straw_error_message(error)
    Native.straw_error_free(error) unless error.null?
    raise Error.new(status, message)
  end

  Chromosome = Struct.new(:index, :name, :length, keyword_init: true)
  ContactRecords = Struct.new(:x, :y, :values, keyword_init: true)
  RawRecords = Struct.new(:x, :y, :counts, :scores, :is_score, keyword_init: true)

  class File
    def initialize(path)
      output = FFI::MemoryPointer.new(:pointer)
      error = FFI::MemoryPointer.new(:pointer)
      HicStraw.check(Native.straw_file_open(path, output, error), error)
      @handle = FileHandle.new(output.read_pointer)
      @closed = false
    end

    def close
      return if @closed
      @handle.autorelease = false
      Native.straw_file_close(@handle)
      @closed = true
    end

    def version
      Native.straw_file_format_version(@handle)
    end

    def genome
      Native.straw_file_genome(@handle)
    end

    def chromosomes
      Native.straw_file_chromosome_count(@handle).times.map do |index|
        Chromosome.new(index: Native.straw_file_chromosome_index(@handle, index),
                       name: Native.straw_file_chromosome_name(@handle, index),
                       length: Native.straw_file_chromosome_length(@handle, index))
      end
    end

    def resolutions(unit = "BP")
      Native.straw_file_resolution_count(@handle, unit).times.map do |index|
        Native.straw_file_resolution(@handle, unit, index)
      end
    end

    def records(first, second, resolution:, normalization: "NONE",
                matrix_type: "observed", unit: "BP")
      output = FFI::MemoryPointer.new(:pointer)
      error = FFI::MemoryPointer.new(:pointer)
      HicStraw.check(Native.straw_query_records_simple(@handle, matrix_type, normalization,
        first, second, unit, resolution, output, error), error)
      owner = output.read_pointer
      begin
        count = Native.straw_records_size(owner)
        ContactRecords.new(x: read_int64(Native.straw_records_x(owner), count),
                           y: read_int64(Native.straw_records_y(owner), count),
                           values: count.zero? ? [] : Native.straw_records_values(owner).read_array_of_float(count))
      ensure
        Native.straw_records_free(owner)
      end
    end

    def raw(first, second, resolution:, unit: "BP")
      output = FFI::MemoryPointer.new(:pointer)
      error = FFI::MemoryPointer.new(:pointer)
      HicStraw.check(Native.straw_query_raw(@handle, first, second, unit, resolution,
                                            output, error), error)
      owner = output.read_pointer
      begin
        count = Native.straw_raw_records_size(owner)
        kinds = count.zero? ? [] : Native.straw_raw_records_kinds(owner).read_array_of_uint8(count)
        RawRecords.new(x: read_uint64(Native.straw_raw_records_x(owner), count),
                       y: read_uint64(Native.straw_raw_records_y(owner), count),
                       counts: read_uint64(Native.straw_raw_records_counts(owner), count),
                       scores: count.zero? ? [] : Native.straw_raw_records_scores(owner).read_array_of_float(count),
                       is_score: kinds.map { |kind| !kind.zero? })
      ensure
        Native.straw_raw_records_free(owner)
      end
    end

    private

    def read_int64(pointer, count)
      count.zero? ? [] : pointer.read_array_of_type(:int64, :read_int64, count)
    end

    def read_uint64(pointer, count)
      count.zero? ? [] : pointer.read_array_of_type(:uint64, :read_uint64, count)
    end
  end
end
