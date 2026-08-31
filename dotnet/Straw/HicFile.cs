using System.Runtime.InteropServices;

namespace HicStraw;

public sealed class StrawException(int status, string message) : Exception(message)
{
    public int Status { get; } = status;
}

public sealed record Chromosome(int Index, string Name, long Length);
public sealed record ContactRecords(long[] X, long[] Y, float[] Values);
public sealed record RawRecords(ulong[] X, ulong[] Y, ulong[] Counts, float[] Scores, bool[] IsScore);

public sealed class HicFile : IDisposable
{
    private readonly SafeFileHandle handle;

    public HicFile(string path)
    {
        int status = Native.straw_file_open(path, out IntPtr native, out IntPtr error);
        Native.ThrowIfError(status, error);
        handle = new SafeFileHandle(native);
    }

    public int FormatVersion => Native.straw_file_format_version(handle);
    public string Genome => Native.Utf8(Native.straw_file_genome(handle));

    public IReadOnlyList<Chromosome> Chromosomes
    {
        get
        {
            int count = checked((int)Native.straw_file_chromosome_count(handle));
            var result = new Chromosome[count];
            for (int i = 0; i < count; ++i)
                result[i] = new Chromosome(
                    Native.straw_file_chromosome_index(handle, (nuint)i),
                    Native.Utf8(Native.straw_file_chromosome_name(handle, (nuint)i)),
                    Native.straw_file_chromosome_length(handle, (nuint)i));
            return result;
        }
    }

    public int[] Resolutions(string unit = "BP")
    {
        int count = checked((int)Native.straw_file_resolution_count(handle, unit));
        var values = new int[count];
        for (int i = 0; i < count; ++i)
            values[i] = Native.straw_file_resolution(handle, unit, (nuint)i);
        return values;
    }

    public ContactRecords Records(string first, string second, int resolution,
        string normalization = "NONE", string matrixType = "observed", string unit = "BP")
    {
        int status = Native.straw_query_records_simple(handle, matrixType, normalization, first,
            second, unit, resolution, out IntPtr native, out IntPtr error);
        Native.ThrowIfError(status, error);
        using var owner = new SafeRecordsHandle(native);
        int count = checked((int)Native.straw_records_size(owner));
        var x = new long[count]; var y = new long[count]; var values = new float[count];
        if (count != 0)
        {
            Marshal.Copy(Native.straw_records_x(owner), x, 0, count);
            Marshal.Copy(Native.straw_records_y(owner), y, 0, count);
            Marshal.Copy(Native.straw_records_values(owner), values, 0, count);
        }
        return new ContactRecords(x, y, values);
    }

    public unsafe RawRecords Raw(string first, string second, int resolution, string unit = "BP")
    {
        int status = Native.straw_query_raw(handle, first, second, unit, resolution,
            out IntPtr native, out IntPtr error);
        Native.ThrowIfError(status, error);
        using var owner = new SafeRawRecordsHandle(native);
        int count = checked((int)Native.straw_raw_records_size(owner));
        var x = new ulong[count]; var y = new ulong[count]; var counts = new ulong[count];
        var scores = new float[count]; var kinds = new bool[count];
        if (count != 0)
        {
            new ReadOnlySpan<ulong>((void*)Native.straw_raw_records_x(owner), count).CopyTo(x);
            new ReadOnlySpan<ulong>((void*)Native.straw_raw_records_y(owner), count).CopyTo(y);
            new ReadOnlySpan<ulong>((void*)Native.straw_raw_records_counts(owner), count).CopyTo(counts);
            new ReadOnlySpan<float>((void*)Native.straw_raw_records_scores(owner), count).CopyTo(scores);
            var nativeKinds = new ReadOnlySpan<byte>((void*)Native.straw_raw_records_kinds(owner), count);
            for (int i = 0; i < count; ++i) kinds[i] = nativeKinds[i] != 0;
        }
        return new RawRecords(x, y, counts, scores, kinds);
    }

    public void Dispose() => handle.Dispose();
}
