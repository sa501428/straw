using System.Runtime.InteropServices;
using Microsoft.Win32.SafeHandles;

namespace HicStraw;

internal static class Native
{
    private const string Library = "straw";

    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern int straw_file_open(
        [MarshalAs(UnmanagedType.LPUTF8Str)] string path, out IntPtr file, out IntPtr error);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void straw_file_close(IntPtr file);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern int straw_file_format_version(SafeFileHandle file);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr straw_file_genome(SafeFileHandle file);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern nuint straw_file_chromosome_count(SafeFileHandle file);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr straw_file_chromosome_name(SafeFileHandle file, nuint index);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern int straw_file_chromosome_index(SafeFileHandle file, nuint index);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern long straw_file_chromosome_length(SafeFileHandle file, nuint index);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern nuint straw_file_resolution_count(
        SafeFileHandle file, [MarshalAs(UnmanagedType.LPUTF8Str)] string unit);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern int straw_file_resolution(
        SafeFileHandle file, [MarshalAs(UnmanagedType.LPUTF8Str)] string unit, nuint index);

    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern int straw_query_records_simple(
        SafeFileHandle file,
        [MarshalAs(UnmanagedType.LPUTF8Str)] string matrixType,
        [MarshalAs(UnmanagedType.LPUTF8Str)] string normalization,
        [MarshalAs(UnmanagedType.LPUTF8Str)] string first,
        [MarshalAs(UnmanagedType.LPUTF8Str)] string second,
        [MarshalAs(UnmanagedType.LPUTF8Str)] string unit,
        int resolution, out IntPtr records, out IntPtr error);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern nuint straw_records_size(SafeRecordsHandle records);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr straw_records_x(SafeRecordsHandle records);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr straw_records_y(SafeRecordsHandle records);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr straw_records_values(SafeRecordsHandle records);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void straw_records_free(IntPtr records);

    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern int straw_query_raw(SafeFileHandle file,
        [MarshalAs(UnmanagedType.LPUTF8Str)] string first,
        [MarshalAs(UnmanagedType.LPUTF8Str)] string second,
        [MarshalAs(UnmanagedType.LPUTF8Str)] string unit,
        int resolution, out IntPtr records, out IntPtr error);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern nuint straw_raw_records_size(SafeRawRecordsHandle records);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr straw_raw_records_x(SafeRawRecordsHandle records);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr straw_raw_records_y(SafeRawRecordsHandle records);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr straw_raw_records_counts(SafeRawRecordsHandle records);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr straw_raw_records_scores(SafeRawRecordsHandle records);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr straw_raw_records_kinds(SafeRawRecordsHandle records);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void straw_raw_records_free(IntPtr records);

    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr straw_error_message(IntPtr error);
    [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void straw_error_free(IntPtr error);

    internal static string Utf8(IntPtr value) =>
        value == IntPtr.Zero ? string.Empty : Marshal.PtrToStringUTF8(value) ?? string.Empty;

    internal static void ThrowIfError(int status, IntPtr error)
    {
        if (status == 0) return;
        string message = error == IntPtr.Zero ? $"libstraw error {status}" : Utf8(straw_error_message(error));
        if (error != IntPtr.Zero) straw_error_free(error);
        throw new StrawException(status, message);
    }
}

internal sealed class SafeFileHandle : SafeHandleZeroOrMinusOneIsInvalid
{
    private SafeFileHandle() : base(true) { }
    internal SafeFileHandle(IntPtr handle) : base(true) => SetHandle(handle);
    protected override bool ReleaseHandle() { Native.straw_file_close(handle); return true; }
}

internal sealed class SafeRecordsHandle : SafeHandleZeroOrMinusOneIsInvalid
{
    private SafeRecordsHandle() : base(true) { }
    internal SafeRecordsHandle(IntPtr handle) : base(true) => SetHandle(handle);
    protected override bool ReleaseHandle() { Native.straw_records_free(handle); return true; }
}

internal sealed class SafeRawRecordsHandle : SafeHandleZeroOrMinusOneIsInvalid
{
    private SafeRawRecordsHandle() : base(true) { }
    internal SafeRawRecordsHandle(IntPtr handle) : base(true) => SetHandle(handle);
    protected override bool ReleaseHandle() { Native.straw_raw_records_free(handle); return true; }
}
