package Hic::Straw;
use strict;
use warnings;
use FFI::Platypus 2.00;
use FFI::Platypus::Buffer qw(buffer_to_scalar);
use File::Basename qw(dirname);
use File::Spec;

our $VERSION = '1.0.0';
my $ffi = FFI::Platypus->new(api => 2);
my $rid = $^O eq 'MSWin32' ? 'win-x64' : $^O eq 'darwin' ? 'osx-arm64' : 'linux-x64';
my $filename = $^O eq 'MSWin32' ? 'straw.dll' : $^O eq 'darwin' ? 'libstraw.dylib' : 'libstraw.so';
my $bundled = File::Spec->catfile(dirname(__FILE__), 'Straw', 'native', $rid, $filename);
$ffi->lib($ENV{LIBSTRAW_PATH} // (-f $bundled ? $bundled : 'straw'));

$ffi->attach(straw_file_open => ['string', 'opaque*', 'opaque*'] => 'int');
$ffi->attach(straw_file_close => ['opaque'] => 'void');
$ffi->attach(straw_file_format_version => ['opaque'] => 'sint32');
$ffi->attach(straw_file_genome => ['opaque'] => 'string');
$ffi->attach(straw_file_chromosome_count => ['opaque'] => 'size_t');
$ffi->attach(straw_file_chromosome_name => ['opaque', 'size_t'] => 'string');
$ffi->attach(straw_file_chromosome_index => ['opaque', 'size_t'] => 'sint32');
$ffi->attach(straw_file_chromosome_length => ['opaque', 'size_t'] => 'sint64');
$ffi->attach(straw_file_resolution_count => ['opaque', 'string'] => 'size_t');
$ffi->attach(straw_file_resolution => ['opaque', 'string', 'size_t'] => 'sint32');
$ffi->attach(straw_query_records_simple =>
    ['opaque', ('string') x 5, 'sint32', 'opaque*', 'opaque*'] => 'int');
$ffi->attach(straw_records_size => ['opaque'] => 'size_t');
$ffi->attach(straw_records_x => ['opaque'] => 'opaque');
$ffi->attach(straw_records_y => ['opaque'] => 'opaque');
$ffi->attach(straw_records_values => ['opaque'] => 'opaque');
$ffi->attach(straw_records_free => ['opaque'] => 'void');
$ffi->attach(straw_error_message => ['opaque'] => 'string');
$ffi->attach(straw_error_free => ['opaque'] => 'void');

sub _check {
    my ($status, $error) = @_;
    return if $status == 0;
    my $message = $error ? straw_error_message($error) : "libstraw error $status";
    straw_error_free($error) if $error;
    die "$message\n";
}

sub open {
    my ($class, $path) = @_;
    my ($handle, $error);
    _check(straw_file_open($path, \$handle, \$error), $error);
    return bless { handle => $handle }, $class;
}

sub close {
    my ($self) = @_;
    if ($self->{handle}) {
        straw_file_close($self->{handle});
        $self->{handle} = undef;
    }
}

sub DESTROY { shift->close }
sub version { straw_file_format_version($_[0]->{handle}) }
sub genome { straw_file_genome($_[0]->{handle}) }

sub chromosomes {
    my ($self) = @_;
    my @values;
    for my $i (0 .. straw_file_chromosome_count($self->{handle}) - 1) {
        push @values, {
            index => straw_file_chromosome_index($self->{handle}, $i),
            name => straw_file_chromosome_name($self->{handle}, $i),
            length => straw_file_chromosome_length($self->{handle}, $i),
        };
    }
    return \@values;
}

sub resolutions {
    my ($self, $unit) = @_; $unit //= 'BP';
    my $count = straw_file_resolution_count($self->{handle}, $unit);
    return [] unless $count;
    return [map { straw_file_resolution($self->{handle}, $unit, $_) } 0 .. $count - 1];
}

sub _array {
    my ($pointer, $bytes, $template) = @_;
    return [] unless $bytes;
    my $data = buffer_to_scalar($pointer, $bytes);
    return [unpack($template, $data)];
}

sub records {
    my ($self, $first, $second, %options) = @_;
    my ($owner, $error);
    _check(straw_query_records_simple($self->{handle}, $options{matrix_type} // 'observed',
        $options{normalization} // 'NONE', $first, $second, $options{unit} // 'BP',
        $options{resolution}, \$owner, \$error), $error);
    my $count = straw_records_size($owner);
    my $result = {
        x => _array(straw_records_x($owner), $count * 8, 'q*'),
        y => _array(straw_records_y($owner), $count * 8, 'q*'),
        values => _array(straw_records_values($owner), $count * 4, 'f*'),
    };
    straw_records_free($owner);
    return $result;
}

1;
