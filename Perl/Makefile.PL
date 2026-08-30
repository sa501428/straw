use ExtUtils::MakeMaker;
my %pm = ('lib/Hic/Straw.pm' => '$(INST_LIB)/Hic/Straw.pm');
for my $native (glob 'lib/Hic/Straw/native/*/*') {
    (my $destination = $native) =~ s{^lib/}{};
    $pm{$native} = '$(INST_LIB)/' . $destination;
}
WriteMakefile(
    NAME => 'Hic::Straw',
    VERSION_FROM => 'lib/Hic/Straw.pm',
    PREREQ_PM => {'FFI::Platypus' => '2.00'},
    LICENSE => 'mit',
    PM => \%pm,
    META_MERGE => { resources => { repository => 'https://github.com/aidenlab/straw' } },
);
