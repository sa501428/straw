use ExtUtils::MakeMaker;
WriteMakefile(
    NAME => 'Hic::Straw',
    VERSION_FROM => 'lib/Hic/Straw.pm',
    PREREQ_PM => {'FFI::Platypus' => '2.00'},
    LICENSE => 'mit',
);
