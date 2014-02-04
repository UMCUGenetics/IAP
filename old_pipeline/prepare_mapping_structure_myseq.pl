#!/usr/bin/perl -w
use strict;

my $genome = $ARGV[0];

die "provide ref genome"  if !$genome;
die "provide ref genome" unless -e "$genome.bwt";

my %folders;
# prepare mappign structure mysqeq run files

foreach my $file (<*fastq.gz>) {
    my $filecore = $file;
    $filecore =~ s/_L00\d_R?[1|2]_001.fastq.gz//;
    print "processing $file\t\t $filecore\n";
    system "mkdir -p $filecore/reads";
    system "mv $file $filecore/reads";
    #system "gunzip $filecore/reads/*.gz";
    $folders{$filecore}++;
}

foreach my $filecore (keys %folders) {
    #execute mapping
    system "cp /hpc/cog_bioinf/data/ies/scripts/illumina_proc/run_bwa_cluster_miseq_v3.1.pl $filecore";
    chdir $filecore or die "cannot change dir $filecore\n";;
    system "perl run_bwa_cluster_miseq_v3.1.pl -genome $genome";
    chdir '../';
}