#!/usr/bin/perl -w
use strict;


# fireoff from raw data project folder



my $genome = shift @ARGV; #can also accept settings for run_bwa_cluster script after genome

die "provide ref genome"  if !$genome;
die "provide ref genome" unless -e "$genome.bwt";

my $pwd = `pwd`;
my @folder=split(/\//, $pwd) ;
my $folder;
my $targetpath;


if ($pwd =~ /hiseq/) {
    $folder=$folder[-3];
    die "you are not in the correct project folder!!\n" unless $folder[-2] eq 'Unaligned';

    $targetpath = "/hpc/cog_bioinf/data/mapping/hiseq";

    print  "Working on $folder\n";
}elsif ($pwd =~ /miseq/) {
    $targetpath = "/hpc/cog_bioinf/data/mapping/miseq";
    $folder = $folder[-1];
    print $folder,"\n";
}



system "mkdir $targetpath/$folder";
my @files = `find \$PWD -iname \"*fastq.gz\"`;
chdir "$targetpath/$folder";
system "ln -s $_" foreach @files;




my %folders;
# prepare mappign structure mysqeq run files

foreach my $file (<*fastq.gz>) {
    my $filecore = $file;
    $filecore =~ s/_L00\d_R?[1|2]_001.fastq.gz//;
    print "processing $file\t\t $filecore\n";
    system "mkdir -p $filecore/reads";
    system "mv $file $filecore/reads";
    $folders{$filecore}++;
}

foreach my $filecore (keys %folders) {
    #execute mapping
    system "cp /hpc/cog_bioinf/data/ies/scripts/illumina_proc/run_bwa_cluster_miseq_v3.2.pl $filecore";
    chdir $filecore or die "cannot change dir $filecore\n";;
    system "perl run_bwa_cluster_miseq_v3.2.pl -genome $genome @ARGV";
    chdir '../';
}