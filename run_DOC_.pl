#!/usr/bin/perl -w
use strict;

my $request = '/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/sorted_Homo_sapiens.GRCh37.74_nopseudo_noRNA_CDS_GATK.list';



die "cannot open the request file!:$!\n" unless -e $request;

my @BAM = `find -maxdepth 2 -iname "*realigned*.bam"`;
print scalar @BAM, " bams found\n";

chomp ($_) foreach @BAM;

my $command = "java -Xmx12g -jar /hpc/cog_bioinf/common_scripts/GenomeAnalysisTK-2.8-1/GenomeAnalysisTK.jar -T DepthOfCoverage ".
    "-o DepthOfCoverage -R /hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta ".
    "-geneList /hpc/cog_bioinf/ENRICH/GATK/sorted_hg19_ensembl_rod ".
    "-L $request ".
    "--minMappingQuality 1 -ct 1 -ct 10 -ct 15 -ct 20 -ct 40 -ct 60 -I ".
    join(" -I ", @BAM);
    
    
#system $command;

print $command,"\n\n";
