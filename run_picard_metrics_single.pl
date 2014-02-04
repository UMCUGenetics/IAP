#!/usr/bin/perl -w
use strict;
use POSIX qw(tmpnam);


#collect several picard metrics files
#	single bams files
#

my $picard = "java -Xmx16G -jar /hpc/cog_bioinf/common_scripts/picard-tools-1.98";
my $targets = '';
my $baits ='';

my $ref ='/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta';
$targets = '/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/sorted_Homo_sapiens.GRCh37.74_nopseudo_noRNA_CDS_picard.bed';
$baits = '/hpc/cog_bioinf/ENRICH/PICARD/sorted_SS_exome_v5_S04380110_Covered_picard.bed';



die "invalid $ref:$!\n" unless -e $ref;

if ( ($targets ne '') or ($baits ne '')) {
    die "invalid $targets:$!\n" unless -e $targets;
    die "invalid $baits:$1\n" unless -e $baits;
}

my $bam = $ARGV[0];
die "$bam does not exists\n" unless -e $bam;


my $pwd = `pwd`;
chomp($pwd);
print "working on $pwd/$bam...\n";
my $command;
my $jid;
my @jid;
my $cfn = $bam;
$cfn =~ s/.bam//;

# CollectInsertSizeMetrics
#$command = $picard."/CollectInsertSizeMetrics.jar R=$ref ASSUME_SORTED=TRUE METRIC_ACCUMULATION_LEVEL=SAMPLE HISTOGRAM_FILE=$pwd/$cfn\_InsertSizeHisto.pdf OUTPUT=$pwd/$cfn\_InsertSizeMetrics.txt INPUT=$pwd/$bam";
#cluster($command,$pwd);

# collect multipleMetrics
$command = $picard."/CollectMultipleMetrics.jar R=$ref ASSUME_SORTED=TRUE OUTPUT=$pwd/$cfn\_MultipleMetrics.txt INPUT=$pwd/$bam PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=QualityScoreDistribution\n";
$command .= "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$cfn\_Picard_Metrics_merged.pdf $pwd/$cfn\_MultipleMetrics.txt.insert_size_histogram.pdf $pwd/$cfn\_MultipleMetrics.txt.quality_by_cycle.pdf $pwd/$cfn\_MultipleMetrics.txt.quality_distribution.pdf";

if (!-e "$pwd/$cfn\_MultipleMetrics.txt.insert_size_histogram.pdf") {
    $jid = cluster($command,$pwd);
    push (@jid,$jid);
}

if ( ($targets ne '') and ($baits ne '') and (!-e "$pwd/$cfn\_HSMetrics.txt" )) {
    #collect enrichment (hybrid NOT amplicon) stats
    $command = $picard."/CalculateHsMetrics.jar R=$ref OUTPUT=$pwd/$cfn\_HSMetrics.txt INPUT=$pwd/$bam BAIT_INTERVALS=$baits TARGET_INTERVALS=$targets METRIC_ACCUMULATION_LEVEL=SAMPLE"; #PER_TARGET_COVERAGE=$pwd/$cfn\_TargetCoverages.txt"; this works only without metric accumulation
    $jid = cluster($command,$pwd);
    push (@jid,$jid);
}


#collect EstimateLibraryComplexity
if (!-e "$pwd/$cfn\_LibComplexity.txt" ) {
    $command = $picard."/EstimateLibraryComplexity.jar OUTPUT=$pwd/$cfn\_LibComplexity.txt INPUT=$pwd/$bam";
    $jid = cluster($command,$pwd);
    push (@jid,$jid);
}

if ((!-e "$pwd/$cfn\_Picard_Metrics_merged.pdf") and (-e "$pwd/$cfn\_MultipleMetrics.txt.insert_size_histogram.pdf")) {
    print " Merging pdfs..\n";
    system "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$cfn\_Picard_Metrics_merged.pdf $pwd/$cfn\_MultipleMetrics.txt.insert_size_histogram.pdf $pwd/$cfn\_MultipleMetrics.txt.quality_by_cycle.pdf $pwd/$cfn\_MultipleMetrics.txt.quality_distribution.pdf";
}





sub cluster {
    my $jid2 = tmpnam();
    $jid2=~s/\/tmp\/file//;
    my $comm = shift;
    my $pwd = shift;
    
    open OUT, ">$pwd/PICARD_$jid2.sh" or die "cannot open file $pwd/PICARD_$jid2.sh\n\n";
    print OUT "#!/bin/bash\n\n";
    print OUT "cd $pwd\n";
    print OUT "$comm\n";
    
    system "qsub -q veryshort -pe threaded 2 -o $pwd/reads -e $pwd/reads -N PICARD_$jid2 $pwd/PICARD_$jid2.sh"; #require two slots for memory reasons
    return $jid2;


}