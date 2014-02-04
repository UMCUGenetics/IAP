#!/bin/perl -w
use strict;
use POSIX qw(tmpnam);


#collect several picard metrics files
#	single bams files
#

my $picard = "java -Xmx16G -jar /hpc/cog_bioinf/common_scripts/picard-tools-1.98";
my $targets = '';
my $baits ='';


my $ref = "/hpc/cog_bioinf/GENOMES/rat_GATK_illumina_rnor_50/Rn_Rn05_ill_gatk_sorted.fasta";
my $targets = '/hpc/cog_bioinf/ENRICH/PICARD/UCSC_coords_Rno5_full_chr_regios.5col.bed'; #ie exome
my $baits = '/hpc/cog_bioinf/ENRICH/PICARD/OID41455_Rn05_Match5_probes_loc_5col.bed'; #ie sureselect

#my $ref ='/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta';
#$targets = '/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/sorted_Homo_sapiens.GRCh37.74_nopseudo_noRNA_CDS_picard.bed';
#$baits = '/hpc/cog_bioinf/ENRICH/PICARD/sorted_SS_exome_v5_S04380110_Covered_picard.bed';



die "invalid $ref:$!\n" unless -e $ref;

if ( ($targets ne '') or ($baits ne '')) {
    die "invalid $targets:$!\n" unless -e $targets;
    die "invalid $baits:$1\n" unless -e $baits;
}



#BAIT_INTERVALS=/hpc/cog_bioinf/ENRICH/PICARD/OID41455_Rn05_Match5_probes_loc_5col.bed TARGET_INTERVALS=/hpc/cog_bioinf/ENRICH/PICARD/UCSC_coords_Rno5_full_chr_regios.5col.bed 


system "cd \$PWD";

foreach my $dir (<*>) {
    next unless -d $dir;
    chdir $dir;
    foreach my $bam (<*MERGED*.bam>) {
	my $pwd = `pwd`;
	chomp($pwd);
	print "working on $pwd/$bam...\n";
	my $command;
	my $cfn = $bam;
	$cfn =~ s/.bam//;
	
	# CollectInsertSizeMetrics
	#$command = $picard."/CollectInsertSizeMetrics.jar R=$ref ASSUME_SORTED=TRUE METRIC_ACCUMULATION_LEVEL=SAMPLE HISTOGRAM_FILE=$pwd/$cfn\_InsertSizeHisto.pdf OUTPUT=$pwd/$cfn\_InsertSizeMetrics.txt INPUT=$pwd/$bam";
	#cluster($command,$pwd);
	
	# collect multipleMetrics
	$command = $picard."/CollectMultipleMetrics.jar R=$ref ASSUME_SORTED=TRUE OUTPUT=$pwd/$cfn\_MultipleMetrics.txt INPUT=$pwd/$bam PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=QualityScoreDistribution";
	cluster($command,$pwd);
	
	if ( ($targets ne '') and ($baits ne '')) {
	    #collect enrichment (hybrid NOT amplicon) stats
	    $command = $picard."/CalculateHsMetrics.jar R=$ref OUTPUT=$pwd/$cfn\_HSMetrics.txt INPUT=$pwd/$bam BAIT_INTERVALS=$baits TARGET_INTERVALS=$targets METRIC_ACCUMULATION_LEVEL=SAMPLE"; #PER_TARGET_COVERAGE=$pwd/$cfn\_TargetCoverages.txt"; this works only without metric accumulation
	    cluster($command,$pwd);
	}
	
	
	#collect EstimateLibraryComplexity
	$command = $picard."/EstimateLibraryComplexity.jar OUTPUT=$pwd/$cfn\_LibComplexity.txt INPUT=$pwd/$bam";
	cluster($command,$pwd);
    }
    chdir "../";
}



sub cluster {
    my $jid = tmpnam();
    $jid=~s/\/tmp\/file//;
    my $comm = shift;
    my $pwd = shift;
    
    open OUT, ">$pwd/PRE_$jid.sh" or die "cannot open file $pwd/PRE_$jid.sh\n\n";
    print OUT "#!/bin/bash\n\n";
    print OUT "cd $pwd\n";
    print OUT "$comm\n";
    
    system "qsub -q veryshort -pe threaded 2 -o $pwd/reads -e $pwd/reads -N PRE_$jid $pwd/PRE_$jid.sh"; #require two slots for memory reasons


}