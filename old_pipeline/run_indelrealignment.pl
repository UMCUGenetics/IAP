#!/usr/bin/perl -w
use strict;

system "cd \$PWd";

my $q = 'veryshort';
$q = $ARGV[0] if $ARGV[0];


print "submitting to queue $q use first parameter after script to override queue\n";

my $scalaScript="/hpc/cog_bioinf/common_scripts/GATK_v2.7/ies/preProcessing_indelRealign.scala";
my $inputKnownIndelFile="/hpc/cog_bioinf/common_scripts/GATK_v2.7/bundle/Mills_and_1000G_gold_standard.indels.b37.vcf -known /hpc/cog_bioinf/common_scripts/GATK_v2.7/bundle/Mills_and_1000G_gold_standard.indels.b37.vcf ";
my $inputReference="/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta";
#my $inputReference="/hpc/cog_bioinf/GENOMES/rat_GATK_illumina_rnor_50/Rn_Rn05_ill_gatk_sorted.fasta";


foreach my $dir (<*>) {
    next unless -d $dir;
    chdir $dir;
    my $run = 0;
    
    if (-e "$dir\_MERGED_sorted_dedup.realigned.bam") {
	print "realignment already done for $dir, skipping...\n";
	chdir '../';
	next;
    }
    foreach my $bam (<*MERGED*.bam>) {
	my $run = 1;
	print "working on $bam...\n";
	my $outputCoreName = $bam;
	$outputCoreName =~ s/.bam//;
	
	next if (-e "$outputCoreName\.realigned.bam");
	
	
	
	my $pwd = `pwd`;
	chomp($pwd);
	open OUT, ">reads/REALIGN.sh" or die "cannot open realign.sh:$!\n";
	print OUT "#!/bin/bash\n\n";
	print OUT "cd $pwd\n";
	print OUT "java -Xmx5G -Xms2G -jar /hpc/cog_bioinf/common_scripts/GATK_v2.7/Queue.jar -jobRunner GridEngine -jobQueue $q -jobEnv \"threaded 1\" -S $scalaScript -I $pwd/$bam -known $inputKnownIndelFile -R $inputReference -O $pwd/$outputCoreName -run\n";
	close OUT;
	
	system "qsub -q short -M inijman\@umcutrecht.nl -m ase -o $pwd/reads/realign_out -e $pwd/reads/realign_err -N realign_Q reads/REALIGN.sh";
    }
    chdir '../';
}

