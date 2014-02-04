#!/usr/bin/perl -w
use strict;
use POSIX qw(tmpnam);

my $HCscalaScript = "/hpc/cog_bioinf/common_scripts/GATK_v2.7/ies/discoverVariants_HaplotypeCaller_multiSample_noRef.scala";
my $HFscalaScript = "/hpc/cog_bioinf/common_scripts/GATK_v2.7/ies/discoverVariants_hardfilter_ies.scala";

my $inputReference="/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta";
my $dbsnpFile="/hpc/cog_bioinf/common_scripts/GATK_v2.7/bundle/dbsnp_137.b37.vcf";
my $targetList = "/hpc/cog_bioinf/ENRICH/GATK/SS_exome_v5_S04380110_Covered_GATK.interval_list";


my $email = "i.nijman\@umcutrecht.nl";


my $pwd = `pwd`;
chomp($pwd);
my $path = $pwd;
my @dirs = split(/\//,$path);

my $outputCoreName = $dirs[-1];

print "outputting to $outputCoreName\n"; 

die "nothing to do..Forgot to supply BAMS?\n" if scalar(@ARGV) < 1;

print "selected BAMS:\n";
print "$_\n" foreach @ARGV;

my $bam_string;
$bam_string  .= " -I $_ " foreach @ARGV;

#print $bam_string,"\n" ;

die "bam files invalid\n" if $bam_string =~ /\*/;

my $command;

#run Haplotype caller
$command = "java -Xmx22G -Xms2G -jar /hpc/cog_bioinf/common_scripts/GATK_v2.7/Queue.jar -jobRunner GridEngine -jobQueue veryshort -jobEnv \"threaded 4\" -S $HCscalaScript -R $inputReference -D $dbsnpFile -O $outputCoreName -L $targetList $bam_string -run\n";

#add hard filters
$command .= "java -Xmx12G -Xms2G -jar /hpc/cog_bioinf/common_scripts/GATK_v2.7/Queue.jar -jobRunner GridEngine -jobQueue veryshort -jobEnv \"threaded 4\" -S $HFscalaScript -R $inputReference -V $outputCoreName.raw_variants.vcf -O $outputCoreName -run\n";

#cleanup
$command .="rm *SNPS.vcf* *INDELS.vcf*\n";
$command .="bgzip *.vcf\n";
$command .="tabix *.vcf.gz -p vcf\n";


# execute
cluster($command,$pwd);


# call vcf filter



#===========================================================================================================================
sub cluster {
    my $jid2 = tmpnam();
    $jid2=~s/\/tmp\/file//;
    my $comm = shift;
    my $pwd = shift;
    
    open OUT, ">$pwd/HC_$jid2.sh" or die "cannot open file $pwd/HC_$jid2.sh\n\n";
    print OUT "#!/bin/bash\n\n";
    print OUT "cd $pwd\n";
    print OUT "$comm\n";
    
    system "qsub -q short -pe threaded 2 -M $email -m beas -o $pwd -e $pwd -N HC_$jid2 $pwd/HC_$jid2.sh"; 
    return $jid2;


}