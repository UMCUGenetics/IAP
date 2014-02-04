#!/usr/bin/perl -w
use strict;
use POSIX qw(tmpnam);

my $inputReference="/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta";
my $dbsnpFile="/hpc/cog_bioinf/common_scripts/GATK_v2.7/bundle/dbsnp_137.b37.vcf";
my $targetList = "/hpc/cog_bioinf/ENRICH/GATK/SS_exome_v5_S04380110_Covered_GATK.interval_list";


my $email = "i.nijman\@umcutrecht.nl";
my $snpEff_path = "/hpc/cog_bioinf/data/ies/src/snpEff";
my $snpEff_db = "GRCh37.74";
my $dbNSFP = "/hpc/cog_bioinf/data/ies/src/snpEff/data/dbNSFP2.1.txt";


my $fn_file = "/hpc/cog_bioinf/data/ies/src/snpEff/data/annotation_fields_ies.txt";
if ($ARGV[1]) {
    $fn_file = $ARGV[1] if -e $ARGV[1];
}



my $vcf = $ARGV[0];
die "vcf file invalid or non-existent...\n" unless -e $vcf;

my ($field_names,$command);

open FL, $fn_file;
warn "field_names file invalid or non-existent...using default $fn_file\n" unless -e $fn_file;


my $pwd = `pwd`;
chomp($pwd);
my $path = $pwd;

my $jid = tmpnam();
$jid=~s/\/tmp\/file//;


while (my $field=<FL>) {
    chomp($field);
    $field =~ s/\r//;
    $field_names .= "$field,";
}

my $gonl_freq = "/hpc/cog_bioinf/common_dbs/GoNL/GoNL_freq_annotation2.txt.gz";
die "no GoNl freq file\n" unless -e $gonl_freq;

my $gonl_freq_descr = "/hpc/cog_bioinf/common_dbs/GoNL/GoNL_freq_annotation_descriptions.txt";
die "no GoNl freq descr file\n" unless -e $gonl_freq_descr;

print "Going to annotate $vcf with the following fields:\n";
print $field_names,"\n\n";

open OUT, ">ANNOTATE\_$jid.sh" or die "cannot create annotate.sh file: $!\n";


my $outvcf = $vcf;
my $sift_vcf = $outvcf;
$outvcf =~ s/.vcf$/_snpEff.vcf/;
$outvcf =~ s/.vcf.gz$/_snpEff.vcf/;

$sift_vcf =~ s/.vcf/_Sift.vcf/; 
$sift_vcf =~ s/.gz$//;

my $gonl_vcf = $sift_vcf;
$gonl_vcf =~ s/_Sift/_Sift_GoNL/;



print OUT "#!/bin/bash\n\n";
print OUT "cd $pwd\n";

#run basic eff predication and annotation
$command = "java -Xmx5g -jar $snpEff_path/snpEff.jar -c $snpEff_path/snpEff.config $snpEff_db -v $vcf -hgvs -lof -o gatk -no-downstream -no-intergenic -no-upstream > $outvcf";
print OUT $command,"\n";

#run detailed annotation (sift, polyphen gerp, phyloP
$command = "java -Xmx5g -jar $snpEff_path/SnpSift.jar dbnsfp -f $field_names -v $dbNSFP $outvcf > $sift_vcf";
print OUT $command,"\n";

print OUT "if [ -f $sift_vcf ]\nthen\n rm $outvcf\nfi\n\n";


#Add GoNL annotation
$command = "vcf-annotate -a $gonl_freq -c CHROM,FROM,TO,-,GoNL_CHR,GoNL_FREQ1,GoNL_FREQ2 -d $gonl_freq_descr $sift_vcf > $gonl_vcf";
print OUT $command,"\n";

print OUT "if [ -f $gonl_vcf ]\nthen\n rm $sift_vcf\nfi\n\n";

#system "time sh ANNOTATE.sh";
system "qsub -q veryshort -o $pwd -e $pwd -N ANNOTATE_$jid ANNOTATE_$jid.sh";