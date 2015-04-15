#!/usr/bin/perl -w
##################################################################################################################################################
###This script is designed to run the various parts of the illumina 'pipeline' using a single .conf configuration file.
###
###
###Author: S.W.Boymans
###Latest change: ini parsing
###
###TODO:
##################################################################################################################################################

#### Load common perl modules ####
use strict;
use POSIX qw(tmpnam);
use Getopt::Long;
use FindBin;
use File::Path qw(make_path);

#### Load pipeline modules ####
use lib "$FindBin::Bin"; #locates pipeline directory
use illumina_prestats;
use illumina_mapping;
use illumina_poststats;
use illumina_realign;
use illumina_baseRecal;
use illumina_calling;
use illumina_filterVariants;
use illumina_somaticVariants;
use illumina_copyNumber;
use illumina_annotateVariants;
use illumina_check;

### Check correct usage
die usage() if @ARGV == 0;

### initiate opt hash with settings
my %opt;
my $configurationFile;

%opt = (
    'RUNNING_JOBS'		=> {}, #do not use in .conf or .ini
    'BAM_FILES'			=> {}, #do not use in .conf or .ini
    'SAMPLES'			=> undef #do not use in .conf or .ini
);

############ READ RUN SETTINGS FORM .conf FILE ############
$configurationFile = $ARGV[0];

open (CONFIGURATION, "<$configurationFile") or die "Couldn't open .conf file: $configurationFile\n";
while(<CONFIGURATION>){
    chomp;
    next if m/^#/ or ! $_;
    my ($key, $val) = split("\t",$_,2);
    #parse ini file
    if($key eq 'INIFILE') {
	$opt{$key} = $val;
	open (INI, "<$val") or die "Couldn't open .ini file $val\n";
	while(<INI>){
	    chomp;
	    next if m/^#/ or ! $_;
	    my ($key, $val) = split("\t",$_,2);
	    $opt{$key} = $val;
	}
	close INI;
    #parse other config attributes
    } elsif($key eq 'FASTQ' || $key eq 'BAM') {
        $opt{$key}->{$val} = 1;
    } else {
        $opt{$key} = $val;
    }

}
close CONFIGURATION;

############ START PIPELINE  ############

### Check config file
checkConfig();

###Parse samples from FASTQ or BAM files
getSamples();
createOutputDirs();

### Copy ini file to logs dir
system "cp $opt{INIFILE} $opt{OUTPUT_DIR}/logs";

### Start pipeline components
my $opt_ref;

### Mapping or bam input
if(! ($opt{BAM} || $opt{VCF}) ){
    if($opt{PRESTATS} eq "yes"){
	print "###SCHEDULING PRESTATS###\n";
	illumina_prestats::runPreStats(\%opt);
    }

    if($opt{MAPPING} eq "yes"){
	print "\n###SCHEDULING MAPPING###\n";
	$opt_ref = illumina_mapping::runMapping(\%opt);
	%opt = %$opt_ref;
    }

} elsif ( $opt{BAM} ) {
    print "\n###SCHEDULING BAM PREP###\n";
    $opt_ref = illumina_mapping::runBamPrep(\%opt);
    %opt = %$opt_ref;
}

### Post mapping
if (! $opt{VCF} ){
    if($opt{POSTSTATS} eq "yes"){
	print "\n###SCHEDULING POSTSTATS###\n";
	my $postStatsJob = illumina_poststats::runPostStats(\%opt);
	$opt{RUNNING_JOBS}->{'postStats'} = $postStatsJob;
    }

    if($opt{INDELREALIGNMENT} eq "yes"){
	print "\n###SCHEDULING INDELREALIGNMENT###\n";
	$opt_ref = illumina_realign::runRealignment(\%opt);
	%opt = %$opt_ref;
    }

    if($opt{BASEQUALITYRECAL} eq "yes"){
	print "\n###SCHEDULING BASERECALIBRATION###\n";
	$opt_ref = illumina_baseRecal::runBaseRecalibration(\%opt);
	%opt = %$opt_ref;
    }

### Variant Caller
    ### Somatic variant callers
    if($opt{SOMATIC_VARIANTS} eq "yes"){
	print "\n###SCHEDULING SOMATIC VARIANT CALLERS####\n";
	$opt_ref = illumina_somaticVariants::parseSamples(\%opt);
	%opt = %$opt_ref;
	my $somVar_jobs = illumina_somaticVariants::runSomaticVariantCallers(\%opt);
	$opt{RUNNING_JOBS}->{'somVar'} = $somVar_jobs;
    }
    if($opt{COPY_NUMBER} eq "yes"){
	print "\n###SCHEDULING COPY NUMBER TOOLS####\n";
	$opt_ref = illumina_copyNumber::parseSamples(\%opt);
	%opt = %$opt_ref;
	my $cnv_jobs = illumina_copyNumber::runCopyNumberTools(\%opt);
	$opt{RUNNING_JOBS}->{'CNV'} = $cnv_jobs;
    }
    ### GATK
    if($opt{VARIANT_CALLING} eq "yes"){
	print "\n###SCHEDULING VARIANT CALLING####\n";
	$opt_ref = illumina_calling::runVariantCalling(\%opt);
	%opt = %$opt_ref;
    }

} elsif ( $opt{VCF} ) {
    print "\n###RUNNING VCF PREP###\n";
    $opt_ref = illumina_calling::runVcfPrep(\%opt);
    %opt = %$opt_ref;
}

### Filter variants
if($opt{FILTER_VARIANTS} eq "yes"){
    print "\n###SCHEDULING VARIANT FILTRATION####\n";
    my $FVJob = illumina_filterVariants::runFilterVariants(\%opt);
    
    foreach my $sample (@{$opt{SAMPLES}}){
	push (@{$opt{RUNNING_JOBS}->{$sample}} , $FVJob);
    }
}

### Annotate variants
if($opt{ANNOTATE_VARIANTS} eq "yes"){
    print "\n###SCHEDULING VARIANT ANNOTATION####\n";
    my $AVJob = illumina_annotateVariants::runAnnotateVariants(\%opt);
    
    foreach my $sample (@{$opt{SAMPLES}}){
	push (@{$opt{RUNNING_JOBS}->{$sample}} , $AVJob);
    }
}

if($opt{CHECKING} eq "yes"){
    print "\n###SCHEDULING CHECK AND CLEAN####\n";
    illumina_check::runCheck(\%opt);
}

############ SUBROUTINES  ############
sub getSamples{
    my %samples;

    #parse fastq files
    if ($opt{FASTQ}){
	foreach my $input (keys %{$opt{FASTQ}}){
	    my $fastqFile = (split("/", $input))[-1];
	    my $sampleName = (split("_", $fastqFile))[0];
	    $samples{$sampleName} ++;
	    @{$opt{RUNNING_JOBS}->{$sampleName}} = ();
	}
    }

    #parse bam files
    elsif ($opt{BAM}){
	foreach my $input (keys %{$opt{BAM}}){
	    my $bamFile = (split("/", $input))[-1];
	    my $sampleName = $bamFile;
	    $sampleName =~ s/\.bam//g;
	    $samples{$sampleName} ++;
	    @{$opt{RUNNING_JOBS}->{$sampleName}} = ();
	}
    }
    
    @{$opt{SAMPLES}} = keys(%samples);
}

sub createOutputDirs{
    ### Create main output directories
    if(! -e $opt{OUTPUT_DIR}){
	make_path($opt{OUTPUT_DIR}) or die "Couldn't create directory: $opt{OUTPUT_DIR}\n";
    }
    if(! -e "$opt{OUTPUT_DIR}/QCStats"){
	mkdir("$opt{OUTPUT_DIR}/QCStats") or die "Couldn't create directory: $opt{OUTPUT_DIR}/QCStats\n";
    }
    if(! -e "$opt{OUTPUT_DIR}/jobs"){
	mkdir("$opt{OUTPUT_DIR}/jobs") or die "Couldn't create directory: $opt{OUTPUT_DIR}/jobs\n";
    }
    if(! -e "$opt{OUTPUT_DIR}/logs"){
	mkdir("$opt{OUTPUT_DIR}/logs") or die "Couldn't create directory: $opt{OUTPUT_DIR}/logs\n";
    }
    if(! -e "$opt{OUTPUT_DIR}/tmp"){
	mkdir("$opt{OUTPUT_DIR}/tmp") or die "Couldn't create directory: $opt{OUTPUT_DIR}/tmp\n";
    }

    ### Create sample specific output directories
    foreach my $sample (@{$opt{SAMPLES}}){
	if(! -e "$opt{OUTPUT_DIR}/$sample"){
	    mkdir("$opt{OUTPUT_DIR}/$sample") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sample\n";
	}
	if(! -e "$opt{OUTPUT_DIR}/$sample/mapping"){
	    mkdir("$opt{OUTPUT_DIR}/$sample/mapping") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sample/mapping\n";
	}
	if(! -e "$opt{OUTPUT_DIR}/$sample/QCStats"){
	    mkdir("$opt{OUTPUT_DIR}/$sample/QCStats") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sample/QCStats\n";
	}
	if(! -e "$opt{OUTPUT_DIR}/$sample/jobs"){
	    mkdir("$opt{OUTPUT_DIR}/$sample/jobs") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sample/jobs\n";
	}
	if(! -e "$opt{OUTPUT_DIR}/$sample/logs"){
	    mkdir("$opt{OUTPUT_DIR}/$sample/logs") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sample/logs\n";
	}
	if(! -e "$opt{OUTPUT_DIR}/$sample/tmp"){
	    mkdir("$opt{OUTPUT_DIR}/$sample/tmp") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sample/tmp\n";
	}
    }
}

sub usage{
    warn <<END;
    Usage: perl illumina_pipeline.pl configurationFile.conf
END
    exit;
}

sub get_job_id {
    my $id = tmpnam();
    $id =~ s/\/tmp\/file//;
    return $id;
}

sub checkConfig{
    my $checkFailed = 0;
    ### Input and Output
    if(! $opt{INIFILE}){ print "ERROR: No INIFILE option found in config files.\n"; $checkFailed = 1; }
    if(! $opt{OUTPUT_DIR}){ print "ERROR: No OUTPUT_DIR found in config files.\n"; $checkFailed = 1; }
    if(! ($opt{FASTQ} || $opt{BAM} || $opt{VCF}) ){ print "ERROR: No FASTQ/BAM/VCF files found in config files.\n"; $checkFailed = 1; }
    if(! $opt{MAIL}){ print "ERROR: No MAIL address specified in config files.\n"; $checkFailed = 1; }

    ### Cluster settings
    if(! $opt{CLUSTER_PATH}){ print "ERROR: No CLUSTER_PATH option found in config files.\n"; $checkFailed = 1; }
    if(! $opt{CLUSTER_TMP}){ print "ERROR: No CLUSTER_TMP option found in config files.\n"; $checkFailed = 1; }
    if(! $opt{CLUSTER_RESERVATION}){ print "ERROR: No CLUSTER_RESERVATION option found in config files.\n"; $checkFailed = 1; }

    ### Module yes or No
    if(! $opt{PRESTATS}){ print "ERROR: No PRESTATS option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{MAPPING}){ print "ERROR: No MAPPING option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{POSTSTATS}){ print "ERROR: No POSTSTATS option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{INDELREALIGNMENT}){ print "ERROR: No INDELREALIGNMENT option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{BASEQUALITYRECAL}){ print "ERROR: No BASEQUALITYRECAL option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{VARIANT_CALLING}){ print "ERROR: No VARIANT_CALLING option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{FILTER_VARIANTS}){ print "ERROR: No FILTER_VARIANTS option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{SOMATIC_VARIANTS}){ print "ERROR: No SOMATIC_VARIANTS option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{COPY_NUMBER}){ print "ERROR: No COPY_NUMBER option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{ANNOTATE_VARIANTS}){ print "ERROR: No ANNOTATE_VARIANTS option found in config files. \n"; $checkFailed = 1; }
    if(! $opt{CHECKING}){ print "ERROR: No CHECKING option found in config files. \n"; $checkFailed = 1; }

    ### Module Settings / tools
    if(! $opt{GENOME}){ print "ERROR: No GENOME option found in config files.\n"; $checkFailed = 1; }
    elsif(! -e $opt{GENOME}){ print"ERROR: $opt{GENOME} does Not exist\n"}
    if(! $opt{SAMBAMBA_PATH}){ print "ERROR: No SAMBAMBA_PATH option found in config files.\n"; $checkFailed = 1; }
    if(! $opt{QUEUE_PATH}){ print "ERROR: No PICARD_PATH option found in config files.\n"; $checkFailed = 1; }
    ## PRESTATS
    if($opt{PRESTATS} eq "yes"){
	if(! $opt{FASTQC_PATH}){ print "ERROR: No FASTQC_PATH option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{PRESTATS_THREADS}){ print "ERROR: No PRESTATS_THREADS option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{PRESTATS_MEM}){ print "ERROR: No PRESTATS_MEM option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{PRESTATS_QUEUE}){ print "ERROR: No PRESTATS_QUEUE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{PRESTATS_PROJECT}){ print "ERROR: No PRESTATS_PROJECT option found in config files.\n"; $checkFailed = 1; }
    }
    ## MAPPING
    if($opt{MAPPING} eq "yes"){
	if(! $opt{BWA_PATH}){ print "ERROR: No BWA_PATH option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{MAPPING_THREADS}){ print "ERROR: No MAPPING_THREADS option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{MAPPING_MEM}){ print "ERROR: No MAPPING_MEM option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{MAPPING_QUEUE}){ print "ERROR: No MAPPING_QUEUE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{MAPPING_PROJECT}){ print "ERROR: No MAPPING_PROJECT option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{MAPPING_MODE}){ print "ERROR: No MAPPING_MODE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{MAPPING_MARKDUP}){ print "ERROR: No MAPPING_MARKDUP option found in config files.\n"; $checkFailed = 1; }
	if($opt{MAPPING_MODE} eq 'single'){
	    if(! $opt{FLAGSTAT_QUEUE}){ print "ERROR: No FLAGSTAT_QUEUE option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{FLAGSTAT_THREADS}){ print "ERROR: No FLAGSTAT_THREADS option found in config files.\n"; $checkFailed = 1; }
	}
    }
    ## POSTSTATS
    if($opt{POSTSTATS} eq "yes"){
	if(! $opt{PICARD_PATH}){ print "ERROR: No PICARD_PATH option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{POSTSTATS_THREADS}){ print "ERROR: No POSTSTATS_THREADS option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{POSTSTATS_MEM}){ print "ERROR: No POSTSTATS_MEM option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{POSTSTATS_QUEUE}){ print "ERROR: No POSTSTATS_THREADS option found in config files.\n"; $checkFailed = 1; }
	if(! ($opt{POSTSTATS_TARGETS}) && ! ($opt{POSTSTATS_BAITS}) ){
	    if(! $opt{POSTSTATS_COVERAGECAP}){ print "ERROR: No POSTSTATS_COVERAGECAP or (POSTSTATS_TARGETS & POSTSTATS_BAITS) options found in config files.\n"; $checkFailed = 1; }
	}
	if( $opt{POSTSTATS_TARGETS} && ! -e $opt{POSTSTATS_TARGETS}){ print "ERROR: $opt{POSTSTATS_TARGETS} does Not exist\n"; $checkFailed = 1; }
	if( $opt{POSTSTATS_BAITS} && ! -e $opt{POSTSTATS_BAITS}){ print "ERROR: $opt{POSTSTATS_BAITS} does Not exist\n"; $checkFailed = 1; }
    }
    ## INDELREALIGNMENT
    if($opt{INDELREALIGNMENT} eq "yes"){
	if(! $opt{REALIGNMENT_MASTERQUEUE}){ print "ERROR: No REALIGNMENT_MASTERQUEUE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{REALIGNMENT_MASTERTHREADS}){ print "ERROR: No REALIGNMENT_MASTERTHREADS option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{REALIGNMENT_QUEUE}){ print "ERROR: No REALIGNMENT_QUEUE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{REALIGNMENT_PROJECT}){ print "ERROR: No REALIGNMENT_PROJECT option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{REALIGNMENT_THREADS}){ print "ERROR: No REALIGNMENT_THREADS option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{REALIGNMENT_MERGETHREADS}){ print "ERROR: No REALIGNMENT_MERGETHREADS option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{REALIGNMENT_MEM}){ print "ERROR: No REALIGNMENT_MEM option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{REALIGNMENT_SCALA}){ print "ERROR: No REALIGNMENT_SCALA option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{REALIGNMENT_SCATTER}){ print "ERROR: No REALIGNMENT_SCATTER option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{REALIGNMENT_MODE}){ print "ERROR: No REALIGNMENT_MODE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{QUEUE_RETRY}){ print "ERROR: No QUEUE_RETRY option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{FLAGSTAT_QUEUE}){ print "ERROR: No FLAGSTAT_QUEUE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{FLAGSTAT_THREADS}){ print "ERROR: No FLAGSTAT_THREADS option found in config files.\n"; $checkFailed = 1; }
    }
    ## BASEQUALITYRECAL
    if($opt{BASEQUALITYRECAL} eq "yes"){
	if(! $opt{BASERECALIBRATION_MASTERQUEUE}){ print "ERROR: No BASERECALIBRATION_QUEUE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{BASERECALIBRATION_MASTERTHREADS}){ print "ERROR: No BASERECALIBRATION_THREADS option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{BASERECALIBRATION_QUEUE}){ print "ERROR: No BASERECALIBRATION_QUEUE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{BASERECALIBRATION_THREADS}){ print "ERROR: No BASERECALIBRATION_THREADS option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{BASERECALIBRATION_MEM}){ print "ERROR: No BASERECALIBRATION_MEM option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{BASERECALIBRATION_SCALA}){ print "ERROR: No BASERECALIBRATION_SCALA option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{BASERECALIBRATION_SCATTER}){ print "ERROR: No BASERECALIBRATION_SCATTER option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{QUEUE_RETRY}){ print "ERROR: No QUEUE_RETRY option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{FLAGSTAT_QUEUE}){ print "ERROR: No FLAGSTAT_QUEUE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{FLAGSTAT_THREADS}){ print "ERROR: No FLAGSTAT_THREADS option found in config files.\n"; $checkFailed = 1; }
    }
    ## VARIANT_CALLING
    if($opt{VARIANT_CALLING} eq "yes"){
	if(! $opt{CALLING_MASTERQUEUE}){ print "ERROR: No CALLING_MASTERQUEUE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{CALLING_MASTERTHREADS}){ print "ERROR: No CALLING_MASTERTHREADS option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{CALLING_QUEUE}){ print "ERROR: No CALLING_QUEUE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{CALLING_THREADS}){ print "ERROR: No CALLING_THREADS option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{CALLING_MEM}){ print "ERROR: No CALLING_QUEUE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{CALLING_SCATTER}){ print "ERROR: No CALLING_SCATTER option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{CALLING_GVCF}){ print "ERROR: No CALLING_GVCF option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{CALLING_SCALA}){ print "ERROR: No CALLING_SCALA option found in config files.\n"; $checkFailed = 1; }
	if($opt{CALLING_UGMODE}){ 
	    if($opt{CALLING_UGMODE} ne "SNP" and $opt{CALLING_UGMODE} ne "INDEL" and $opt{CALLING_UGMODE} ne "BOTH"){ print "ERROR: UGMODE: $opt{CALLING_UGMODE} does Not exist use SNP, INDEL or BOTH\n"; $checkFailed = 1; }
	}
	if(! $opt{CALLING_STANDCALLCONF}){ print "ERROR: No CALLING_STANDCALLCONF option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{CALLING_STANDEMITCONF}){ print "ERROR: No CALLING_STANDEMITCONF option found in config files.\n"; $checkFailed = 1; }
	if( $opt{CALLING_TARGETS} && ! -e $opt{CALLING_TARGETS}) { print"ERROR: $opt{CALLING_TARGETS} does Not exist\n"; $checkFailed = 1; }
	if( $opt{CALLING_DBSNP} && ! -e $opt{CALLING_DBSNP}) { print"ERROR: $opt{CALLING_DBSNP} does Not exist\n"; $checkFailed = 1; }
	if(! $opt{QUEUE_RETRY}){ print "ERROR: No QUEUE_RETRY option found in config files.\n"; $checkFailed = 1; }
    }
    ## FILTER_VARIANTS
    if($opt{FILTER_VARIANTS} eq "yes"){
	if(! $opt{FILTER_MASTERQUEUE}){ print "ERROR: No FILTER_MASTERQUEUE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{FILTER_MASTERTHREADS}){ print "ERROR: No FILTER_MASTERTHREADS option found in config files.\n"; $checkFailed = 1; }	
	if(! $opt{FILTER_QUEUE}){ print "ERROR: No FILTER_QUEUE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{FILTER_THREADS}){ print "ERROR: No FILTER_THREADS option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{FILTER_MEM}){ print "ERROR: No FILTER_QUEUE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{FILTER_SCATTER}){ print "ERROR: No FILTER_SCATTER option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{FILTER_SCALA}){ print "ERROR: No FILTER_SCALA option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{FILTER_MODE}){ print "ERROR: No FILTER_MODE  option found in config files.\n"; $checkFailed = 1; }
	if($opt{FILTER_MODE} ne "SNP" and $opt{FILTER_MODE} ne "INDEL" and $opt{FILTER_MODE} ne "BOTH"){ print "ERROR: FILTER_MODE $opt{FILTER_MODE} does Not exist use SNP, INDEL or BOTH\n"; $checkFailed = 1; }
	if ($opt{FILTER_MODE} eq "SNP" || $opt{FILTER_MODE} eq "BOTH") {
	    if(! $opt{FILTER_SNPNAME}){ print "ERROR: No FILTER_SNPNAME option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{FILTER_SNPEXPR}){ print "ERROR: No FILTER_SNPEXPR  option found in config files.\n"; $checkFailed = 1; }
	}
	if ($opt{FILTER_MODE} eq "INDEL" || $opt{FILTER_MODE} eq "BOTH") {
	    if(! $opt{FILTER_INDELNAME}){ print "ERROR: No FILTER_INDELNAME option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{FILTER_INDELEXPR}){ print "ERROR: No FILTER_INDELEXPR option found in config files.\n"; $checkFailed = 1; }
	}
	if(! $opt{QUEUE_RETRY}){ print "ERROR: No QUEUE_RETRY option found in config files.\n"; $checkFailed = 1; }
    }
    ## SOMATIC_VARIANTS
    if($opt{SOMATIC_VARIANTS} eq "yes"){
	if(! $opt{SAMTOOLS_PATH}){ print "ERROR: No SAMTOOLS_PATH option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{SOMVAR_TARGETS}){ print "ERROR: No SOMVAR_TARGETS option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{SOMVAR_STRELKA}){ print "ERROR: No SOMVAR_STRELKA option found in config files.\n"; $checkFailed = 1; }
	if($opt{SOMVAR_STRELKA} eq "yes"){
	    if(! $opt{STRELKA_PATH}){ print "ERROR: No STRELKA_PATH option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{STRELKA_INI}){ print "ERROR: No STRELKA_INI option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{STRELKA_QUEUE}){ print "ERROR: No STRELKA_QUEUE option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{STRELKA_THREADS}){ print "ERROR: No STRELKA_THREADS option found in config files.\n"; $checkFailed = 1; }
	}
	if(! $opt{SOMVAR_VARSCAN}){ print "ERROR: No SOMVAR_VARSCAN option found in config files.\n"; $checkFailed = 1; }
	if($opt{SOMVAR_VARSCAN} eq "yes"){
	    if(! $opt{VARSCAN_PATH}){ print "ERROR: No VARSCAN_PATH option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{VARSCAN_QUEUE}){ print "ERROR: No VARSCAN_QUEUE option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{VARSCAN_THREADS}){ print "ERROR: No VARSCAN_THREADS option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{VARSCAN_SETTINGS}){ print "ERROR: No VARSCAN_SETTINGS option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{VARSCAN_POSTSETTINGS}){ print "ERROR: No VARSCAN_POSTSETTINGS option found in config files.\n"; $checkFailed = 1; }
	}
	if(! $opt{SOMVAR_FREEBAYES}){ print "ERROR: No SOMVAR_FREEBAYES option found in config files.\n"; $checkFailed = 1; }
	if($opt{SOMVAR_FREEBAYES} eq "yes"){
	    if(! $opt{FREEBAYES_PATH}){ print "ERROR: No FREEBAYES_PATH option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{VCFSAMPLEDIFF_PATH}){ print "ERROR: No VCFSAMPLEDIFF_PATH option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{BIOVCF_PATH}){ print "ERROR: No BIOVCF_PATH option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{FREEBAYES_QUEUE}){ print "ERROR: No FREEBAYES_QUEUE option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{FREEBAYES_THREADS}){ print "ERROR: No FREEBAYES_THREADS option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{FREEBAYES_SETTINGS}){ print "ERROR: No FREEBAYES_SETTINGS option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{FREEBAYES_SOMATICFILTER}){ print "ERROR: No FREEBAYES_SOMATICFILTER option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{FREEBAYES_GERMLINEFILTER}){ print "ERROR: No FREEBAYES_GERMLINEFILTER option found in config files.\n"; $checkFailed = 1; }
	}
	if(! $opt{SOMVARMERGE_QUEUE}){ print "ERROR: No SOMVARMERGE_QUEUE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{SOMVARMERGE_THREADS}){ print "ERROR: No SOMVARMERGE_THREADS option found in config files.\n"; $checkFailed = 1; }
    }
    ## COPY_NUMBER
    if($opt{COPY_NUMBER} eq "yes"){
	if(! $opt{CNVCHECK_QUEUE} ) { print "ERROR: No CNVCHECK_QUEUE in config files.\n"; $checkFailed = 1; }
	if(! $opt{CNVCHECK_THREADS} ) { print "ERROR: No CNVCHECK_THREADS  in config files.\n"; $checkFailed = 1; }
	if(! $opt{CNV_CONTRA}){ print "ERROR: No CNV_CONTRA  in config files.\n"; $checkFailed = 1; }
	if($opt{CNV_CONTRA} eq "yes"){
	    if(! $opt{CONTRA_PATH}){ print "ERROR: No CONTRA_PATH option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{CONTRA_QUEUE}){ print "ERROR: No CONTRA_QUEUE option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{CONTRA_THREADS}){ print "ERROR: No CONTRA_THREADS option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{CONTRA_TARGETS}){ print "ERROR: No CONTRA_TARGETS option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{CONTRA_FLAGS}){ print "ERROR: No CONTRA_FLAGS option found in config files.\n"; $checkFailed = 1; }
	}
    }
    ## ANNOTATE_VARIANTS
    if($opt{ANNOTATE_VARIANTS} eq "yes"){
	if(! $opt{SNPEFF_PATH}){ print "ERROR: No SNPEFF_PATH option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{IGVTOOLS_PATH}){ print "ERROR: No IGVTOOLS_PATH option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{ANNOTATE_QUEUE}){ print "ERROR: No ANNOTATE_QUEUE option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{ANNOTATE_THREADS}){ print "ERROR: No ANNOTATE_THREADS option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{ANNOTATE_MEM}){ print "ERROR: No ANNOTATE_MEM option found in config files.\n"; $checkFailed = 1; }
	if(! $opt{ANNOTATE_SNPEFF}){ print "ERROR: No ANNOTATE_SNPEFF option found in config files.\n"; $checkFailed = 1; }
	if($opt{ANNOTATE_SNPEFF} eq "yes"){
	    if(! $opt{ANNOTATE_DB}){ print "ERROR: No ANNOTATE_DB option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{ANNOTATE_FLAGS}){ print "ERROR: No ANNOTATE_FLAGS option found in config files.\n"; $checkFailed = 1; }
	}
	if(! $opt{ANNOTATE_SNPSIFT}){ print "ERROR: No ANNOTATE_SNPSIFT option found in config files.\n"; $checkFailed = 1; }
	if($opt{ANNOTATE_SNPSIFT} eq "yes"){
	    if(! $opt{ANNOTATE_DBNSFP}){ print "ERROR: No ANNOTATE_DBNSFP option found in config files.\n"; $checkFailed = 1; }
	    elsif( $opt{ANNOTATE_DBNSFP} && ! -e $opt{ANNOTATE_DBNSFP}) { print"ERROR: $opt{ANNOTATE_DBNSFP} does Not exist\n"; $checkFailed = 1; }
	    if(! $opt{ANNOTATE_FIELDS}){ print "ERROR: No ANNOTATE_FIELDS option found in config files.\n"; $checkFailed = 1; }
	}
	if(! $opt{ANNOTATE_FREQUENCIES}){ print "ERROR: No ANNOTATE_FREQUENCIES option found in config files.\n"; $checkFailed = 1; }
	if($opt{ANNOTATE_FREQUENCIES} eq "yes"){
	    if(! $opt{ANNOTATE_FREQNAME}){ print "ERROR: No ANNOTATE_FREQNAME option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{ANNOTATE_FREQDB}){ print "ERROR: No ANNOTATE_FREQDB option found in config files.\n"; $checkFailed = 1; }
	    elsif( $opt{ANNOTATE_FREQDB} && ! -e $opt{ANNOTATE_FREQDB}) { print"ERROR: $opt{ANNOTATE_FREQDB} does Not exist\n"; $checkFailed = 1; }
	    if(! $opt{ANNOTATE_FREQINFO}){ print "ERROR: No ANNOTATE_FREQINFO option found in config files.\n"; $checkFailed = 1; }
	}
	if(! $opt{ANNOTATE_IDFIELD}){ print "ERROR: No ANNOTATE_IDFIELD option found in config files.\n"; $checkFailed = 1; }
	if($opt{ANNOTATE_IDFIELD} eq "yes"){
	    if(! $opt{ANNOTATE_IDNAME}){ print "ERROR: No ANNOTATE_IDNAME option found in config files.\n"; $checkFailed = 1; }
	    if(! $opt{ANNOTATE_IDDB}){ print "ERROR: No ANNOTATE_IDDB option found in config files.\n"; $checkFailed = 1; }
	}
    }
    if ($checkFailed) { die "One or more options not found in config files." }
}

1;