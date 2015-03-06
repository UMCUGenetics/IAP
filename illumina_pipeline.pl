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
      $id=~s/\/tmp\/file//;
   return $id;
}

sub checkConfig{
    ### Input and Output
    if(! $opt{INIFILE}){ die "ERROR: No INIFILE found in .conf file\n" }
    if(! $opt{OUTPUT_DIR}){ die "ERROR: No OUTPUT_DIR found in .conf file\n" }
    if(! ($opt{FASTQ} || $opt{BAM} || $opt{VCF}) ){ die "ERROR: No FASTQ/BAM/VCF files specified\n" }
    if(! $opt{MAIL}){ die "ERROR: No MAIL address specified in .conf file\n" }
    
    ### Cluster settings
    if(! $opt{CLUSTER_PATH}){ die "ERROR: No CLUSTER_PATH found in .ini file\n" }
    if(! $opt{CLUSTER_TMP}){ die "ERROR: No CLUSTER_TMP found in .ini file\n" }
    if(! $opt{CLUSTER_RESERVATION}){ die "ERROR: No CLUSTER_RESERVATION found in .ini file\n" }
    
    ### Module yes or no
    if(! $opt{PRESTATS}){ die "ERROR: No PRESTATS option in .conf file \n" }
    if(! $opt{MAPPING}){ die "ERROR: No MAPPING option in .conf file \n" }
    if(! $opt{POSTSTATS}){ die "ERROR: No POSTSTATS option in .conf file \n" }
    if(! $opt{INDELREALIGNMENT}){ die "ERROR: No INDELREALIGNMENT option in .conf file \n" }
    if(! $opt{BASEQUALITYRECAL}){ die "ERROR: No BASEQUALITYRECAL option in .conf file \n" }
    if(! $opt{VARIANT_CALLING}){ die "ERROR: No VARIANT_CALLING option in .conf file \n" }
    if(! $opt{FILTER_VARIANTS}){ die "ERROR: No FILTER_VARIANTS option in .conf file \n" }
    if(! $opt{SOMATIC_VARIANTS}){ die "ERROR: No SOMATIC_VARIANTS option in .conf file \n" }
    if(! $opt{COPY_NUMBER}){ die "ERROR: No COPY_NUMBER option in .conf file \n" }
    if(! $opt{ANNOTATE_VARIANTS}){ die "ERROR: No ANNOTATE_VARIANTS option in .conf file \n" }
    if(! $opt{CHECKING}){ die "ERROR: No CHECKING option in .conf file \n" }
    
    ### Module Settings / tools
    if(! $opt{GENOME}){ die "ERROR: No GENOME found in .ini file\n" }
    elsif(! -e $opt{GENOME}){ die"ERROR: $opt{GENOME} does not exist\n"}
    if(! $opt{SAMBAMBA_PATH}){ die "ERROR: No SAMBAMBA_PATH found in .ini file\n" }
    if(! $opt{QUEUE_PATH}){ die "ERROR: No PICARD_PATH found in .ini file\n" }
    ## PRESTATS
    if($opt{PRESTATS} eq "yes"){
	if(! $opt{FASTQC_PATH}){ die "ERROR: No FASTQC_PATH found in .ini file\n" }
	if(! $opt{PRESTATS_THREADS}){ die "ERROR: No PRESTATS_THREADS found in .ini file\n" }
	if(! $opt{PRESTATS_MEM}){ die "ERROR: No PRESTATS_MEM found in .ini file\n" }
	if(! $opt{PRESTATS_QUEUE}){ die "ERROR: No PRESTATS_QUEUE found in .ini file\n" }
	if(! $opt{PRESTATS_PROJECT}){ die "ERROR: No PRESTATS_PROJECT found in .ini file\n" }
    }
    ## MAPPING
    if($opt{MAPPING} eq "yes"){
	if(! $opt{BWA_PATH}){ die "ERROR: No BWA_PATH found in .ini file\n" }
	if(! $opt{MAPPING_THREADS}){ die "ERROR: No MAPPING_THREADS found in .ini file\n" }
	if(! $opt{MAPPING_MEM}){ die "ERROR: No MAPPING_MEM found in .ini file\n" }
	if(! $opt{MAPPING_QUEUE}){ die "ERROR: No MAPPING_QUEUE found in .ini file\n" }
	if(! $opt{MAPPING_PROJECT}){ die "ERROR: No MAPPING_PROJECT found in .ini file\n" }
	if(! $opt{MAPPING_MODE}){ die "ERROR: No MAPPING_MODE found in .ini file\n" }
	if(! $opt{MAPPING_MARKDUP}){ die "ERROR: No MAPPING_MARKDUP found in .ini file\n" }
    }
    ## POSTSTATS
    if($opt{POSTSTATS} eq "yes"){
	if(! $opt{PICARD_PATH}){ die "ERROR: No PICARD_PATH found in .ini file\n" }
	if(! $opt{POSTSTATS_THREADS}){ die "ERROR: No POSTSTATS_THREADS found in .ini file\n" }
	if(! $opt{POSTSTATS_MEM}){ die "ERROR: No POSTSTATS_MEM found in .ini file\n" }
	if(! $opt{POSTSTATS_QUEUE}){ die "ERROR: No POSTSTATS_THREADS found in .ini file\n" }
	if(! ($opt{POSTSTATS_TARGETS}) && ! ($opt{POSTSTATS_BAITS}) ){
	    if(! $opt{POSTSTATS_COVERAGECAP}){ die "ERROR: No POSTSTATS_COVERAGECAP found in .ini file\n" }
	}
	if( $opt{POSTSTATS_TARGETS} && ! -e $opt{POSTSTATS_TARGETS}){ die "ERROR: $opt{POSTSTATS_TARGETS} does not exist\n" }
	if( $opt{POSTSTATS_BAITS} && ! -e $opt{POSTSTATS_BAITS}){ die "ERROR: $opt{POSTSTATS_BAITS} does not exist\n" }
    }
    ## INDELREALIGNMENT
    if($opt{INDELREALIGNMENT} eq "yes"){
	if(! $opt{REALIGNMENT_MASTERQUEUE}){ die "ERROR: No REALIGNMENT_MASTERQUEUE found in .ini file\n" }
	if(! $opt{REALIGNMENT_MASTERTHREADS}){ die "ERROR: No REALIGNMENT_MASTERTHREADS found in .ini file\n" }
	if(! $opt{REALIGNMENT_QUEUE}){ die "ERROR: No REALIGNMENT_QUEUE found in .ini file\n" }
	if(! $opt{REALIGNMENT_PROJECT}){ die "ERROR: No REALIGNMENT_PROJECT found in .ini file\n" }
	if(! $opt{REALIGNMENT_THREADS}){ die "ERROR: No REALIGNMENT_THREADS found in .ini file\n" }
	if(! $opt{REALIGNMENT_MERGETHREADS}){ die "ERROR: No REALIGNMENT_MERGETHREADS found in .ini file\n" }
	if(! $opt{REALIGNMENT_MEM}){ die "ERROR: No REALIGNMENT_MEM found in .ini file\n" }
	if(! $opt{REALIGNMENT_SCALA}){ die "ERROR: No REALIGNMENT_SCALA found in .ini file\n" }
	if(! $opt{REALIGNMENT_SCATTER}){ die "ERROR: No REALIGNMENT_SCATTER found in .ini file\n" }
	if(! $opt{REALIGNMENT_MODE}){ die "ERROR: No REALIGNMENT_MODE found in .ini file\n" }
	if(! $opt{QUEUE_RETRY}){ die "ERROR: No QUEUE_RETRY found in .ini file\n" }
    }
    ## BASEQUALITYRECAL
    if($opt{BASEQUALITYRECAL} eq "yes"){
	if(! $opt{BASERECALIBRATION_MASTERQUEUE}){ die "ERROR: No BASERECALIBRATION_QUEUE found in .ini file\n" }
	if(! $opt{BASERECALIBRATION_MASTERTHREADS}){ die "ERROR: No BASERECALIBRATION_THREADS found in .ini file\n" }
	if(! $opt{BASERECALIBRATION_QUEUE}){ die "ERROR: No BASERECALIBRATION_QUEUE found in .ini file\n" }
	if(! $opt{BASERECALIBRATION_THREADS}){ die "ERROR: No BASERECALIBRATION_THREADS found in .ini file\n" }
	if(! $opt{BASERECALIBRATION_MEM}){ die "ERROR: No BASERECALIBRATION_MEM found in .ini file\n" }
	if(! $opt{BASERECALIBRATION_SCALA}){ die "ERROR: No BASERECALIBRATION_SCALA found in .ini file\n" }
	if(! $opt{BASERECALIBRATION_SCATTER}){ die "ERROR: No BASERECALIBRATION_SCATTER found in .ini file\n" }
	if(! $opt{QUEUE_RETRY}){ die "ERROR: No QUEUE_RETRY found in .ini file\n" }
    }
    ## VARIANT_CALLING
    if($opt{VARIANT_CALLING} eq "yes"){
	if(! $opt{CALLING_MASTERQUEUE}){ die "ERROR: No CALLING_MASTERQUEUE found in .ini file\n" }
	if(! $opt{CALLING_MASTERTHREADS}){ die "ERROR: No CALLING_MASTERTHREADS found in .ini file\n" }
	if(! $opt{CALLING_QUEUE}){ die "ERROR: No CALLING_QUEUE found in .ini file\n" }
	if(! $opt{CALLING_THREADS}){ die "ERROR: No CALLING_THREADS found in .ini file\n" }
	if(! $opt{CALLING_MEM}){ die "ERROR: No CALLING_QUEUE found in .ini file\n" }
	if(! $opt{CALLING_SCATTER}){ die "ERROR: No CALLING_SCATTER found in .ini file\n" }
	if(! $opt{CALLING_GVCF}){ die "ERROR: No CALLING_GVCF found in .ini file\n"}
	if(! $opt{CALLING_SCALA}){ die "ERROR: No CALLING_SCALA found in .ini file\n" }
	if($opt{CALLING_UGMODE}){ 
	    if($opt{CALLING_UGMODE} ne "SNP" and $opt{CALLING_UGMODE} ne "INDEL" and $opt{CALLING_UGMODE} ne "BOTH"){ die "ERROR: UGMODE: $opt{CALLING_UGMODE} does not exist use SNP, INDEL or BOTH\n"}
	}
	if(! $opt{CALLING_STANDCALLCONF}){ die "ERROR: No CALLING_STANDCALLCONF found in .ini file\n" }
	if(! $opt{CALLING_STANDEMITCONF}){ die "ERROR: No CALLING_STANDEMITCONF found in .ini file\n" }
	if( $opt{CALLING_TARGETS} && ! -e $opt{CALLING_TARGETS}) { die"ERROR: $opt{CALLING_TARGETS} does not exist\n" }
	if( $opt{CALLING_DBSNP} && ! -e $opt{CALLING_DBSNP}) { die"ERROR: $opt{CALLING_DBSNP} does not exist\n" }
	if(! $opt{QUEUE_RETRY}){ die "ERROR: No QUEUE_RETRY found in .ini file\n" }
    }
    ## FILTER_VARIANTS
    if($opt{FILTER_VARIANTS} eq "yes"){
	if(! $opt{FILTER_MASTERQUEUE}){ die "ERROR: No FILTER_MASTERQUEUE found in .ini file\n" }
	if(! $opt{FILTER_MASTERTHREADS}){ die "ERROR: No FILTER_MASTERTHREADS found in .ini file\n" }    
	if(! $opt{FILTER_QUEUE}){ die "ERROR: No FILTER_QUEUE found in .ini file\n" }
	if(! $opt{FILTER_THREADS}){ die "ERROR: No FILTER_THREADS found in .ini file\n" }
	if(! $opt{FILTER_MEM}){ die "ERROR: No FILTER_QUEUE found in .ini file\n" }
	if(! $opt{FILTER_SCATTER}){ die "ERROR: No FILTER_SCATTER found in .ini file\n" }
	if(! $opt{FILTER_SCALA}){ die "ERROR: No FILTER_SCALA found in .ini file\n" }
	if(! $opt{FILTER_MODE}){ die "ERROR: No FILTER_MODE  found in .ini file\n" }
	if($opt{FILTER_MODE} ne "SNP" and $opt{FILTER_MODE} ne "INDEL" and $opt{FILTER_MODE} ne "BOTH"){ die "ERROR: FILTER_MODE $opt{FILTER_MODE} does not exist use SNP, INDEL or BOTH\n"}
	if ($opt{FILTER_MODE} eq "SNP" || $opt{FILTER_MODE} eq "BOTH") {
	    if(! $opt{FILTER_SNPNAME}){ die "ERROR: No FILTER_SNPNAME found in .ini file\n" }
	    if(! $opt{FILTER_SNPEXPR}){ die "ERROR: No FILTER_SNPEXPR  found in .ini file\n" }
	}
	if ($opt{FILTER_MODE} eq "INDEL" || $opt{FILTER_MODE} eq "BOTH") {
	    if(! $opt{FILTER_INDELNAME}){ die "ERROR: No FILTER_INDELNAME found in .ini file\n" }
	    if(! $opt{FILTER_INDELEXPR}){ die "ERROR: No FILTER_INDELEXPR found in .ini file\n" }
	}
	if(! $opt{QUEUE_RETRY}){ die "ERROR: No QUEUE_RETRY found in .ini file\n" }
    }
    ## SOMATIC_VARIANTS
    if($opt{SOMATIC_VARIANTS} eq "yes"){
	if(! $opt{SAMTOOLS_PATH}){ die "ERROR: NO SAMTOOLS_PATH found in .ini file\n" }
	if(! $opt{SOMVAR_TARGETS}){ die "ERROR: NO SOMVAR_TARGETS found in .ini file\n" }
	if(! $opt{SOMVAR_STRELKA}){ die "ERROR: NO SOMVAR_STRELKA found in .ini file\n" }
	if($opt{SOMVAR_STRELKA} eq "yes"){
	    if(! $opt{STRELKA_PATH}){ die "ERROR: NO STRELKA_PATH found in .ini file\n" }
	    if(! $opt{STRELKA_INI}){ die "ERROR: NO STRELKA_INI found in .ini file\n" }
	    if(! $opt{STRELKA_QUEUE}){ die "ERROR: NO STRELKA_QUEUE found in .ini file\n" }
	    if(! $opt{STRELKA_THREADS}){ die "ERROR: NO STRELKA_THREADS found in .ini file\n" }
	}
	if(! $opt{SOMVAR_VARSCAN}){ die "ERROR: NO SOMVAR_VARSCAN found in .ini file\n" }
	if($opt{SOMVAR_VARSCAN} eq "yes"){
	    if(! $opt{VARSCAN_PATH}){ die "ERROR: NO VARSCAN_PATH found in .ini file\n" }
	    if(! $opt{VARSCAN_QUEUE}){ die "ERROR: NO VARSCAN_QUEUE found in .ini file\n" }
	    if(! $opt{VARSCAN_THREADS}){ die "ERROR: NO VARSCAN_THREADS found in .ini file\n" }
	    if(! $opt{VARSCAN_SETTINGS}){ die "ERROR: NO VARSCAN_SETTINGS found in .ini file\n" }
	    if(! $opt{VARSCAN_POSTSETTINGS}){ die "ERROR: NO VARSCAN_POSTSETTINGS found in .ini file\n" }
	}
	if(! $opt{SOMVAR_FREEBAYES}){ die "ERROR: NO SOMVAR_FREEBAYES found in .ini file\n" }
	if($opt{SOMVAR_FREEBAYES} eq "yes"){
	    if(! $opt{FREEBAYES_PATH}){ die "ERROR: NO FREEBAYES_PATH found in .ini file\n" }
	    if(! $opt{VCFSAMPLEDIFF_PATH}){ die "ERROR: NO VCFSAMPLEDIFF_PATH found in .ini file\n" }
	    if(! $opt{BIOVCF_PATH}){ die "ERROR: NO BIOVCF_PATH found in .ini file\n" }
	    if(! $opt{FREEBAYES_QUEUE}){ die "ERROR: NO FREEBAYES_QUEUE found in .ini file\n" }
	    if(! $opt{FREEBAYES_THREADS}){ die "ERROR: NO FREEBAYES_THREADS found in .ini file\n" }
	    if(! $opt{FREEBAYES_SETTINGS}){ die "ERROR: NO FREEBAYES_SETTINGS found in .ini file\n" }
	    if(! $opt{FREEBAYES_SOMATICFILTER}){ die "ERROR: NO FREEBAYES_SOMATICFILTER found in .ini file\n" }
	    if(! $opt{FREEBAYES_GERMLINEFILTER}){ die "ERROR: NO FREEBAYES_GERMLINEFILTER found in .ini file\n" }
	}
	if(! $opt{SOMVARMERGE_QUEUE}){ die "ERROR: NO SOMVARMERGE_QUEUE found in .ini file\n" }
	if(! $opt{SOMVARMERGE_THREADS}){ die "ERROR: NO SOMVARMERGE_THREADS found in .ini file\n" }
    }
    ## COPY_NUMBER
    if($opt{COPY_NUMBER} eq "yes"){
	if(! $opt{CNVCHECK_QUEUE} ) { die "ERROR: $opt{CNVCHECK_QUEUE} does not exist\n"}
	if(! $opt{CNVCHECK_THREADS} ) { die "ERROR: $opt{CNVCHECK_THREADS} does not exist\n"}
	if(! $opt{CNV_CONTRA}){ die "ERROR: NO CNV_CONTRA found in .ini file\n" }
	if($opt{CNV_CONTRA} eq "yes"){
	    if(! $opt{CONTRA_PATH}){ die "ERROR: NO CONTRA_PATH found in .ini file\n" }
	    if(! $opt{CONTRA_QUEUE}){ die "ERROR: NO CONTRA_QUEUE found in .ini file\n" }
	    if(! $opt{CONTRA_THREADS}){ die "ERROR: NO CONTRA_THREADS found in .ini file\n" }
	    if(! $opt{CONTRA_TARGETS}){ die "ERROR: NO CONTRA_TARGETS found in .ini file\n" }
	    if(! $opt{CONTRA_FLAGS}){ die "ERROR: NO CONTRA_FLAGS found in .ini file\n" }
	}
    }
    ## ANNOTATE_VARIANTS
    if($opt{ANNOTATE_VARIANTS} eq "yes"){
	if(! $opt{SNPEFF_PATH}){ die "ERROR: No SNPEFF_PATH found in .ini file\n" }
	if(! $opt{IGVTOOLS_PATH}){ die "ERROR: No IGVTOOLS_PATH found in .ini file\n" }
	if(! $opt{ANNOTATE_QUEUE}){ die "ERROR: No ANNOTATE_QUEUE found in .ini file\n" }
	if(! $opt{ANNOTATE_THREADS}){ die "ERROR: No ANNOTATE_THREADS found in .ini file\n" }
	if(! $opt{ANNOTATE_MEM}){ die "ERROR: No ANNOTATE_MEM found in .ini file\n" }
	if(! $opt{ANNOTATE_SNPEFF}){ die "ERROR: No ANNOTATE_SNPEFF found in .ini file\n" }
	if($opt{ANNOTATE_SNPEFF} eq "yes"){
	    if(! $opt{ANNOTATE_DB}){ die "ERROR: No ANNOTATE_DB found in .ini file\n" }
	    if(! $opt{ANNOTATE_FLAGS}){ die "ERROR: No ANNOTATE_FLAGS found in .ini file\n" }
	}
	if(! $opt{ANNOTATE_SNPSIFT}){ die "ERROR: No ANNOTATE_SNPSIFT found in .ini file\n" }
	if($opt{ANNOTATE_SNPSIFT} eq "yes"){
	    if(! $opt{ANNOTATE_DBNSFP}){ die "ERROR: No ANNOTATE_DBNSFP found in .ini file\n" }
	    elsif( $opt{ANNOTATE_DBNSFP} && ! -e $opt{ANNOTATE_DBNSFP}) { die"ERROR: $opt{ANNOTATE_DBNSFP} does not exist\n" }
	    if(! $opt{ANNOTATE_FIELDS}){ die "ERROR: No ANNOTATE_FIELDS found in .ini file\n" }
	}
	if(! $opt{ANNOTATE_FREQUENCIES}){ die "ERROR: No ANNOTATE_FREQUENCIES found in .ini file\n" }
	if($opt{ANNOTATE_FREQUENCIES} eq "yes"){
	    if(! $opt{ANNOTATE_FREQNAME}){ die "ERROR: No ANNOTATE_FREQNAME found in .ini file\n" }
	    if(! $opt{ANNOTATE_FREQDB}){ die "ERROR: No ANNOTATE_FREQDB found in .ini file\n" }
	    elsif( $opt{ANNOTATE_FREQDB} && ! -e $opt{ANNOTATE_FREQDB}) { die"ERROR: $opt{ANNOTATE_FREQDB} does not exist\n" }
	    if(! $opt{ANNOTATE_FREQINFO}){ die "ERROR: No ANNOTATE_FREQINFO found in .ini file\n" }
	}
	if(! $opt{ANNOTATE_IDFIELD}){ die "ERROR: No ANNOTATE_IDFIELD found in .ini file\n" }
	if($opt{ANNOTATE_IDFIELD} eq "yes"){
	    if(! $opt{ANNOTATE_IDNAME}){ die "ERROR: No ANNOTATE_IDNAME found in .ini file\n" }
	    if(! $opt{ANNOTATE_IDDB}){ die "ERROR: No ANNOTATE_IDDB found in .ini file\n" }
	}
    }
}
