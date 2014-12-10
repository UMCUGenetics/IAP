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
if(! $opt{INIFILE}){ die "ERROR: No INIFILE found in .conf file\n" }
if(! $opt{OUTPUT_DIR}){ die "ERROR: No OUTPUT_DIR found in .conf file\n" }
if(! ($opt{FASTQ} || $opt{BAM} || $opt{VCF}) ){ die "ERROR: No FASTQ/BAM/VCF files specified\n" }
if(! $opt{PRESTATS}){ die "ERROR: No PRESTATS option in .conf file \n" }
if(! $opt{MAPPING}){ die "ERROR: No MAPPING option in .conf file \n" }
if(! $opt{POSTSTATS}){ die "ERROR: No POSTSTATS option in .conf file \n" }
if(! $opt{INDELREALIGNMENT}){ die "ERROR: No INDELREALIGNMENT option in .conf file \n" }
if(! $opt{BASEQUALITYRECAL}){ die "ERROR: No BASEQUALITYRECAL option in .conf file \n" }
if(! $opt{VARIANT_CALLING}){ die "ERROR: No VARIANT_CALLING option in .conf file \n" }
if(! $opt{FILTER_VARIANTS}){ die "ERROR: No FILTER_VARIANTS option in .conf file \n" }
if(! $opt{SOMATIC_VARIANTS}){ die "ERROR: No SOMATIC_VARIANTS option in .conf file \n" }
if(! $opt{ANNOTATE_VARIANTS}){ die "ERROR: No ANNOTATE_VARIANTS option in .conf file \n" }
if(! $opt{CHECKING}){ die "ERROR: No CHECKING option in .conf file \n" }

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
    
### Variant Caller or vcf input
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

### Somatic variant caller
if($opt{SOMATIC_VARIANTS} eq "yes"){
    print "\n###SCHEDULING SOMATIC VARIANT CALLERS####\n";
    $opt_ref = illumina_somaticVariants::parseSamples(\%opt);
    %opt = %$opt_ref;
    
    illumina_somaticVariants::runSomaticVariantCallers(\%opt);

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
	#if(! -e "$opt{OUTPUT_DIR}/$sample/variants"){
	#    mkdir("$opt{OUTPUT_DIR}/$sample/variants") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sample/variants\n";
	#}
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
