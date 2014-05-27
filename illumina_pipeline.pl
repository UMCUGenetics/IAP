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
use illumina_annotateVariants;
use illumina_check;

### Check correct usage
die usage() if @ARGV == 0;

### initiate opt hash with settings
my %opt;
my $configurationFile;

%opt = (
    'RUNNING_JOBS'		=> {}, #do not use in .conf or .ini
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
    } elsif($key eq 'FASTQ') {
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
if(! $opt{FASTQ}){ die "ERROR: No FASTQ files specified\n" }
if(! $opt{PRESTATS}){ die "ERROR: No PRESTATS option in .conf file \n" }
if(! $opt{MAPPING}){ die "ERROR: No MAPPING option in .conf file \n" }
if(! $opt{POSTSTATS}){ die "ERROR: No POSTSTATS option in .conf file \n" }
if(! $opt{INDELREALIGNMENT}){ die "ERROR: No INDELREALIGNMENT option in .conf file \n" }
if(! $opt{BASEQUALITYRECAL}){ die "ERROR: No BASEQUALITYRECAL option in .conf file \n" }
if(! $opt{VARIANT_CALLING}){ die "ERROR: No VARIANT_CALLING option in .conf file \n" }
if(! $opt{FILTER_VARIANTS}){ die "ERROR: No FILTER_VARIANTS option in .conf file \n" }
if(! $opt{ANNOTATE_VARIANTS}){ die "ERROR: No ANNOTATE_VARIANTS option in .conf file \n" }
if(! $opt{CHECKING}){ die "ERROR: No CHECKING option in .conf file \n" }

###Read samples from FASTQ's
getSamples();
createOutputDirs();

### Copy ini file to logs dir
system "cp $opt{INIFILE} $opt{OUTPUT_DIR}/logs";

### Start pipeline components
if($opt{PRESTATS} eq "yes"){
    print "###SCHEDULING PRESTATS###\n";
    illumina_prestats::runPreStats(\%opt);
}

if($opt{MAPPING} eq "yes"){
    print "\n###SCHEDULING MAPPING###\n";
    my $mappingJobs = illumina_mapping::runMapping(\%opt);
    
    foreach my $sample (keys %{$mappingJobs}){
	push (@{$opt{RUNNING_JOBS}->{$sample}} , $mappingJobs->{$sample});
    }

}

if($opt{POSTSTATS} eq "yes"){
    print "\n###SCHEDULING POSTSTATS###\n";
    my $postStatsJob = illumina_poststats::runPostStats(\%opt);
    push (@{$opt{RUNNING_JOBS}->{'postStats'}} , $postStatsJob);
}

if($opt{INDELREALIGNMENT} eq "yes"){
    print "\n###SCHEDULING INDELREALIGNMENT###\n";
    my $realignJobs = illumina_realign::runRealignment(\%opt);
    
    foreach my $sample (keys %{$realignJobs}){
	push (@{$opt{RUNNING_JOBS}->{$sample}} , $realignJobs->{$sample});
    }
}

if($opt{BASEQUALITYRECAL} eq "yes"){
    print "\n###SCHEDULING BASERECALIBRATION###\n";
    my %baseRecalJobs = illumina_baseRecal::runBaseRecalibration(\%opt);
    
    foreach my $sample (keys %baseRecalJobs){
	push (@{$opt{RUNNING_JOBS}->{$sample}} , $baseRecalJobs{$sample});
    }
    
}

if($opt{VARIANT_CALLING} eq "yes"){
    print "\n###SCHEDULING VARIANT CALLING####\n";
    my $VCJob = illumina_calling::runVariantCalling(\%opt);
    
    foreach my $sample (@{$opt{SAMPLES}}){
	push (@{$opt{RUNNING_JOBS}->{$sample}} , $VCJob);
    }
}

if($opt{FILTER_VARIANTS} eq "yes"){
    print "\n###SCHEDULING VARIANT FILTRATION####\n";
    my $FVJob = illumina_filterVariants::runFilterVariants(\%opt);
    
    foreach my $sample (@{$opt{SAMPLES}}){
	push (@{$opt{RUNNING_JOBS}->{$sample}} , $FVJob);
    }
}

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
    foreach my $input (keys %{$opt{FASTQ}}){
	my $fastqFile = (split("/", $input))[-1];
	my $sampleName =  (split("_", $fastqFile))[0];
	$samples{$sampleName} ++;
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


=cut
############
sub get_job_id {
   my $id = tmpnam(); 
      $id=~s/\/tmp\/file//;
   return $id;
}
