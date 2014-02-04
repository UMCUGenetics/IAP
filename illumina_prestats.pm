#!/usr/bin/perl -w


##################################################################################################################################################
###This script is designed to run FastQC on a set of fastq.gz files. Currently this script can only be run using the illumina_pipeline.pl script.
###In a later version I also plan to add per-sample summary PDFs, this is currrently not supported yet.
###
###Author: S.W.Boymans
###Latest change: Created first version
###
###TODO: Add system call actually firing jobs off
##################################################################################################################################################



package illumina_prestats;

use strict;
use POSIX qw(tmpnam);



sub runPreStats {
    
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
    my $jobIds = {};
    
    my $mainJobID = "$opt{OUTPUT_DIR}/".get_job_id()."_prestats_qsub.sh";

    open (QSUB,">$mainJobID") or die "ERROR: Couldn't create $mainJobID\n";
    print QSUB "\#!/bin/sh\n\n. $opt{CLUSTER_PATH}/settings.sh\n\n";
    print "Creating FASTQC report for the following fastq.gz files:\n";

    foreach my $input (keys %{$opt{FASTQ}}){
	my $coreName = undef;
	$coreName = (split("/", $input))[-1];
	$coreName =~ s/\.fastq.gz//;
	my ($sampleName) =  split("_", $coreName);
	print "\t$input\n";
    

	if(! -e "$opt{OUTPUT_DIR}/$sampleName"){
    	    mkdir("$opt{OUTPUT_DIR}/$sampleName") or die "ERROR: Couldn't create directory: $opt{OUTPUT_DIR}/$sampleName\n";
	}
	if(! -e "$opt{OUTPUT_DIR}/$sampleName/fastqc"){
    	    mkdir("$opt{OUTPUT_DIR}/$sampleName/fastqc") or die "ERROR: Couldn't create directory: $opt{OUTPUT_DIR}/$sampleName/fastqc\n";
	}
	if(! -e "$opt{OUTPUT_DIR}/$sampleName/jobs"){
    	    mkdir("$opt{OUTPUT_DIR}/$sampleName/jobs") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sampleName/jobs\n";
	}
	if(! -e "$opt{OUTPUT_DIR}/$sampleName/logs"){
    	    mkdir("$opt{OUTPUT_DIR}/$sampleName/logs") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sampleName/logs\n";
	}	
    
	if(! -e "$opt{OUTPUT_DIR}/$sampleName/fastqc/$input.done"){

	    my $preStatsJobId = "PRESTATS_$coreName\_".get_job_id();
	    push(@{$jobIds->{$sampleName}}, $preStatsJobId);
	    open PS,">$opt{OUTPUT_DIR}/$sampleName/jobs/$preStatsJobId.sh";
	    print PS "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
	    print PS "cd $opt{OUTPUT_DIR}/$sampleName\n\n";
	    print PS "uname -n > logs/$preStatsJobId.host\n";
	    print PS "echo \"FastQC\t\" `date` >> logs/$preStatsJobId.host\n";
	    print PS "$opt{FASTQC_PATH}/fastqc $input -o fastqc\n";
	    print PS "touch fastqc/$input.done\n";
	    close PS;
	    
	    print QSUB "qsub -pe threaded $opt{PRESTATS_THREADS} -q $opt{PRESTATS_QUEUE} -P $opt{PRESTATS_PROJECT} -o $opt{OUTPUT_DIR}/$sampleName/logs -e $opt{OUTPUT_DIR}/$sampleName/logs -N $preStatsJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$preStatsJobId.sh\n";
	}

    }

    close QSUB;

    #system();
    
}


sub readConfiguration {
    my $configuration = shift;
    my %opt = (
        'FASTQC_PATH'      	=> undef,
        'CLUSTER_PATH'  	=> undef,
	'CLUSTER_TMP'		=> undef,
	'PRESTATS_THREADS'	=> undef,
        'PRESTATS_MEM'		=> undef,
        'PRESTATS_QUEUE'	=> undef,
	'PRESTATS_PROJECT'	=> undef,
        'OUTPUT_DIR'		=> undef,
	'FASTQ'			=> []
    );

    foreach my $key (keys %{$configuration}){
	    $opt{$key} = $configuration->{$key}; 
    }

    if(! $opt{FASTQC_PATH}){ die "ERROR: No BWA_PATH found in .conf file\n" }
    
    if(! $opt{PRESTATS_THREADS}){ die "ERROR: No PRESTATS_THREADS found in .ini file\n" }
    if(! $opt{PRESTATS_MEM}){ die "ERROR: No PRESTATS_MEM found in .ini file\n" }
    if(! $opt{PRESTATS_QUEUE}){ die "ERROR: No PRESTATS_QUEUE found in .ini file\n" }
    if(! $opt{PRESTATS_PROJECT}){ die "ERROR: No PRESTATS_PROJECT found in .ini file\n" }
    if(! $opt{CLUSTER_PATH}){ die "ERROR: No CLUSTER_PATH found in .conf file\n" }
    if(! $opt{CLUSTER_TMP}){ die "ERROR: No CLUSTER_TMP found in .conf file\n" }
    if(! $opt{OUTPUT_DIR}){ die "ERROR: No OUTPUT_DIR found in .conf file\n" }
    if(! $opt{FASTQ}){ die "ERROR: No FASTQ files specified\n" }


    return \%opt;
}


############
sub get_job_id {
   my $id = tmpnam(); 
      $id=~s/\/tmp\/file//;
   return $id;
}
############ 
1;
