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
    
    my $mainJobID = "$opt{OUTPUT_DIR}/jobs/PreStatsMainJob_".get_job_id().".sh";

    open (QSUB,">$mainJobID") or die "ERROR: Couldn't create $mainJobID\n";
    print QSUB "\#!/bin/sh\n\n. $opt{CLUSTER_PATH}/settings.sh\n\n";
    print "Creating FASTQC report for the following fastq.gz files:\n";

    foreach my $input (keys %{$opt{FASTQ}}){
	my $coreName = undef;
	$coreName = (split("/", $input))[-1];
	$coreName =~ s/\.fastq.gz//;
	my ($sampleName) =  split("_", $coreName);
	print "\t$input\n"; #print fastq filename

	if(! -e "$opt{OUTPUT_DIR}/$sampleName/logs/PreStats_$sampleName.done"){

	    my $preStatsJobId = "PreStat_$coreName\_".get_job_id();
	    push(@{$jobIds->{$sampleName}}, $preStatsJobId);
	    open PS,">$opt{OUTPUT_DIR}/$sampleName/jobs/$preStatsJobId.sh";
	    print PS "\#!/bin/sh\n\n";
	    print PS "cd $opt{OUTPUT_DIR}/$sampleName\n\n";
	    print PS "echo \"Start PreStats\t\" `date` \"\t$coreName\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	    print PS "$opt{FASTQC_PATH}/fastqc $input -o QCStats --noextract\n";
	    print PS "touch logs/PreStats_$sampleName.done\n";
	    print PS "echo \"End PreStats\t\" `date` \"\t$coreName\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	    close PS;

	    print QSUB "qsub -pe threaded $opt{PRESTATS_THREADS} -m abe -M $opt{MAIL} -q $opt{PRESTATS_QUEUE} -P $opt{PRESTATS_PROJECT} -o $opt{OUTPUT_DIR}/$sampleName/logs/PreStat_$coreName.out -e $opt{OUTPUT_DIR}/$sampleName/logs/PreStats_$coreName.err -N $preStatsJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$preStatsJobId.sh\n";
	} else {
	    warn "\t WARNING: FASTQC report for $input already exists, skipping.\n";
	}

    }

    close QSUB;

    system("sh $mainJobID");
}

sub readConfiguration {
    my $configuration = shift;
    my %opt;

    foreach my $key (keys %{$configuration}){
	    $opt{$key} = $configuration->{$key}; 
    }

    if(! $opt{FASTQC_PATH}){ die "ERROR: No FASTQC_PATH found in .conf file\n" }
    if(! $opt{PRESTATS_THREADS}){ die "ERROR: No PRESTATS_THREADS found in .conf file\n" }
    if(! $opt{PRESTATS_MEM}){ die "ERROR: No PRESTATS_MEM found in .conf file\n" }
    if(! $opt{PRESTATS_QUEUE}){ die "ERROR: No PRESTATS_QUEUE found in .conf file\n" }
    if(! $opt{PRESTATS_PROJECT}){ die "ERROR: No PRESTATS_PROJECT found in .conf file\n" }
    if(! $opt{CLUSTER_PATH}){ die "ERROR: No CLUSTER_PATH found in .conf file\n" }
    if(! $opt{OUTPUT_DIR}){ die "ERROR: No OUTPUT_DIR found in .conf file\n" }
    if(! $opt{FASTQ}){ die "ERROR: No FASTQ files specified\n" }
    if(! $opt{MAIL}){ die "ERROR: No MAIL adress specified\n" }
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