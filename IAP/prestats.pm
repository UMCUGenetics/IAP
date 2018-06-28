#!/usr/bin/perl -w

###################################
### illumina_prestats.pm
### - Run fastqc for each fastq
###
### Author: S.W.Boymans & H.H.D.Kerstens
###################################

package IAP::prestats;

use strict;
use File::Temp;
use File::Basename;
use IAP::sge;

sub runPreStats {
    ###
    # Run fastqc
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my $jobIds = {};

    print "Creating FASTQC report for the following fastq.gz files:\n";

    foreach my $input (keys %{$opt{FASTQ}}){
	my $coreName = undef;
	$coreName = (split("/", $input))[-1];
	$coreName =~ s/\.fastq.gz//;
	my ($sampleName) =  split("_", $coreName);
	print "\t$input\n"; #print fastq filename

	if(! -e "$opt{OUTPUT_DIR}/$sampleName/logs/PreStats_$coreName.done"){

	    my $preStatsJobId = "PreStat_$coreName\_".get_job_id();
	    open PS,">$opt{OUTPUT_DIR}/$sampleName/jobs/$preStatsJobId.sh";
	    print PS "\#!/bin/sh\n\n";
	    print PS "cd $opt{OUTPUT_DIR}/$sampleName\n\n";
	    print PS "echo \"Start PreStats\t\" `date` \"\t$coreName\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	    print PS "$opt{FASTQC_PATH}/fastqc -o QCStats --noextract -t $opt{PRESTATS_THREADS} $input\n";
	    print PS "if [ -s $opt{OUTPUT_DIR}/$sampleName/QCStats/$coreName\_fastqc.zip -a -s $opt{OUTPUT_DIR}/$sampleName/QCStats/$coreName\_fastqc.html ]\n";
	    print PS "then\n";
	    print PS "\ttouch logs/PreStats_$coreName.done\n";
	    print PS "fi\n";
	    print PS "echo \"End PreStats\t\" `date` \"\t$coreName\t\" `uname -n` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sampleName.log\n";
	    close PS;

	    my $qsub = &qsubTemplate(\%opt,"PRESTATS");
	    system $qsub." -o ".$opt{OUTPUT_DIR}."/".$sampleName."/logs/PreStat_".$coreName.".out -e ".$opt{OUTPUT_DIR}."/".$sampleName."/logs/PreStats_".$coreName.".err -N ".$preStatsJobId." ".$opt{OUTPUT_DIR}."/".$sampleName."/jobs/".$preStatsJobId.".sh";
	    push(@{$opt{RUNNING_JOBS}->{"preStats"}}, $preStatsJobId);
	} else {
	    print "\t WARNING: FASTQC report for $input already exists, skipping.\n";
	}

    }
    return \%opt;
}

############
sub get_job_id {
   my $id = tmpnam(); 
   $id=basename($id);
   return $id;
}
############
1;
