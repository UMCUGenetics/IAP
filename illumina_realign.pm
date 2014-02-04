#!/usr/bin/perl -w

##################################################################################################################################################
###This script is designed to run GATK indelrealignment using GATK Queue
###
###
###Author: S.W.Boymans
###Latest change: Created skeleton
###
###TODO: A lot
##################################################################################################################################################


package illumina_poststats;

use strict;
use POSIX qw(tmpnam);


sub runPostStats {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};


    print "Running post mapping statistics for the following BAM-files:\n";
    
    foreach my $sample (keys %{$opt{RUNNING_JOBS}}){
	print "\t$opt{OUTPUT_DIR}/mapping/$sample\_dedup.bam\n";
	print "\tDEBUG: waiting for job(s): ".join(",",@{$opt{RUNNING_JOBS}->{$sample}}). " to finish first\n";
	
	
	####REST OF CODE HERE
	
	
	
	
	
	
    }

}

sub readConfiguration{
    my $configuration = shift;
    
    my %opt = (
	
	'QUALIMAP_PATH'		=> undef,
	'SAMBAMBA_PATH'		=> undef,
	'CLUSTER_PATH'  	=> undef,
	'POSTSTATS_THREADS'	=> undef,
	'POSTSTATS_MEM'		=> undef,
	'POSTSTATS_QUEUE'	=> undef,
	'POSTSTATS_PROJECT'	=> undef,
	'CLUSTER_TMP'		=> undef,
	'GENOME'		=> undef,
	'OUTPUT_DIR'		=> undef,
	'RUNNING_JOBS'		=> {} #do not use in .conf file
    );

    foreach my $key (keys %{$configuration}){
	$opt{$key} = $configuration->{$key};
    }

    if(! $opt{QUALIMAP_PATH}){ die "ERROR: No QUALIMAP_PATH found in .counf file\n" }
    if(! $opt{SAMBAMBA_PATH}){ die "ERROR: No SAMBAMBA_PATH found in .conf file\n" }
    if(! $opt{POSTSTATS_PROJECT}){ die "ERROR: No POSTSTATS_PROJECT found in .ini file\n" }
    if(! $opt{POSTSTATS_THREADS}){ die "ERROR: No POSTSTATS_THREADS found in .ini file\n" }
    if(! $opt{POSTSTATS_MEM}){ die "ERROR: No POSTSTATS_MEM found in .ini file\n" }
    if(! $opt{CLUSTER_PATH}){ die "ERROR: No CLUSTER_PATH found in .conf file\n" }
    if(! $opt{CLUSTER_TMP}){ die "ERROR: No CLUSTER_TMP found in .conf file\n" }
    if(! $opt{GENOME}){ die "ERROR: No GENOME found in .conf file\n" }
    if(! $opt{OUTPUT_DIR}){ die "ERROR: No OUTPUT_DIR found in .conf file\n" }
    if(! $opt{RUNNING_JOBS}){ die "ERROR: No RUNNING_JOBS found in .conf file\n" }

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