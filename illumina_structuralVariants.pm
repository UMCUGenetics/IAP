#!/usr/bin/perl -w

##################################################################
### illumina_structuralVariants.pm
### - Run structural variant caller Delly
###
### Author: R.F.Ernst & M. van Roosmalen
##################################################################

package illumina_structuralVariants;

use strict;
use POSIX qw(tmpnam);

sub runDelly {
    ###
    # Run structural variant caller Delly
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my @sampleBams;
    my @runningJobs;
    my $jobID = "SV_".get_job_id();

    ### Skip sv calling if .done file exists
    if (-e "$opt{OUTPUT_DIR}/logs/StructuralVariants.done"){
	print "WARNING: $opt{OUTPUT_DIR}/logs/StructuralVariants.done exists, skipping \n";
	return \%opt;
    }

    ### Application logic goes here.
    ### Important variables available:
    # $opt{GENOME}
    # @{$opt{SAMPLES}}
    # $runName -> can/should be used for naming vcf file
    
    # get bam files and running jobs per sample
    # foreach my $sample (@{$opt{SAMPLES}}){
    # 	my $sampleBam = "$opt{OUTPUT_DIR}/$sample/mapping/$opt{BAM_FILES}->{$sample}";
    ## Running jobs
    #	if ( @{$opt{RUNNING_JOBS}->{$sample}} ){
    #		push( @runningJobs, @{$opt{RUNNING_JOBS}->{$sample}} );
    #	}
    # }

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