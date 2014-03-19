#!/usr/bin/perl -w

##################################################################################################################################################
###This script is designed check the final output of the illumina pipeline
###
###
###Author: R.F.Ernst
###Latest change: skeleton
###TODO:
##################################################################################################################################################

package illumina_check;

use strict;
use POSIX qw(tmpnam);

sub runCheck {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};

}

sub readConfiguration{
    my $configuration = shift;

#    my %opt = (
#	'CHECKING_QUEUE'		=> undef,
#	'CHECKING_THREADS'	=> undef,
#	'CHECKING_MEM'		=> undef,
#	'OUTPUT_DIR'		=> undef,
#	'RUNNING_JOBS'		=> {} #do not use in .conf file
#    );

#    foreach my $key (keys %{$configuration}){
#	$opt{$key} = $configuration->{$key};
#    }

    if(! $opt{CHECKING_QUEUE}){ die "ERROR: No CALLING_QUEUE found in .conf file\n" }
    if(! $opt{CHECKING_THREADS}){ die "ERROR: No CALLING_THREADS found in .conf file\n" }
    if(! $opt{CHECKING_MEM}){ die "ERROR: No CALLING_QUEUE found in .conf file\n" }
    if(! $opt{OUTPUT_DIR}){ die "ERROR: No OUTPUT_DIR found in .conf file\n" }

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