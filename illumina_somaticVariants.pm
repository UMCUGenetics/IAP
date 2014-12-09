#!/usr/bin/perl -w

##################################################################################################################################################
###This script is designed to run Somatic Variant callers and combine the result.
###
###
###Author: R.F.Ernst
###Latest change:
###TODO: CPCT sample name parsing, add callers, combine the result
##################################################################################################################################################

package illumina_somaticVariants;

use strict;
use POSIX qw(tmpnam);

### Parse sample names
# Expects CPCT samples (CPCT........T/R)
sub parseSamples {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
}

### Somatic Variant Callers
sub runStrelka {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
}

sub runVarscan {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
}

sub runFreeBayes {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
}

### Merge results
sub mergeSomaticVariants {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
}

### Read config and check input
sub readConfiguration{
    my $configuration = shift;
    my %opt;

    foreach my $key (keys %{$configuration}){
	$opt{$key} = $configuration->{$key};
    }


    if(! $opt{GENOME}){ die "ERROR: No GENOME found in .ini file\n" }
    elsif(! -e $opt{GENOME}){ die"ERROR: $opt{GENOME} does not exist\n"}
    if(! $opt{OUTPUT_DIR}){ die "ERROR: No OUTPUT_DIR found in .conf file\n" }
    if(! $opt{MAIL}){ die "ERROR: No MAIL address specified in .conf file\n" }
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
