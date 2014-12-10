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
use File::Path qw(make_path);

### Parse sample names
# Expects CPCT samples (CPCT........T/R)
# Creates directory structure
sub parseSamples {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
    my %somatic_samples;

    foreach my $sample (@{$opt{SAMPLES}}){
	# Parse cpct samples based on expected naming
	my ($cpct_name,$origin) = ($sample =~ /(CPCT\d{8})([TR].*)/);

	# Reference sample
	if ($origin =~ m/R.*/){
	    if ($somatic_samples{$cpct_name}{"ref"}){
		warn "\t WARNING: $cpct_name has multiple reference samples, using: $somatic_samples{$cpct_name}{'ref'} \n";
	    } else {
		$somatic_samples{$cpct_name}{"ref"} = $sample;
	    }
	}

	# Tumor samples
	elsif ($origin =~ m/T.*/){
	    push(@{$somatic_samples{$cpct_name}{"tumor"}},$sample);
	}
    }

    $opt{SOMATIC_SAMPLES} = {%somatic_samples};
    return \%opt;
}

### Somatic Variant Callers
sub runSomaticVariantCallers {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
    my @somvar_jobs;
    
    if($opt{SOMVAR_STRELKA} eq "yes"){
	print "\n###SCHEDULING STRELKA####\n";
	my @strelka_jobs = @{runStrelka(\%opt)};
	push(@somvar_jobs, @strelka_jobs);
    }
    if($opt{SOMVAR_VARSCAN} eq "yes"){
	print "\n###SCHEDULING VARSCAN####\n";
    }
    if($opt{SOMVAR_FREEBAYES} eq "yes"){
	print "\n###SCHEDULING FREEBAYES####\n";
    }
}


sub runStrelka {
    my %opt = %{$_[0]};
    my @strelka_jobs;
    # Create strelka output, job and log dirs
    my $job_dir = "$opt{OUTPUT_DIR}/somaticVariants/strelka/jobs";
    my $log_dir = "$opt{OUTPUT_DIR}/somaticVariants/strelka/logs";
    if(! -e "$opt{OUTPUT_DIR}/somaticVariants/strelka"){
	make_path("$opt{OUTPUT_DIR}/somaticVariants/strelka") or die "Couldn't create directory: $opt{OUTPUT_DIR}/somaticVariants/strelka\n";
    }
    if(! -e $job_dir){
	make_path($job_dir) or die "Couldn't create directory: $job_dir\n";
    }
    if(! -e $log_dir){
	make_path($log_dir) or die "Couldn't create directory: $log_dir\n";
    }

    foreach my $sample (keys(%{$opt{SOMATIC_SAMPLES}})){
	foreach my $sample_tumor (@{$opt{SOMATIC_SAMPLES}{$sample}{'tumor'}}){
	    # Lookup bams and running jobs
	    my $sample_tumor_bam = "$opt{OUTPUT_DIR}/$sample_tumor/mapping/$opt{BAM_FILES}->{$sample_tumor}";
	    my @running_jobs;
	    if ( @{$opt{RUNNING_JOBS}->{$sample_tumor}} ){
		push(@running_jobs, @{$opt{RUNNING_JOBS}->{$sample_tumor}});
	    }

	    my $sample_ref_bam = "$opt{OUTPUT_DIR}/$opt{SOMATIC_SAMPLES}{$sample}{'ref'}/mapping/$opt{BAM_FILES}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}";
	    if ( @{$opt{RUNNING_JOBS}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}} ){
		push(@running_jobs, @{$opt{RUNNING_JOBS}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}});
	    }

	    print "$sample \t $sample_ref_bam \t $sample_tumor_bam \n";

	    # Create strelka bash script
	    my $job_id = "STR_".$sample_tumor."_".get_job_id();
	    my $bash_file = $job_dir."/".$job_id.".sh";
	    my $sample_tumor_ref_out = "$opt{OUTPUT_DIR}/somaticVariants/strelka/$opt{SOMATIC_SAMPLES}{$sample}{'ref'}\_$sample_tumor";

	    open STRELKA_SH, ">$bash_file" or die "cannot open file $bash_file \n";
	    print STRELKA_SH "#!/bin/bash\n\n";
	    print STRELKA_SH "$opt{STRELKA_PATH}/bin/configureStrelkaWorkflow.pl --tumor $sample_tumor_bam --normal $sample_ref_bam --ref $opt{GENOME} --config $opt{STRELKA_INI} --output-dir $sample_tumor_ref_out\n";
	    print STRELKA_SH "cd $sample_tumor_ref_out\n";
	    print STRELKA_SH "make -j 8\n";

	    # Run job
	    if ( @running_jobs ){
		system "qsub -q $opt{STRELKA_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{STRELKA_THREADS} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
		push(@strelka_jobs, $job_id);
	    } else {
		system "qsub -q $opt{STRELKA_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{STRELKA_THREADS} -o $log_dir -e $log_dir -N $job_id $bash_file";
		push(@strelka_jobs, $job_id);
	    }
	}
    }
    return \@strelka_jobs;
}

sub runVarscan {

}

sub runFreeBayes {

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
