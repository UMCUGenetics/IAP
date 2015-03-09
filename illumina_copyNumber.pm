#!/usr/bin/perl -w

##################################################################################################################################################
###This script is designed to run copy number tools.
###
###
###Author: R.F.Ernst
###Latest change:
###TODO:
##################################################################################################################################################

package illumina_copyNumber;

use strict;
use POSIX qw(tmpnam);
use File::Path qw(make_path);

### Parse sample names
# Expects CPCT samples (CPCT........T/R)
sub parseSamples {
    my $configuration = shift;
    my %opt = %{$configuration};
    my %somatic_samples;

    foreach my $sample (@{$opt{SAMPLES}}){
	# Parse cpct samples based on expected naming
	my ($cpct_name,$origin) = ($sample =~ /(CPCT\d{8})([TR][IVX]*$)/);
	
	if ( (! $cpct_name) || (! $origin) ){
	    warn "WARNING: $sample is not passing somatic samplename parsing, skipping \n\n";
	    next;
	}
	
	# Reference sample
	if ($origin =~ m/R.*/){
	    if ($somatic_samples{$cpct_name}{"ref"}){
		warn "\t WARNING: $cpct_name has multiple reference samples, using: $somatic_samples{$cpct_name}{'ref'} \n\n";
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

### Run and merge
sub runCopyNumberTools {
    my $configuration = shift;
    my %opt = %{$configuration};
    my @check_cnv_jobs;
    ### Loop over tumor samples
    foreach my $sample (keys(%{$opt{SOMATIC_SAMPLES}})){

	# Check correct sample ref
	if (! $opt{SOMATIC_SAMPLES}{$sample}{'ref'}){
	    warn "WARNING: No ref sample for $sample, skipping \n";
	    next;
	}

	foreach my $sample_tumor (@{$opt{SOMATIC_SAMPLES}{$sample}{'tumor'}}){
	    my @cnv_jobs;
	    ## Create output, log and job directories
	    my $sample_tumor_name = "$opt{SOMATIC_SAMPLES}{$sample}{'ref'}\_$sample_tumor";
	    my $sample_tumor_out_dir = "$opt{OUTPUT_DIR}/copyNumber/$sample_tumor_name";
	    my $sample_tumor_log_dir = "$sample_tumor_out_dir/logs/";
	    my $sample_tumor_job_dir = "$sample_tumor_out_dir/jobs/";

	    if(! -e $sample_tumor_out_dir){
		make_path($sample_tumor_out_dir) or die "Couldn't create directory:  $sample_tumor_out_dir\n";
	    }
	    if(! -e $sample_tumor_job_dir){
		make_path($sample_tumor_job_dir) or die "Couldn't create directory: $sample_tumor_job_dir\n";
	    }
	    if(! -e $sample_tumor_log_dir){
		make_path($sample_tumor_log_dir) or die "Couldn't create directory: $sample_tumor_log_dir\n";
	    }

	    ## Lookup running jobs and bams
	    my $sample_tumor_bam = "$opt{OUTPUT_DIR}/$sample_tumor/mapping/$opt{BAM_FILES}->{$sample_tumor}";
	    my @running_jobs;
	    if ( @{$opt{RUNNING_JOBS}->{$sample_tumor}} ){
		push(@running_jobs, @{$opt{RUNNING_JOBS}->{$sample_tumor}});
	    }
	    my $sample_ref_bam = "$opt{OUTPUT_DIR}/$opt{SOMATIC_SAMPLES}{$sample}{'ref'}/mapping/$opt{BAM_FILES}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}";
	    if ( @{$opt{RUNNING_JOBS}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}} ){
		push(@running_jobs, @{$opt{RUNNING_JOBS}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}});
	    }

	    ## Print sample and bam info
	    print "\n$sample \t $sample_ref_bam \t $sample_tumor_bam \n";

	    ## Skip Copy number tools if .done file exist
	    if (-e "$sample_tumor_log_dir/$sample_tumor_name.done"){
		warn "WARNING: $sample_tumor_log_dir/$sample_tumor_name.done, skipping \n";
		next;
	    }

	    ## Run CNV callers
	    if($opt{CNV_CONTRA} eq "yes"){
		print "\n###SCHEDULING CONTRA####\n";
		my $contra_job = runContra($sample_tumor, $sample_tumor_out_dir, $sample_tumor_job_dir, $sample_tumor_log_dir, $sample_tumor_bam, $sample_ref_bam, \@running_jobs, \%opt);
		if($contra_job){push(@cnv_jobs, $contra_job)};
	    }
	    
	    ## Check copy number analysis
	    my $job_id = "CHECK_".$sample_tumor."_".get_job_id();
	    my $bash_file = $sample_tumor_job_dir."/".$job_id.".sh";

	    open CHECK_SH, ">$bash_file" or die "cannot open file $bash_file \n";
	    print CHECK_SH "#!/bin/bash\n\n";
	    print CHECK_SH "echo \"Start Check\t\" `date` `uname -n` >> $sample_tumor_log_dir/check.log\n\n";
	    print CHECK_SH "if [ -f $sample_tumor_log_dir/contra.done ]\n";
	    print CHECK_SH "then\n";
	    print CHECK_SH "\ttouch $sample_tumor_log_dir/$sample_tumor_name.done\n";
	    print CHECK_SH "fi\n\n";
	    print CHECK_SH "echo \"End Check\t\" `date` `uname -n` >> $sample_tumor_log_dir/check.log\n";
	    close CHECK_SH;
	    
	    if ( @cnv_jobs ){
		system "qsub -q $opt{CNVCHECK_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{CNVCHECK_THREADS} -o $sample_tumor_log_dir -e $sample_tumor_log_dir -N $job_id -hold_jid ".join(",",@cnv_jobs)." $bash_file";
	    } else {
		system "qsub -q $opt{CNVCHECK_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{CNVCHECK_THREADS} -o $sample_tumor_log_dir -e $sample_tumor_log_dir -N $job_id $bash_file";
	    }
	    push(@check_cnv_jobs, $job_id);
	    
	}
    }

    return \@check_cnv_jobs;
}

### Somatic Variant Callers
sub runContra {
    my ($sample_tumor, $out_dir, $job_dir, $log_dir, $sample_tumor_bam, $sample_ref_bam, $running_jobs, $opt) = (@_);
    my @running_jobs = @{$running_jobs};
    my %opt = %{$opt};
    my $contra_out_dir = "$out_dir/contra";

    ## Skip Strelka if .done file exist
    if (-e "$log_dir/contra.done"){
	warn "WARNING: $log_dir/contra.done, skipping \n";
	return;
    }

    ## Create strelka bash script
    my $job_id = "CNTR_".$sample_tumor."_".get_job_id();
    my $bash_file = $job_dir."/".$job_id.".sh";

    open CONTRA_SH, ">$bash_file" or die "cannot open file $bash_file \n";
    print CONTRA_SH "#!/bin/bash\n\n";
    print CONTRA_SH "if [ -f $sample_tumor_bam -a -f $sample_ref_bam ]\n";
    print CONTRA_SH "then\n";
    print CONTRA_SH "\techo \"Start Contra\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/contra.log\n\n";

    # Run CONTRA
    print CONTRA_SH "\t$opt{CONTRA_PATH}/contra.py -s $sample_tumor_bam -c $sample_ref_bam -f $opt{GENOME} -t $opt{CONTRA_TARGETS} -o $contra_out_dir/ --sampleName $sample_tumor $opt{CONTRA_FLAGS} \n\n";
    print CONTRA_SH "\tcd $contra_out_dir\n";
    # Check contra completed
    print CONTRA_SH "\tif [ -f $contra_out_dir/table/$sample_tumor*.vcf ]\n";
    print CONTRA_SH "\tthen\n";
    print CONTRA_SH "\t\ttouch $log_dir/contra.done\n";
    print CONTRA_SH "\tfi\n\n";
    
    print CONTRA_SH "\techo \"End Contra\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/contra.log\n\n";
    
    print CONTRA_SH "else\n";
    print CONTRA_SH "\techo \"ERROR: $sample_tumor_bam or $sample_ref_bam does not exist.\" >&2\n";
    print CONTRA_SH "fi\n";
    
    close CONTRA_SH;

    ## Run job
    if ( @running_jobs ){
	system "qsub -q $opt{CONTRA_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{CONTRA_THREADS} -R $opt{CLUSTER_RESERVATION} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
    } else {
	system "qsub -q $opt{CONTRA_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{CONTRA_THREADS} -R $opt{CLUSTER_RESERVATION} -o $log_dir -e $log_dir -N $job_id $bash_file";
    }

    return $job_id;
}

############
sub get_job_id {
    my $id = tmpnam();
    $id=~s/\/tmp\/file//;
    return $id;
}
############

1;
