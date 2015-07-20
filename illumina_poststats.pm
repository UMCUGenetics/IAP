#!/usr/bin/perl -w

##################################################################################################################################################
###Author: R.F.Ernst & S.W.Boymans
##################################################################################################################################################

package illumina_poststats;

use strict;
use POSIX qw(tmpnam);
use FindBin;

sub runPostStats {
    my $configuration = shift;
    my %opt = %{$configuration};
    my @runningJobs; #internal job array
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my $jobID = "PostStats_".get_job_id();
    my $jobIDCheck = "PostStats_Check_".get_job_id();

    ### Run bamMetrics
    if(! -e "$opt{OUTPUT_DIR}/logs/PostStats.done"){
	my $command = "perl $opt{BAMMETRICS_PATH}/bamMetrics.pl ";
	foreach my $sample (@{$opt{SAMPLES}}){
	    my $sampleBam = "$opt{OUTPUT_DIR}/$sample/mapping/$opt{BAM_FILES}->{$sample}";
	    $command .= "-bam $sampleBam ";
	    if (@{$opt{RUNNING_JOBS}->{$sample}}) {
		push(@runningJobs, join(",",@{$opt{RUNNING_JOBS}->{$sample}}));
	    }
	}
	$command .= "-output_dir $opt{OUTPUT_DIR}/QCStats/ ";
	$command .= "-run_name $runName ";
	$command .= "-genome $opt{GENOME} ";
	$command .= "-queue $opt{POSTSTATS_QUEUE} ";
	$command .= "-queue_threads $opt{POSTSTATS_THREADS} ";
	$command .= "-queue_mem $opt{POSTSTATS_MEM} ";
	$command .= "-queue_project $opt{CLUSTER_PROJECT} ";
	$command .= "-picard_path $opt{PICARD_PATH} ";
	#$command .= "-debug ";
	
	if ( ($opt{POSTSTATS_TARGETS}) && ($opt{POSTSTATS_BAITS}) ) {
	    $command .= "-capture ";
	    $command .= "-targets $opt{POSTSTATS_TARGETS} ";
	    $command .= "-baits $opt{POSTSTATS_BAITS} ";
	} else {
	    $command .= "-wgs ";
	    $command .= "-coverage_cap 250 "; # add option to ini file
	}
	
	if ( $opt{SINGLE_END} ) {
	    $command .= "-single_end ";
	}
	
	if ( $opt{CLUSTER_RESERVATION} eq "yes") {
	    $command .= "-queue_reserve ";
	}
	
	my $bashFile = $opt{OUTPUT_DIR}."/jobs/".$jobID.".sh";
	my $logDir = $opt{OUTPUT_DIR}."/logs";
        
	open BM_SH, ">$bashFile" or die "cannot open file $bashFile\n";
	print BM_SH "#!/bin/bash\n\n";
	print BM_SH "cd $opt{OUTPUT_DIR}\n";
	print BM_SH "echo \"Start poststats\t\" `date` \"\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n";
	print BM_SH "$command\n";
	print BM_SH "qalter -hold_jid bamMetrics_report_".$runName." $jobIDCheck\n"; #hack to make sure check does not start before bamMetrics ends.
	close BM_SH;
	
	if (@runningJobs){
	    system "qsub -q $opt{POSTSTATS_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{POSTSTATS_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $logDir/PostStats_$runName.out -e $logDir/PostStats_$runName.err -N $jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
	} else {
	    system "qsub -q $opt{POSTSTATS_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{POSTSTATS_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $logDir/PostStats_$runName.out -e $logDir/PostStats_$runName.err -N $jobID $bashFile";
	}
	
	### Check BamMetrics result
	my $bashFileCheck = $opt{OUTPUT_DIR}."/jobs/".$jobIDCheck.".sh";
	open BMCHECK_SH, ">$bashFileCheck" or die "cannot open file $bashFileCheck\n";
	print BMCHECK_SH "cd $opt{OUTPUT_DIR}\n";
	print BMCHECK_SH "if [ -f QCStats/*.bamMetrics.pdf -a -f QCStats/*.bamMetrics.html ]\nthen\n";
	print BMCHECK_SH "\ttouch logs/PostStats.done \n";
	print BMCHECK_SH "\techo \"Finished poststats\t\" `date` \"\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n";
	print BMCHECK_SH "fi\n";
	close BMCHECK_SH;

	system "qsub -q $opt{POSTSTATS_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{POSTSTATS_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $logDir/PostStats_$runName.out -e $logDir/PostStats_$runName.err -N $jobIDCheck -hold_jid bamMetrics_report_".$runName.",$jobID $bashFileCheck";
	return $jobIDCheck;

    } else {
	print "WARNING: $opt{OUTPUT_DIR}/logs/PostStats.done exists, skipping\n";
    }
}

############
sub get_job_id {
   my $id = tmpnam(); 
      $id=~s/\/tmp\/file//;
   return $id;
}

sub bashAndSubmit {
    my $command = shift;
    my $sample = shift;
    my %opt = %{shift()};
    
    my $jobID = "PostStats_".$sample."_".get_job_id();
    my $bashFile = $opt{OUTPUT_DIR}."/".$sample."/jobs/PICARD_".$sample."_".$jobID.".sh";
    my $logDir = $opt{OUTPUT_DIR}."/".$sample."/logs";
    
    open OUT, ">$bashFile" or die "cannot open file $bashFile\n";
    print OUT "#!/bin/bash\n\n";
    print OUT "cd $opt{OUTPUT_DIR}\n";
    print OUT "$command\n";
    
    if ( @{$opt{RUNNING_JOBS}->{$sample}} ){
	system "qsub -q $opt{POSTSTATS_QUEUE} -pe threaded $opt{POSTSTATS_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $logDir/PostStats_".$sample."_".$jobID.".out -e $logDir/PostStats_".$sample."_".$jobID.".err -N $jobID -hold_jid ".join(",",@{$opt{RUNNING_JOBS}->{$sample} })." $bashFile";
    } else {
	system "qsub -q $opt{POSTSTATS_QUEUE} -pe threaded $opt{POSTSTATS_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $logDir/PostStats_".$sample."_".$jobID.".out -e $logDir/PostStats_".$sample."_".$jobID.".err -N $jobID $bashFile";
    }
    return $jobID;
}

############ 

1;