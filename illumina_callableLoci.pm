#!/usr/bin/perl -w

##################################################################
### illumina_callableLoci.pm
### - Run gatk CallableLoci
###
### Authors: R.F.Ernst 
##################################################################

package illumina_callableLoci;

use strict;
use POSIX qw(tmpnam);
use lib "$FindBin::Bin"; #locates pipeline directory
use illumina_sge;

sub runCallableLoci {
    ###
    # Run CallableLoci per sample
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my @callable_loci_jobs;

    foreach my $sample (@{$opt{SAMPLES}}){
	###
	# Setup sample variables
	###
	my $sample_bam = "$opt{OUTPUT_DIR}/$sample/mapping/$opt{BAM_FILES}->{$sample}";
	my $log_dir = $opt{OUTPUT_DIR}."/".$sample."/logs/";
	my $tmp_dir = $opt{OUTPUT_DIR}."/".$sample."/tmp/";
	my $job_dir = $opt{OUTPUT_DIR}."/".$sample."/jobs/";
	my $output_dir = $opt{OUTPUT_DIR}."/".$sample."/";
	my $command;
	my @running_jobs;
	
	if (-e "$log_dir/CallableLoci_$sample.done"){
	    print "WARNING: $log_dir/CallableLoci_$sample.done exists, skipping CallableLoci for $sample \n";
	} else {
	    ## Setup CallableLoci sh script
	    my $jobID = "CL_$sample\_".get_job_id();
	    my $bashFile = $job_dir.$jobID.".sh";
	    my $output_summary = $sample."_CallableLoci.txt";
	    my $output_bed = $sample."_CallableLoci.bed";

	    open CALLABLE_LOCI_SH, ">$bashFile" or die "cannot open file $bashFile \n";
	    print CALLABLE_LOCI_SH "#!/bin/bash\n\n";
	    print CALLABLE_LOCI_SH "bash $opt{CLUSTER_PATH}/settings.sh\n\n";
	    print CALLABLE_LOCI_SH "cd $tmp_dir\n\n";

	    ## Running jobs
	    if ( @{$opt{RUNNING_JOBS}->{$sample}} ){
		push( @running_jobs, @{$opt{RUNNING_JOBS}->{$sample}} );
	    }

	    ###
	    # Run CallableLoci
	    ###
	    ### Build gatk command
	    $command = "java -Djava.io.tmpdir=$opt{OUTPUT_DIR}/tmp/ -Xmx".$opt{CALLABLE_LOCI_MEM}."G -jar $opt{QUEUE_PATH}/GenomeAnalysisTK.jar ";
	    $command .= "-T CallableLoci ";
	    $command .= "-R $opt{GENOME} ";
	    $command .= "-I $sample_bam ";
	    $command .= "-o $output_bed ";
	    $command .= "-summary $output_summary ";
	    $command .= "--minBaseQuality $opt{CALLABLE_LOCI_BASEQUALITY} ";
	    $command .= "--minMappingQuality $opt{CALLABLE_LOCI_MAPQUALITY} ";
	    $command .= "--minDepth $opt{CALLABLE_LOCI_DEPTH} ";
	    $command .= "--minDepthForLowMAPQ $opt{CALLABLE_LOCI_DEPTHLOWMAPQ} ";
	    
	    if ( $opt{CALLING_TARGETS} ) {
		$command .= "-L $opt{CALLING_TARGETS} ";
		if ( $opt{CALLING_INTERVALPADDING} ) {
		    $command .= "-ip $opt{CALLING_INTERVALPADDING} ";
		}
	    }
	    #Create UG bash script
	    print CALLABLE_LOCI_SH "echo \"Start CallableLoci\t\" `date` \"\t\" `uname -n` >> $log_dir/CallableLoci_$sample.log\n";
	    print CALLABLE_LOCI_SH "if [ -s $sample_bam ]\n";
	    print CALLABLE_LOCI_SH "then\n";
	    print CALLABLE_LOCI_SH "\t$command\n";
	    print CALLABLE_LOCI_SH "else\n";
	    print CALLABLE_LOCI_SH "\techo \"ERROR: Sample bam file do not exist.\" >&2\n";
	    print CALLABLE_LOCI_SH "fi\n\n";

	    print CALLABLE_LOCI_SH "if [ -s $output_bed -a -s $output_summary ]\n";
	    print CALLABLE_LOCI_SH "then\n";
	    print CALLABLE_LOCI_SH "\tmv $output_bed $output_dir\n";
	    print CALLABLE_LOCI_SH "\tmv $output_summary $output_dir\n";
	    print CALLABLE_LOCI_SH "\ttouch $log_dir/CallableLoci_$sample.done\n";
	    print CALLABLE_LOCI_SH "fi\n";
	    print CALLABLE_LOCI_SH "echo \"Finished CallableLoci\t\" `date` \"\t\" `uname -n` >> $log_dir/CallableLoci_$sample.log\n\n";
	
	    ###
	    # Submit CallableLoci JOB
	    ###
	    my $qsub = &qsubJava(\%opt,"CALLABLE_LOCI");
	    if (@running_jobs){
		system "$qsub -o $log_dir/CallableLoci_$sample.out -e $log_dir/CallableLoci_$sample.err -N $jobID -hold_jid ".join(",",@running_jobs)." $bashFile";
	    } else {
		system "$qsub -o $log_dir/CallableLoci_$sample.out -e $log_dir/CallableLoci_$sample.err -N $jobID $bashFile";
	    }
	    push(@callable_loci_jobs, $jobID);
	}
    }
    return \@callable_loci_jobs;
}

############
sub get_job_id {
    my $id = tmpnam();
    $id=~s/\/tmp\/file//;
    return $id;
}
############

1;
