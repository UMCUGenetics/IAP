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


package illumina_realign;

use strict;
use POSIX qw(tmpnam);


sub runRealignment {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
    my $realignJobs = {};

    print "Running $opt{REALIGNMENT_MODE} sample indel realignment for the following BAM-files:\n";
    
    my $mainJobID = "$opt{OUTPUT_DIR}/".get_job_id()."_realign_qsub.sh";

    open (QSUB,">$mainJobID") or die "ERROR: Couldn't create $mainJobID\n";
    print QSUB "\#!/bin/sh\n\n. $opt{CLUSTER_PATH}/settings.sh\n\n";
    
    #MULTI SAMPLE - MULTI OUTPUT
    if($opt{REALIGNMENT_MODE} eq 'multi'){
	my $jobId = "REALIGN_".get_job_id();
	my $mergeJobs = "";
	my @waitFor = ();
	
	open REALIGN_SH,">$opt{OUTPUT_DIR}/$jobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$jobId.sh\n";
	print REALIGN_SH "\#!/bin/sh\n\n. $opt{CLUSTER_PATH}/settings.sh\n\n";
	print REALIGN_SH "mkdir $opt{CLUSTER_TMP}/$jobId/\n";
	print REALIGN_SH "cd $opt{CLUSTER_TMP}/$jobId/ \n\n";
	print REALIGN_SH "uname -n > $opt{OUTPUT_DIR}/$jobId.host\n";
	print REALIGN_SH "echo \"Starting indel realignment\t\" `date` >> $opt{OUTPUT_DIR}/$jobId.host\n"; 
	print REALIGN_SH "java -Djava.io.tmpdir=$opt{CLUSTER_TMP}/$jobId/ -jar $opt{QUEUE_PATH}/Queue.jar -R $opt{GENOME} -S $opt{REALIGNMENT_SCALA} -threads $opt{REALIGNMENT_THREADS} -maxmem $opt{REALIGNMENT_MEM} -scatter $opt{REALIGNMENT_SCATTER} -jobRunner GridEngine -jobNative \"-pe threaded $opt{REALIGNMENT_THREADS}\" ";
	print REALIGN_SH "mv -f $opt{CLUSTER_TMP}/$jobId/* $opt{OUTPUT_DIR}";
	
	
	foreach my $sample (keys %{$opt{RUNNING_JOBS}}){
	    print "\t$opt{OUTPUT_DIR}/mapping/$sample\_dedup.bam\n";
	    
	    push(@waitFor, join(",",@{$opt{RUNNING_JOBS}->{$sample}}));
	    print REALIGN_SH "-I $opt{OUTPUT_DIR}/mapping/$sample\_dedup.bam ";
	    
    	    my $mergeJobId = "REALIGN_MERGE_$sample\_".get_job_id();
	    
	    open MERGE_SH, ">$opt{OUTPUT_DIR}/$sample/jobs/$mergeJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sample/jobs/$mergeJobId.sh\n";
    	    print MERGE_SH "cd $opt{OUTPUT_DIR}/.queue/\n";
	    print MERGE_SH "CHUNKS=`find \$PWD -name '*$sample\_dedup_realigned.bam' | sort | xargs`\n";
	    print MERGE_SH "$opt{SAMBAMBA_PATH}/sambamba merge -t $opt{REALIGNMENT_THREADS} $opt{OUTPUT_DIR}/mapping/$sample\_dedup_realigned.bam \$CHUNKS 1>>$opt{OUTPUT_DIR}/$sample/logs/realn_merge.log 2>>$opt{OUTPUT_DIR}/$sample/logs/realn_merge.err\n";
	    print MERGE_SH "$opt{SAMBAMBA_PATH}/sambamba index -t $opt{REALIGNMENT_THREADS} $opt{OUTPUT_DIR}/mapping/$sample\_dedup_realigned.bam\n";
	    close MERGE_SH;
	    
	    $mergeJobs .= "qsub -q $opt{REALIGNMENT_QUEUE} -p 100 -P $opt{REALIGNMENT_PROJECT} -pe threaded $opt{REALIGNMENT_THREADS} -o $opt{OUTPUT_DIR}/$sample/logs -e $opt{OUTPUT_DIR}/$sample/logs -N $mergeJobId -hold_jid $jobId $opt{OUTPUT_DIR}/$sample/jobs/$mergeJobId.sh\n";
	    $realignJobs->{$sample} = $mergeJobId;
	}
	
	print REALIGN_SH "1>>$opt{OUTPUT_DIR}/$jobId.host 2>>$opt{OUTPUT_DIR}/$jobId.host\n";
	
		
	#print REALIGN_SH "if [ -f $opt{OUTPUT_DIR}/mapping/$sample\_dedup_realigned.bam ];then\n";
	print REALIGN_SH "rm -r $opt{CLUSTER_TMP}/$jobId/\n";
	#print REALIGN_SH "fi\n";
	close REALIGN_SH;
	
	print QSUB "qsub -q $opt{REALIGNMENT_QUEUE} -P $opt{REALIGNMENT_PROJECT} -pe threaded $opt{REALIGNMENT_THREADS} -o $opt{OUTPUT_DIR} -e $opt{OUTPUT_DIR} -N $jobId -hold_jid ".join(",", @waitFor)." $opt{OUTPUT_DIR}/$jobId.sh\n";
	print QSUB $mergeJobs."\n";
	
    
    }
    
    
    
    
    #SINGLE SAMPLE - MULTI OUPUT
    elsif($opt{REALIGNMENT_MODE} eq 'single'){
	foreach my $sample (keys %{$opt{RUNNING_JOBS}}){
	    print "\t$opt{OUTPUT_DIR}/mapping/$sample\_dedup.bam\n";
	    my $jobId = "REALIGN_$sample\_".get_job_id();
	    
	    	
	    open REALIGN_SH,">$opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh\n";
	    print REALIGN_SH "\#!/bin/sh\n\n. $opt{CLUSTER_PATH}/settings.sh\n\n";
	    print REALIGN_SH "mkdir $opt{CLUSTER_TMP}/$jobId/\n";
	    print REALIGN_SH "cd $opt{CLUSTER_TMP}/$jobId/ \n\n";
	    print REALIGN_SH "uname -n > $opt{OUTPUT_DIR}/$sample/logs/$jobId.host\n";
	    print REALIGN_SH "echo \"Starting indel realignment\t\" `date` >> $opt{OUTPUT_DIR}/$sample/logs/$jobId.host\n"; 
	    print REALIGN_SH "java -Djava.io.tmpdir=$opt{CLUSTER_TMP}/$jobId/ -jar $opt{QUEUE_PATH}/Queue.jar -R $opt{GENOME} -S $opt{REALIGNMENT_SCALA} -threads $opt{REALIGNMENT_THREADS} -maxmem $opt{REALIGNMENT_MEM} -scatter $opt{REALIGNMENT_SCATTER} -jobRunner GridEngine -jobNative \"-pe threaded $opt{REALIGNMENT_THREADS}\" -I $opt{OUTPUT_DIR}/mapping/$sample\_dedup.bam\n";    

	    print REALIGN_SH "cd $opt{CLUSTER_TMP}/$jobId/.queue/\n";
	    print REALIGN_SH "CHUNKS=`find \$PWD -name '*$sample\_dedup_realigned.bam' | sort | xargs`\n";
	    print REALIGN_SH "$opt{SAMBAMBA_PATH}/sambamba merge -t $opt{REALIGNMENT_THREADS} $opt{OUTPUT_DIR}/mapping/$sample\_dedup_realigned.bam \$CHUNKS 1>>$opt{OUTPUT_DIR}/$sample/logs/realn_merge.log 2>>$opt{OUTPUT_DIR}/$sample/logs/realn_merge.err\n";
	    print REALIGN_SH "$opt{SAMBAMBA_PATH}/sambamba index -t $opt{REALIGNMENT_THREADS} $opt{OUTPUT_DIR}/mapping/$sample\_dedup_realigned.bam\n";
	    print REALIGN_SH "rm -r $opt{CLUSTER_TMP}/$jobId/\n";
	    close REALIGN_SH;
	    
	    print QSUB "qsub -q $opt{REALIGNMENT_QUEUE} -P $opt{REALIGNMENT_PROJECT} -pe threaded $opt{REALIGNMENT_THREADS} -o $opt{OUTPUT_DIR}/$sample/logs -e $opt{OUTPUT_DIR}/$sample/logs -N $jobId -hold_jid ".join(",",@{$opt{RUNNING_JOBS}->{$sample}})." $opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh\n";
	    
	    $realignJobs->{$sample} = $jobId;
	}
	
    }else{
	die "ERROR: Invalid REALIGNMENT_MODE $opt{REALIGNMENT_MODE} , use 'single' or 'multi'\n";
    
    }
    
    system("sh $mainJobID");
    
    return $realignJobs;
}

sub readConfiguration{
    my $configuration = shift;
    
    my %opt = (
	
	'QUALIMAP_PATH'		=> undef,
	'SAMBAMBA_PATH'		=> undef,
	'CLUSTER_PATH'  	=> undef,
	'REALIGNMENT_THREADS'	=> undef,
	'REALIGNMENT_MEM'	=> undef,
	'REALIGNMENT_QUEUE'	=> undef,
	'REALIGNMENT_PROJECT'	=> undef,
	'REALIGNMENT_SCALA'	=> undef,
	'REALIGNMENT_SCATTER'	=> undef,
	'REALIGNMENT_MODE'	=> undef,
	'CLUSTER_TMP'		=> undef,
	'GENOME'		=> undef,
	'OUTPUT_DIR'		=> undef,
	'RUNNING_JOBS'		=> {} #do not use in .conf file
    );

    foreach my $key (keys %{$configuration}){
	$opt{$key} = $configuration->{$key};
    }

    if(! $opt{SAMBAMBA_PATH}){ die "ERROR: No SAMBAMBA_PATH found in .conf file\n" }
    if(! $opt{REALIGNMENT_PROJECT}){ die "ERROR: No REALIGNMENT_PROJECT found in .conf file\n" }
    if(! $opt{REALIGNMENT_THREADS}){ die "ERROR: No REALIGNMENT_THREADS found in .conf file\n" }
    if(! $opt{REALIGNMENT_MEM}){ die "ERROR: No REALIGNMENT_MEM found in .conf file\n" }
    if(! $opt{REALIGNMENT_QUEUE}){ die "ERROR: No REALIGNMENT_QUEUE found in .conf file\n" }
    if(! $opt{REALIGNMENT_SCALA}){ die "ERROR: No REALIGNMENT_SCALA found in .conf file\n" }
    if(! $opt{REALIGNMENT_SCATTER}){ die "ERROR: No REALIGNMENT_SCATTER found in .conf file\n" }
    if(! $opt{REALIGNMENT_MODE}){ die "ERROR: No REALIGNMENT_MODE found in .conf file\n" }
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