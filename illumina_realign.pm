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
    
    if(! -e "$opt{OUTPUT_DIR}/tmp"){
	mkdir("$opt{OUTPUT_DIR}/tmp") or die "Couldn't create directory: $opt{OUTPUT_DIR}/tmp\n";
    }
    if(! -e "$opt{OUTPUT_DIR}/jobs"){
	mkdir("$opt{OUTPUT_DIR}/jobs") or die "Couldn't create directory: $opt{OUTPUT_DIR}/jobs\n";
    }
    if(! -e "$opt{OUTPUT_DIR}/logs"){
	mkdir("$opt{OUTPUT_DIR}/logs") or die "Couldn't create directory: $opt{OUTPUT_DIR}/logs\n";
    }



    open (QSUB,">$mainJobID") or die "ERROR: Couldn't create $mainJobID\n";
    print QSUB "\#!/bin/sh\n\n. $opt{CLUSTER_PATH}/settings.sh\n\n";
    
    #MULTI SAMPLE - MULTI OUTPUT
    if($opt{REALIGNMENT_MODE} eq 'multi'){
	my $jobId = "REALIGN_".get_job_id();
	my $cleanupJobId = "REALIGN_CLEANUP\_".get_job_id();
	my $mergeJobs = "";
	my @waitFor = ();
	
	open REALIGN_SH,">$opt{OUTPUT_DIR}/jobs/$jobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/jobs/$jobId.sh\n";
	print REALIGN_SH "\#!/bin/sh\n\n";
	print REALIGN_SH "#\$ -S /bin/sh\n";
	print REALIGN_SH ". $opt{CLUSTER_PATH}/settings.sh\n\n";
	print REALIGN_SH "cd $opt{OUTPUT_DIR}/tmp\n";
	print REALIGN_SH "uname -n > ../logs/$jobId.host\n";
	print REALIGN_SH "echo \"Starting indel realignment\t\" `date` >> ../logs/$jobId.host\n"; 
	print REALIGN_SH "java -Djava.io.tmpdir=$opt{OUTPUT_DIR}/tmp/ -jar $opt{QUEUE_PATH}/Queue.jar -R $opt{GENOME} -S $opt{REALIGNMENT_SCALA} -nt $opt{REALIGNMENT_THREADS} -mem $opt{REALIGNMENT_MEM} -nsc $opt{REALIGNMENT_SCATTER} -mode $opt{REALIGNMENT_MODE} -jobNative \"-pe threaded $opt{REALIGNMENT_THREADS}\" -run ";
	
	
	open CLEAN_SH, ">$opt{OUTPUT_DIR}/jobs/$cleanupJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/jobs/$cleanupJobId.sh\n";
	print CLEAN_SH "\#!/bin/sh\n\n";
	print CLEAN_SH "#\$ -S /bin/sh\n";	
	print CLEAN_SH "uname -n > $opt{OUTPUT_DIR}/logs/$cleanupJobId.host\n";
	print CLEAN_SH "PASS=0\n";
	
	foreach my $sample (keys %{$opt{RUNNING_JOBS}}){
	    print "\t$opt{OUTPUT_DIR}/mapping/$sample\_dedup.bam\n";
	    
	    push(@waitFor, join(",",@{$opt{RUNNING_JOBS}->{$sample}}));
	    print REALIGN_SH "-I $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.bam ";
	    
    	    my $mergeJobId = "REALIGN_MERGE_$sample\_".get_job_id();

	    
	    open MERGE_SH, ">$opt{OUTPUT_DIR}/$sample/jobs/$mergeJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sample/jobs/$mergeJobId.sh\n";
	    print MERGE_SH "\#!/bin/sh\n\n";
	    print MERGE_SH "#\$ -S /bin/sh\n";
    	    print MERGE_SH "cd $opt{OUTPUT_DIR}/tmp/.queue/\n";
	    print MERGE_SH "CHUNKS=`find \$PWD -name '*$sample\_dedup_realigned.bam' | sort | xargs`\n";
	    print MERGE_SH "if [ ! -f $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.done ]\n";
	    print MERGE_SH "then\n";
	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba merge -t $opt{REALIGNMENT_MERGETHREADS} $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.bam \`echo \$CHUNKS\` 1>>$opt{OUTPUT_DIR}/$sample/logs/realn_merge.log 2>>$opt{OUTPUT_DIR}/$sample/logs/realn_merge.err\n";
	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba index -t $opt{REALIGNMENT_MERGETHREADS} $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.bam\n";
	    print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{REALIGNMENT_MERGETHREADS} $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.bam > $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.flagstat\n\n";
	    print MERGE_SH "fi\n\n";
	    print MERGE_SH "if [ -s $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat ] && [ -s $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.flagstat ]\n";
	    print MERGE_SH "then\n";
	    print MERGE_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	    print MERGE_SH "\tFS2=\`grep -m 1 -P \"\\d+ \" $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	    print MERGE_SH "\tif [ \$FS1 -eq \$FS2 ]\n";
	    print MERGE_SH "\tthen\n";
	    print MERGE_SH "\t\ttouch $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.done\n";
	    print MERGE_SH "\telse\n";
	    print MERGE_SH "\t\techo \"ERROR: $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat and $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.flagstat do not have the same read counts\" >>$opt{OUTPUT_DIR}/$sample/logs/realn_merge.err\n";
	    print MERGE_SH "\tfi\n";
	    print MERGE_SH "else\n";
	    print MERGE_SH "\techo \"ERROR: Either $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat or $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.flagstat is empty.\" >> $opt{OUTPUT_DIR}/$sample/logs/realn_merge.err\n";
	    print MERGE_SH "fi\n";
	    close MERGE_SH;
	    
	    print CLEAN_SH "if [ ! -f $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.done ]\n";
	    print CLEAN_SH "then\n";
	    print CLEAN_SH "\tPASS=1\n";
	    print CLEAN_SH "else\n";
	    print CLEAN_SH "\techo \"ERROR: $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.bam didn't finish properly.\" >> $opt{OUTPUT_DIR}/logs/realn_cleanup.err\n";
	    print CLEAN_SH "fi\n\n";
	    
	    $mergeJobs .= "qsub -q $opt{REALIGNMENT_QUEUE} -p 100 -P $opt{REALIGNMENT_PROJECT} -pe threaded $opt{REALIGNMENT_MERGETHREADS} -o $opt{OUTPUT_DIR}/$sample/logs -e $opt{OUTPUT_DIR}/$sample/logs -N $mergeJobId -hold_jid $jobId $opt{OUTPUT_DIR}/$sample/jobs/$mergeJobId.sh\n";
	    $realignJobs->{$sample} = $mergeJobId;
	}
	
	print REALIGN_SH "-jobRunner GridEngine 1>>$opt{OUTPUT_DIR}/logs/$jobId.host 2>>$opt{OUTPUT_DIR}/logs/$jobId.host\n";
	close REALIGN_SH;
	
	print CLEAN_SH "if [ \$PASS -eq 0 ]\n";
	print CLEAN_SH "then\n";
	print CLEAN_SH "\tmv $opt{OUTPUT_DIR}/tmp/IndelRealigner.jobreport.txt $opt{OUTPUT_DIR}/logs/IndelRealigner.jobreport.txt\n";
	print CLEAN_SH "\trm -r $opt{OUTPUT_DIR}/tmp/\n";
	print CLEAN_SH "fi\n";
	close CLEAN_SH;
	
	print QSUB "qsub -q $opt{REALIGNMENT_QUEUE} -P $opt{REALIGNMENT_PROJECT} -pe threaded $opt{REALIGNMENT_THREADS} -o $opt{OUTPUT_DIR}/logs -e $opt{OUTPUT_DIR}/logs -N $jobId -hold_jid ".join(",", @waitFor)." $opt{OUTPUT_DIR}/jobs/$jobId.sh\n";
	print QSUB "qsub -q $opt{REALIGNMENT_QUEUE} -P $opt{REALIGNMENT_PROJECT} -pe threaded $opt{REALIGNMENT_THREADS} -o $opt{OUTPUT_DIR}/logs -e $opt{OUTPUT_DIR}/logs -N $cleanupJobId -hold_jid $jobId $opt{OUTPUT_DIR}/jobs/$cleanupJobId.sh\n";
	print QSUB $mergeJobs."\n";
	
    
    }
    
    
    #SINGLE SAMPLE - MULTI OUPUT
    elsif($opt{REALIGNMENT_MODE} eq 'single'){
	foreach my $sample (keys %{$opt{RUNNING_JOBS}}){
	    print "\t$opt{OUTPUT_DIR}/mapping/$sample\_dedup.bam\n";
	    my $jobId = "REALIGN_$sample\_".get_job_id();
	    
	    if(! -e "$opt{OUTPUT_DIR}/$sample/tmp"){
    		mkdir("$opt{OUTPUT_DIR}/$sample/tmp") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sample/tmp\n";
	    }
	    
	    
	    	
	    open REALIGN_SH,">$opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh\n";
	    
	    print REALIGN_SH "\#!/bin/sh\n\n";
	    print REALIGN_SH "#\$ -S /bin/sh\n";
	    print REALIGN_SH ". $opt{CLUSTER_PATH}/settings.sh\n\n";
	    print REALIGN_SH "cd $opt{OUTPUT_DIR}/$sample/tmp \n\n";
	    print REALIGN_SH "uname -n > ../logs/$jobId.host\n";
	    print REALIGN_SH "echo \"Starting indel realignment\t\" `date` >> ../logs/$jobId.host\n"; 
	    print REALIGN_SH "java -Djava.io.tmpdir=$opt{OUTPUT_DIR}/$sample/tmp -jar $opt{QUEUE_PATH}/Queue.jar -R $opt{GENOME} -S $opt{REALIGNMENT_SCALA} -nt $opt{REALIGNMENT_THREADS} -mem $opt{REALIGNMENT_MEM} -nsc $opt{REALIGNMENT_SCATTER} -mode $opt{REALIGNMENT_MODE} -jobNative \"-pe threaded $opt{REALIGNMENT_THREADS}\" -run -I $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.bam -jobRunner GridEngine 1>>../logs/realign.log 2>>../logs/realign.err\n";
	    print REALIGN_SH "if [ -f $opt{OUTPUT_DIR}/$sample/tmp/$sample\_dedup.realigned.bam ]\n";
	    print REALIGN_SH "then\n";
	    print REALIGN_SH "\t$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{REALIGNMENT_MERGETHREADS} $opt{OUTPUT_DIR}/$sample/tmp/$sample\_dedup.realigned.bam > $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.flagstat\n";
	    print REALIGN_SH "\tmv $opt{OUTPUT_DIR}/$sample/tmp/$sample\_dedup.realigned.bam $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.bam\n";
	    print REALIGN_SH "\tmv $opt{OUTPUT_DIR}/$sample/tmp/$sample\_dedup.realigned.bai $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.bai\n";
	    print REALIGN_SH "fi\n";
	    
	    print REALIGN_SH "if [ -s $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat ] && [ -s $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.flagstat ]\n";
	    print REALIGN_SH "then\n";
	    print REALIGN_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	    print REALIGN_SH "\tFS2=\`grep -m 1 -P \"\\d+ \" $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	    print REALIGN_SH "\tif [ \$FS1 -eq \$FS2 ]\n";
	    print REALIGN_SH "\tthen\n";
	    print REALIGN_SH "\t\ttouch $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.done\n";
	    print REALIGN_SH "\t\tmv $opt{OUTPUT_DIR}/$sample/tmp/IndelRealigner.jobreport.txt $opt{OUTPUT_DIR}/$sample/logs/IndelRealigner.jobreport.txt\n";
	    print REALIGN_SH "\trm -r $opt{OUTPUT_DIR}/$sample/tmp/\n";
	    print REALIGN_SH "\telse\n";
	    print REALIGN_SH "\t\techo \"ERROR: $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat and $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.flagstat do not have the same read counts\" >>../logs/realign.err\n";
	    print REALIGN_SH "\tfi\n";
	    print REALIGN_SH "else\n";
	    print REALIGN_SH "\techo \"ERROR: Either $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat or $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.flagstat is empty.\" >> ../logs/realign.err\n";
	    print REALIGN_SH "fi\n";
	    close REALIGN_SH;
	    


	    #print REALIGN_SH "cd $opt{OUTPUT_DIR}/$sample/tmp/.queue/\n";
	    #print REALIGN_SH "CHUNKS=`find \$PWD -name '*$sample\_dedup_realigned.bam' | sort | xargs`\n";
	    #print REALIGN_SH "$opt{SAMBAMBA_PATH}/sambamba merge -t $opt{REALIGNMENT_THREADS} $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.bam \`echo \$CHUNKS\` 1>>$opt{OUTPUT_DIR}/$sample/logs/realn_merge.log 2>>$opt{OUTPUT_DIR}/$sample/logs/realn_merge.err\n";
	    #print REALIGN_SH "$opt{SAMBAMBA_PATH}/sambamba index -t $opt{REALIGNMENT_THREADS} $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.bam\n";
	    #print REALIGN_SH "rm -r $opt{CLUSTER_TMP}/$jobId/\n";
	    
	    
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
    if(! $opt{REALIGNMENT_MERGETHREADS}){ die "ERROR: No REALIGNMENT_MERGETHREADS found in .conf file\n" }
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