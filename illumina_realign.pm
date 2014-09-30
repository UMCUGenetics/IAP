#!/usr/bin/perl -w

##################################################################################################################################################
###This script is designed to run GATK indelrealignment using GATK Queue
###
###
###Author: S.W.Boymans
###Latest change: Created skeleton
###
###TODO: Add known snp files!!!
##################################################################################################################################################

package illumina_realign;

use strict;
use POSIX qw(tmpnam);

sub runRealignment {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
    my $realignJobs = {};
    my $javaMem = $opt{REALIGNMENT_MASTERTHREADS} * $opt{REALIGNMENT_MEM};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];

    print "Running $opt{REALIGNMENT_MODE} sample indel realignment for the following BAM-files:\n";
    
    ### Parsing known indel files
    my @knownIndelFiles;
    if($opt{REALIGNMENT_KNOWN}) {
	@knownIndelFiles = split('\t', $opt{REALIGNMENT_KNOWN});
    }
    
    my $mainJobID = "$opt{OUTPUT_DIR}/jobs/RealignMainJob_".get_job_id().".sh";
    open (QSUB,">$mainJobID") or die "ERROR: Couldn't create $mainJobID\n";
    print QSUB "\#!/bin/sh\n\n. $opt{CLUSTER_PATH}/settings.sh\n\n";
    
    #MULTI SAMPLE - MULTI OUTPUT
    if($opt{REALIGNMENT_MODE} eq 'multi'){
	my $jobId = "RE_".get_job_id();
	my $cleanupJobId = "REALIGN_CLEANUP\_".get_job_id();
	my $mergeJobs = "";
	my @waitFor = ();
	
	open REALIGN_SH,">$opt{OUTPUT_DIR}/jobs/$jobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/jobs/$jobId.sh\n";
	print REALIGN_SH "\#!/bin/sh\n\n";
	print REALIGN_SH ". $opt{CLUSTER_PATH}/settings.sh\n\n";
	print REALIGN_SH "cd $opt{OUTPUT_DIR}/tmp\n";
	print REALIGN_SH "uname -n > ../logs/$jobId.host\n";
	print REALIGN_SH "echo \"Start indel realignment\t\" `date` >> ../logs/$runName.log\n";
	print REALIGN_SH "java -Djava.io.tmpdir=$opt{OUTPUT_DIR}/tmp/ -Xmx".$javaMem."G -jar $opt{QUEUE_PATH}/Queue.jar -R $opt{GENOME} -S $opt{REALIGNMENT_SCALA} -jobQueue $opt{REALIGNMENT_QUEUE} -nt $opt{REALIGNMENT_THREADS} -mem $opt{REALIGNMENT_MEM} -nsc $opt{REALIGNMENT_SCATTER} -mode $opt{REALIGNMENT_MODE} -jobNative \"-pe threaded $opt{REALIGNMENT_THREADS}\" -run ";
	
	if($opt{REALIGNMENT_KNOWN}) {
	    foreach my $knownIndelFile (@knownIndelFiles) {
		print REALIGN_SH "-known $knownIndelFile ";
	    }
	}
	
	if($opt{QUEUE_RETRY} eq 'yes'){
	    print REALIGN_SH "-retry 1 ";
	}
	
	open CLEAN_SH, ">$opt{OUTPUT_DIR}/jobs/$cleanupJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/jobs/$cleanupJobId.sh\n";
	print CLEAN_SH "\#!/bin/sh\n\n";
	print CLEAN_SH "uname -n > $opt{OUTPUT_DIR}/logs/$cleanupJobId.host\n";
	print CLEAN_SH "PASS=0\n";
	
	foreach my $sample (@{$opt{SAMPLES}}){
	    print "\t$opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.bam\n";
	    
	    ## Check for realigned bam file, skip sample if realigned bam file already exist.
	    if (-e "$opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done"){
		warn "\t WARNING: $opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done exists, skipping\n";
		next;
	    }

	    push(@waitFor, join(",",@{$opt{RUNNING_JOBS}->{$sample}}));
	    print REALIGN_SH "-I $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.bam ";

	    my $mergeJobId = "REALIGN_MERGE_$sample\_".get_job_id();

	    open MERGE_SH, ">$opt{OUTPUT_DIR}/$sample/jobs/$mergeJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sample/jobs/$mergeJobId.sh\n";
	    print MERGE_SH "\#!/bin/sh\n\n";
    	    print MERGE_SH "cd $opt{OUTPUT_DIR}/tmp/.queue/\n";
	    print MERGE_SH "CHUNKS=`find \$PWD -name '*$sample\_dedup_realigned.bam' | sort | xargs`\n";
	    print MERGE_SH "if [ ! -f $opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done ]\n";
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
	    print MERGE_SH "\t\ttouch $opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done\n";
	    print MERGE_SH "\telse\n";
	    print MERGE_SH "\t\techo \"ERROR: $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat and $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.flagstat do not have the same read counts\" >>$opt{OUTPUT_DIR}/$sample/logs/realn_merge.err\n";
	    print MERGE_SH "\tfi\n";
	    print MERGE_SH "else\n";
	    print MERGE_SH "\techo \"ERROR: Either $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat or $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.flagstat is empty.\" >> $opt{OUTPUT_DIR}/$sample/logs/realn_merge.err\n";
	    print MERGE_SH "fi\n";
	    close MERGE_SH;
	    
	    print CLEAN_SH "if [ ! -f $opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done ]\n";
	    print CLEAN_SH "then\n";
	    print CLEAN_SH "\tPASS=1\n";
	    print CLEAN_SH "else\n";
	    print CLEAN_SH "\techo \"ERROR: $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.bam didn't finish properly.\" >> $opt{OUTPUT_DIR}/logs/realn_cleanup.err\n";
	    print CLEAN_SH "fi\n\n";
	    
	    $mergeJobs .= "qsub -q $opt{REALIGNMENT_QUEUE} -p 100 -P $opt{REALIGNMENT_PROJECT} -m a -M $opt{MAIL} -pe threaded $opt{REALIGNMENT_MERGETHREADS} -o $opt{OUTPUT_DIR}/$sample/logs -e $opt{OUTPUT_DIR}/$sample/logs -N $mergeJobId -hold_jid $jobId $opt{OUTPUT_DIR}/$sample/jobs/$mergeJobId.sh\n";
	    $realignJobs->{$sample} = $mergeJobId;
	}
	
	print REALIGN_SH "-jobRunner GridEngine 1>>$opt{OUTPUT_DIR}/logs/$jobId.host 2>>$opt{OUTPUT_DIR}/logs/$jobId.host\n";
	close REALIGN_SH;
	
	print CLEAN_SH "if [ \$PASS -eq 0 ]\n";
	print CLEAN_SH "then\n";
	print CLEAN_SH "echo \"Finished indel realignment\t\" `date` >> ../logs/$runName.log\n";
	print CLEAN_SH "\tmv $opt{OUTPUT_DIR}/tmp/IndelRealigner.jobreport.txt $opt{OUTPUT_DIR}/logs/IndelRealigner.jobreport.txt\n";
	print CLEAN_SH "fi\n";
	close CLEAN_SH;
	
	print QSUB "qsub -q $opt{REALIGNMENT_MASTERQUEUE} -P $opt{REALIGNMENT_PROJECT} -m a -M $opt{MAIL} -pe threaded $opt{REALIGNMENT_MASTERTHREADS} -o $opt{OUTPUT_DIR}/logs -e $opt{OUTPUT_DIR}/logs -N $jobId -hold_jid ".join(",", @waitFor)." $opt{OUTPUT_DIR}/jobs/$jobId.sh\n";
	print QSUB "qsub -q $opt{REALIGNMENT_MASTERQUEUE} -P $opt{REALIGNMENT_PROJECT} -m a -M $opt{MAIL} -pe threaded $opt{REALIGNMENT_MASTERTHREADS} -o $opt{OUTPUT_DIR}/logs -e $opt{OUTPUT_DIR}/logs -N $cleanupJobId -hold_jid $jobId $opt{OUTPUT_DIR}/jobs/$cleanupJobId.sh\n";
	print QSUB $mergeJobs."\n";
	
    
    }
    
    #SINGLE SAMPLE - MULTI OUPUT
    elsif($opt{REALIGNMENT_MODE} eq 'single'){
	foreach my $sample (@{$opt{SAMPLES}}){
	    print "\t$opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.bam\n";
	    
	    ## Check for realigned bam file, skip sample if realigned bam file already exist.
	    if (-e "$opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done"){
		warn "\t WARNING: $opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done exists, skipping\n";
		next;
	    }
	    
	    my $jobId = "Realign_$sample\_".get_job_id();

	    open REALIGN_SH,">$opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh\n";
	    
	    print REALIGN_SH "\#!/bin/bash\n\n";
	    print REALIGN_SH ". $opt{CLUSTER_PATH}/settings.sh\n\n";
	    print REALIGN_SH "cd $opt{OUTPUT_DIR}/$sample/tmp \n\n";
	    print REALIGN_SH "echo \"Start indel realignment\t\" `date` \"\t$sample\_dedup.bam\t\" `uname -n` >> ../logs/$sample.log\n\n";
	    
	    print REALIGN_SH "if [ -f $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.bam ]\n";
	    print REALIGN_SH "then\n";
	    print REALIGN_SH "\tjava -Djava.io.tmpdir=$opt{OUTPUT_DIR}/$sample/tmp -Xmx".$javaMem."G -jar $opt{QUEUE_PATH}/Queue.jar -R $opt{GENOME} -S $opt{REALIGNMENT_SCALA} -jobQueue $opt{REALIGNMENT_QUEUE} -nt $opt{REALIGNMENT_THREADS} -mem $opt{REALIGNMENT_MEM} -nsc $opt{REALIGNMENT_SCATTER} -mode $opt{REALIGNMENT_MODE} -jobNative \"-pe threaded $opt{REALIGNMENT_THREADS}\" ";
	    
	    if($opt{REALIGNMENT_KNOWN}) {
		foreach my $knownIndelFile (@knownIndelFiles) {
		    if(! -e $knownIndelFile){ die"ERROR: $knownIndelFile does not exist\n" }
		    else { print REALIGN_SH "-known $knownIndelFile " }
		}
	    }
	    
	    if($opt{QUEUE_RETRY} eq 'yes'){
		print REALIGN_SH "-retry 1 ";
	    }
	    
	    print REALIGN_SH "-run -I $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.bam -jobRunner GridEngine\n";
	    print REALIGN_SH "else\n";
	    print REALIGN_SH "echo \"ERROR: $sample\_dedup.bam does not exist.\" >&2\n";
	    print REALIGN_SH "fi\n\n";

	    print REALIGN_SH "if [ -f $opt{OUTPUT_DIR}/$sample/tmp/$sample\_dedup.realigned.bam ]\n";
	    print REALIGN_SH "then\n";
	    print REALIGN_SH "\t$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{REALIGNMENT_THREADS} $opt{OUTPUT_DIR}/$sample/tmp/$sample\_dedup.realigned.bam > $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.flagstat\n";
	    print REALIGN_SH "\tmv $opt{OUTPUT_DIR}/$sample/tmp/$sample\_dedup.realigned.bam $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.bam\n";
	    print REALIGN_SH "\tmv $opt{OUTPUT_DIR}/$sample/tmp/$sample\_dedup.realigned.bai $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.bai\n";
	    print REALIGN_SH "fi\n\n";

	    print REALIGN_SH "if [ -s $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat ] && [ -s $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.flagstat ]\n";
	    print REALIGN_SH "then\n";
	    print REALIGN_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	    print REALIGN_SH "\tFS2=\`grep -m 1 -P \"\\d+ \" $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	    print REALIGN_SH "\tif [ \$FS1 -eq \$FS2 ]\n";
	    print REALIGN_SH "\tthen\n";
	    print REALIGN_SH "\t\ttouch $opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.done\n";
	    print REALIGN_SH "\t\tmv $opt{OUTPUT_DIR}/$sample/tmp/IndelRealigner.jobreport.txt $opt{OUTPUT_DIR}/$sample/logs/IndelRealigner.jobreport.txt\n";
	    print REALIGN_SH "\telse\n";
	    print REALIGN_SH "\t\techo \"ERROR: $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat and $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.flagstat do not have the same read counts\" >>../logs/Realignment_$sample.err\n";
	    print REALIGN_SH "\tfi\n";
	    print REALIGN_SH "else\n";
	    print REALIGN_SH "\techo \"ERROR: Either $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat or $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup_realigned.flagstat is empty.\" >> ../logs/Realignment_$sample.err\n";
	    print REALIGN_SH "fi\n\n";
	    
	    print REALIGN_SH "echo \"End indel realignment\t\" `date` \"\t$sample\_dedup.bam\t\" `uname -n` >> ../logs/$sample.log\n"; 
	    close REALIGN_SH;
	    
	    if ( $opt{RUNNING_JOBS}->{$sample} ){
		print QSUB "qsub -q $opt{REALIGNMENT_MASTERQUEUE} -P $opt{REALIGNMENT_PROJECT} -m a -M $opt{MAIL} -pe threaded $opt{REALIGNMENT_MASTERTHREADS} -o $opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.out -e $opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.err -N $jobId -hold_jid ".join(",",@{$opt{RUNNING_JOBS}->{$sample}})." $opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh\n";
	    } else {
		print QSUB "qsub -q $opt{REALIGNMENT_MASTERQUEUE} -P $opt{REALIGNMENT_PROJECT} -m a -M $opt{MAIL} -pe threaded $opt{REALIGNMENT_MASTERTHREADS} -o $opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.out -e $opt{OUTPUT_DIR}/$sample/logs/Realignment_$sample.err -N $jobId $opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh\n";
	    }
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
    my %opt;

    foreach my $key (keys %{$configuration}){
	$opt{$key} = $configuration->{$key};
    }

    if(! $opt{SAMBAMBA_PATH}){ die "ERROR: No SAMBAMBA_PATH found in .ini file\n" }
    if(! $opt{REALIGNMENT_MASTERQUEUE}){ die "ERROR: No REALIGNMENT_MASTERQUEUE found in .ini file\n" }
    if(! $opt{REALIGNMENT_MASTERTHREADS}){ die "ERROR: No REALIGNMENT_MASTERTHREADS found in .ini file\n" }
    if(! $opt{REALIGNMENT_QUEUE}){ die "ERROR: No REALIGNMENT_QUEUE found in .ini file\n" }
    if(! $opt{REALIGNMENT_PROJECT}){ die "ERROR: No REALIGNMENT_PROJECT found in .ini file\n" }
    if(! $opt{REALIGNMENT_THREADS}){ die "ERROR: No REALIGNMENT_THREADS found in .ini file\n" }
    if(! $opt{REALIGNMENT_MERGETHREADS}){ die "ERROR: No REALIGNMENT_MERGETHREADS found in .ini file\n" }
    if(! $opt{REALIGNMENT_MEM}){ die "ERROR: No REALIGNMENT_MEM found in .ini file\n" }
    if(! $opt{REALIGNMENT_SCALA}){ die "ERROR: No REALIGNMENT_SCALA found in .ini file\n" }
    if(! $opt{REALIGNMENT_SCATTER}){ die "ERROR: No REALIGNMENT_SCATTER found in .ini file\n" }
    if(! $opt{REALIGNMENT_MODE}){ die "ERROR: No REALIGNMENT_MODE found in .ini file\n" }
    if(! $opt{QUEUE_RETRY}){ die "ERROR: No QUEUE_RETRY found in .ini file\n" }
    if(! $opt{CLUSTER_PATH}){ die "ERROR: No CLUSTER_PATH found in .ini file\n" }
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