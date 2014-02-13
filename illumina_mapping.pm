#!/usr/bin/perl -w


##################################################################################################################################################
###This script is designed to run BWA mapping, followed my sambamba markdup. Where possible a per-sample merged BAM file will be generated
###using sambamba merge. Aside from mapping the script also runs sambamba flagstat on each *dedup.bam file.
###
###Author: S.W.Boymans
###Latest change: Created first version
###
###
##################################################################################################################################################


package illumina_mapping;

use strict;
use POSIX qw(tmpnam);

sub runMapping {
    
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
    
    my $FAI = "$opt{GENOME}\.fai";
    die "genome $opt{GENOME} does not exists!!\t$!\n" if !-e "$opt{GENOME}.bwt";
    die "fai file $FAI does not exists!!\n" if !-e $FAI;


    my $mainJobID = "$opt{OUTPUT_DIR}/".get_job_id()."_mapping_qsub.sh";

    open (my $QSUB,">$mainJobID") or die "ERROR: Couldn't create $mainJobID\n";
    print $QSUB "\#!/bin/sh\n\n. $opt{CLUSTER_PATH}/settings.sh\n\n";

    my $samples = {};
    my @jobs_to_wait;
    my $toMap = {};
    
    #Try to search for matching pairs in the input FASTQ files
    foreach my $input (keys %{$opt{FASTQ}}){
	if($input =~ m/\_R1/){
	    my $pairName = $input;
	    $pairName =~ s/\_R1/\_R2/;
	    if(exists ($opt{FASTQ}->{$pairName})){
		$toMap->{$input."-".$pairName} = 1;
	    }else{
		$toMap->{$input} = 1;
	    }
	}elsif($input =~ m/\_R2/){
	    my $pairName = $input;
	    $pairName =~ s/\_R2/\_R1/;
	    if(exists ($opt{FASTQ}->{$pairName})){
		$toMap->{$pairName."-".$input} = 1;
	    }else{
		$toMap->{$input} = 1;
	    }
	}
    }


    foreach my $input (keys %{$toMap}){
	my @files = split("-",$input);
	my $R1 = undef;
        my $R2 = undef;
	my $coreName = undef;
    
	if(scalar(@files) == 2){
	    print "Switching to paired end mode!\n";
	    $R1 = $files[0];
	    $R2 = $files[1];
	    if($R1 !~ m/fastq.gz$/ or $R2 !~ m/fastq.gz$/){
		die "ERROR: Invalid input files:\n\t$R1\n\t$R2\n";
	    }
	
	
	}elsif(scalar(@files) == 1){
	    print "Switching to fragment mode!\n";
	    $R1 = $files[0];

	    if($R1 !~ m/fastq.gz$/){
	        die "ERROR: Invalid input file:\n\t$R1\n";
	    }
	

	}else{
	    die "ERROR: Invalid input pair: $input\n";
	}

	$coreName = (split("/", $R1))[-1];
	my ($sampleName, $flowcellID, $index, $lane, $tag) =  split("_", $coreName);
	$coreName =~ s/\.fastq.gz//;
	$coreName =~ s/\_R1//;
	$coreName =~ s/\_R2//;
    
	my ($RG_PL, $RG_ID, $RG_LB, $RG_SM) = ('ILLUMINA_HISEQ', $coreName, $sampleName, $sampleName);
	
	
	print "Creating $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_dedup.bam with:\n";
		
    
	if (-e "$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_dedup.bam"){
	    warn "WARNING: $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_dedup.bam already exists, skipping\n";
	    next;
        }
    
    
	if(! -e "$opt{OUTPUT_DIR}/$sampleName"){
    	    mkdir("$opt{OUTPUT_DIR}/$sampleName") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sampleName\n";
	}
	if(! -e "$opt{OUTPUT_DIR}/$sampleName/mapping"){
    	    mkdir("$opt{OUTPUT_DIR}/$sampleName/mapping") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sampleName/mapping\n";
	}
	if(! -e "$opt{OUTPUT_DIR}/$sampleName/jobs"){
    	    mkdir("$opt{OUTPUT_DIR}/$sampleName/jobs") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sampleName/jobs\n";
	}
	if(! -e "$opt{OUTPUT_DIR}/$sampleName/logs"){
    	    mkdir("$opt{OUTPUT_DIR}/$sampleName/logs") or die "Couldn't create directory: $opt{OUTPUT_DIR}/$sampleName/logs\n";
	}

	if($opt{MAPPING_MODE} eq 'batch'){
	
	    submitBatchJobs(\%opt,$QSUB,$samples, $sampleName, $coreName, $R1, $R2);
	
	}elsif($opt{MAPPING_MODE} eq 'single'){
	    
	    submitSingleJobs(\%opt,$QSUB, $samples, $sampleName, $coreName, $R1, $R2);
	
	}
	

    }


    my $mergeJobs = {};

    #Create a merging joblist for every sample
    foreach my $sample (keys %{$samples}){
    
	my @bamList = ();
	my @jobIds = ();
	my $pass = 1;
	foreach my $chunk (@{$samples->{$sample}}){
	    push(@bamList, $chunk->{'file'});
	    push(@jobIds, $chunk->{'jobId'});
	}
	
	
	my $jobId = "MERGING_$sample\_".get_job_id();
	$mergeJobs->{$sample} = $jobId;
	
	open MERGE_SH,">$opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh\n";
	print MERGE_SH "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
	print MERGE_SH "cd $opt{OUTPUT_DIR}/$sample\n";
	print MERGE_SH "uname -n > logs/$jobId.host\n";
	
	print MERGE_SH "BAMS=( ".join(" ",@bamList)." )\n";
	print MERGE_SH "PASS=0\n";
	print MERGE_SH "for i in \"\$\{BAMS\[\@\]\}\"\n";
	print MERGE_SH "do\n";
	print MERGE_SH "\tDONEFILE=\`echo \$i | sed 's/\.bam/\.done/'\`\n";
	print MERGE_SH "\tif [ ! -f \$DONEFILE ]\n";
	print MERGE_SH "\tthen\n";
	print MERGE_SH "\t\techo \"ERROR: \$i is probably incomplete, no .done file found for it\" >> logs/merge.err\n";
	print MERGE_SH "\t\tPASS=1\n";
	print MERGE_SH "\tfi\n";
	print MERGE_SH "done\n\n";
	print MERGE_SH "if [ \$PASS -eq 1 ]\n";
	print MERGE_SH "then\n";
	print MERGE_SH "\techo \"ERROR: merging failed due to incomplete BAM-file(s)\" >> logs/merge.err\n";
	print MERGE_SH "else\n";
	
	print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba merge -t $opt{MAPPING_THREADS} mapping/$sample\_dedup.bam ".join(" ",@bamList)." 1>>logs/merge.log 2>>logs/merge.err\n";
	print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba index -t $opt{MAPPING_THREADS} mapping/$sample\_dedup.bam\n";
	print MERGE_SH "\t$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MAPPING_THREADS} mapping/$sample\_dedup.bam > mapping/$sample\_dedup.flagstat\n";
	print MERGE_SH "fi\n\n";
	
	print MERGE_SH "TOTALREADS=0\n";
	print MERGE_SH "for i in \$( find \$PWD/mapping -name '*sorted_dedup.flagstat')\n";
	print MERGE_SH "do\n";
	print MERGE_SH "\tVAL=\`grep -m 1 -P \"\\d+\" \$i | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
        print MERGE_SH "\tTOTALREADS=\$((\$TOTALREADS + \$VAL))\n";
	print MERGE_SH "done\n\n";
	
	print MERGE_SH "if [ -s mapping/$sample\_dedup.flagstat ]\n";
	print MERGE_SH "then\n";
	print MERGE_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" mapping/$sample\_dedup.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	print MERGE_SH "\tif [ \$FS1 -eq \$TOTALREADS ]\n";
	print MERGE_SH "\tthen\n";
        print MERGE_SH "\t\tfor i in \$( find \$PWD/mapping -name '*sorted_dedup.bam')\n";
        print MERGE_SH "\t\tdo\n";
        print MERGE_SH "\t\t\trm \$i\n";
        print MERGE_SH "\t\tdone\n";
	print MERGE_SH "\telse\n";
        print MERGE_SH "\t\techo \"ERROR: read counts from *sorted_dedup.flagstat files and mapping/$sample\_dedup.flagstat do not match\" >> logs/$sample\_cleanup.err\n";
        print MERGE_SH "\tfi\n";
	print MERGE_SH "else\n";
	print MERGE_SH "\techo \"ERROR: mapping/$sample\_dedup.flagstat is empty.\" >> logs/$sample\_cleanup.err\n";
	print MERGE_SH "fi\n\n";
	
	print MERGE_SH "if [ ! -s logs/$sample\_cleanup.err ]\n";
	print MERGE_SH "then\n";
	print MERGE_SH "\ttouch mapping/$sample\_dedup.done\n";
	print MERGE_SH "fi\n\n";
	
	close MERGE_SH;
    
	print $QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{MAPPING_PROJECT} -pe threaded $opt{MAPPING_THREADS} -o $opt{OUTPUT_DIR}/$sample/logs -e $opt{OUTPUT_DIR}/$sample/logs -N $jobId -hold_jid ".join(",",@jobIds)." $opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh\n\n";
    
    }

    
    close $QSUB;
    

    system("sh $mainJobID");
    
    return $mergeJobs;
}

sub submitBatchJobs{
    
    my ($opt,$QSUB ,$samples, $sampleName, $coreName, $R1, $R2) = @_;
    my %opt = %$opt;
    my $jobId = "MAPPING_$coreName\_".get_job_id();
    my ($RG_PL, $RG_ID, $RG_LB, $RG_SM) = ('ILLUMINA_HISEQ', $coreName, $sampleName, $sampleName);
    

    open BWA_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$jobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$jobId.sh\n";
    print BWA_SH "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
    print BWA_SH "mkdir $opt{CLUSTER_TMP}/$jobId/\n";
    print BWA_SH "cd $opt{CLUSTER_TMP}/$jobId/ \n";
    print BWA_SH "uname -n > $opt{OUTPUT_DIR}/$sampleName/logs/$jobId.host\n\n";
    print BWA_SH "echo \"Mapping pair\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$jobId.host\n";
    
    if($R2){
        print "\t$R1\n\t$R2\n";
        print BWA_SH "$opt{BWA_PATH}/bwa mem -t $opt{MAPPING_THREADS} -c 100 -M -R \"\@RG\tID:$RG_ID\tSM:$RG_SM\tPL:$RG_PL\tLB:$RG_LB\" $opt{GENOME} $R1 $R2 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa_log | $opt{SAMBAMBA_PATH}/sambamba view -t $opt{MAPPING_THREADS} --format=bam -S -o $coreName.bam /dev/stdin\n";
    }else{
        print "\t$R1\n";
        print BWA_SH "$opt{BWA_PATH}/bwa mem -t $opt{MAPPING_THREADS} -c 100 -M -R \"\@RG\tID:$RG_ID\tSM:$RG_SM\tPL:$RG_PL\tLB:$RG_LB\" $opt{GENOME} $R1 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa_log | $opt{SAMBAMBA_PATH}/sambamba view -t $opt{MAPPING_THREADS} --format=bam -S -o $coreName.bam /dev/stdin\n";
    }

    print BWA_SH "echo \"mapping finished\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$jobId.host\n";
    print BWA_SH "$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MAPPING_THREADS} $coreName.bam > $coreName.flagstat\n";
    print BWA_SH "echo \"flagstat $coreName.bam finished\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$jobId.host\n\n";
    
    print BWA_SH "$opt{SAMBAMBA_PATH}/sambamba sort -m $opt{MAPPING_MEM}GB -t $opt{MAPPING_THREADS} -o $coreName\_sorted.bam $coreName.bam 1>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_sort.log 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_sort.err\n";
    print BWA_SH "echo \"sortbam finished\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$jobId.host\n";
    print BWA_SH "$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MAPPING_THREADS} $coreName\_sorted.bam > $coreName\_sorted.flagstat\n";
    print BWA_SH "echo \"flagstat $coreName\_sorted.bam finished\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$jobId.host\n\n";
    
    print BWA_SH "$opt{SAMBAMBA_PATH}/sambamba index -t $opt{MAPPING_THREADS} $coreName\_sorted.bam $coreName\_sorted.bai 1>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_index.log 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_index.err\n";
    print BWA_SH "echo \"indexing finished\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$jobId.host\n\n";
        
    print BWA_SH "$opt{SAMBAMBA_PATH}/sambamba markdup -t $opt{MAPPING_THREADS} $coreName\_sorted.bam $coreName\_sorted_dedup.bam 1>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_dedup.log 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_dedup.err\n";
    print BWA_SH "echo \"markdup finished\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$jobId.host\n";
    print BWA_SH "$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MAPPING_THREADS} $coreName\_sorted_dedup.bam > $coreName\_sorted_dedup.flagstat\n";
    print BWA_SH "echo \"flagstat $coreName\_sorted_dedup.bam finished\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$jobId.host\n\n";

    print BWA_SH "if [ -s $coreName.flagstat ] && [ -s $coreName\_sorted.flagstat ]\n";
    print BWA_SH "then\n";
    print BWA_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" $coreName.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
    print BWA_SH "\tFS2=\`grep -m 1 -P \"\\d+ \" $coreName\_sorted.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
    print BWA_SH "\tif [ \$FS1 -eq \$FS2 ]\n";
    print BWA_SH "\tthen\n";
    print BWA_SH "\t\trm $coreName.bam\n";
    print BWA_SH "\telse\n";
    print BWA_SH "\t\techo \"ERROR: $coreName.flagstat and $coreName\_sorted.flagstat do not have the same read counts\" >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n";
    print BWA_SH "\tfi\n";
    print BWA_SH "else\n";
    print BWA_SH "\techo \"ERROR: Either $coreName.flagstat or $coreName\_sorted.flagstat is empty.\" >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n";
    print BWA_SH "fi\n\n";
    
    print BWA_SH "if [ -s $coreName\_sorted.flagstat ] && [ -s $coreName\_sorted_dedup.flagstat ]\n";
    print BWA_SH "then\n";
    print BWA_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" $coreName\_sorted.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
    print BWA_SH "\tFS2=\`grep -m 1 -P \"\\d+ \" $coreName\_sorted_dedup.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
    print BWA_SH "\tif [ \$FS1 -eq \$FS2 ]\n";
    print BWA_SH "\tthen\n";
    print BWA_SH "\t\trm $coreName\_sorted.bam\n";
    print BWA_SH "\t\trm $coreName\_sorted.bai\n";
    print BWA_SH "\telse\n";
    print BWA_SH "\t\techo \"ERROR: $coreName\_sorted.flagstat and $coreName\_sorted_dedup.flagstat do not have the same read counts\" >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n";
    print BWA_SH "\tfi\n";
    print BWA_SH "else\n";
    print BWA_SH "\techo \"ERROR: Either $coreName\_sorted.flagstat or $coreName\_sorted_dedup.flagstat is empty.\" >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n";
    print BWA_SH "fi\n\n";

    print BWA_SH "if [ ! -s $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err ]\n";
    print BWA_SH "then\n";
    print BWA_SH "\ttouch $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.done\n";
    print BWA_SH "fi\n\n";

    print BWA_SH "mv -f $opt{CLUSTER_TMP}/$jobId/* $opt{OUTPUT_DIR}/$sampleName/mapping/\n";
    print BWA_SH "if [ -f $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.bam ];then\n";
    print BWA_SH "\trm -r $opt{CLUSTER_TMP}/$jobId/\n";
    print BWA_SH "fi\n";
    print BWA_SH "echo \"moving results finished\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$jobId.host\n";

    close BWA_SH;
	
    push(@{$samples->{$sampleName}}, {'jobId'=>$jobId, 'file'=>"$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.bam"});
    print $QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{MAPPING_PROJECT} -pe threaded $opt{MAPPING_THREADS} -o $opt{OUTPUT_DIR}/$sampleName/logs -e $opt{OUTPUT_DIR}/$sampleName/logs -N $jobId $opt{OUTPUT_DIR}/$sampleName/jobs/$jobId.sh\n\n";

}

sub submitSingleJobs{

    my ($opt,$QSUB ,$samples, $sampleName, $coreName, $R1, $R2) = @_;
    my %opt = %$opt;
    my ($RG_PL, $RG_ID, $RG_LB, $RG_SM) = ('ILLUMINA_HISEQ', $coreName, $sampleName, $sampleName);
    
    my $mappingJobId = "MAPPING_$coreName\_".get_job_id();
    my $mappingFSJobId = "MAPPING_FS_$coreName\_".get_job_id();
    my $sortJobId = "SORT_$coreName\_".get_job_id();
    my $sortFSJobId = "SORT_FS_$coreName\_".get_job_id();
    my $indexJobId = "INDEX_$coreName\_".get_job_id();
    my $markdupJobId = "MARKDUP_$coreName\_".get_job_id();
    my $markdupFSJobId = "MARKDUP_FS_$coreName\_".get_job_id();
    my $cleanupJobId = "CLEANUP_$coreName\_".get_job_id();
    
    ###############BWA JOB###############
    open BWA_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$mappingJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$mappingJobId.sh\n";
    print BWA_SH "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
    print BWA_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
    print BWA_SH "uname -n > $opt{OUTPUT_DIR}/$sampleName/logs/$mappingJobId.host\n";
    print BWA_SH "echo \"Mapping pair\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$mappingJobId.host\n\n";
    
    if($R2){
        print "\t$R1\n\t$R2\n";
        print BWA_SH "$opt{BWA_PATH}/bwa mem -t $opt{MAPPING_THREADS} -c 100 -M -R \"\@RG\tID:$RG_ID\tSM:$RG_SM\tPL:$RG_PL\tLB:$RG_LB\" $opt{GENOME} $R1 $R2 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa_log | $opt{SAMBAMBA_PATH}/sambamba view -t $opt{MAPPING_THREADS} --format=bam -S -o $coreName.bam /dev/stdin\n";
    }else{
        print "\t$R1\n";
        print BWA_SH "$opt{BWA_PATH}/bwa mem -t $opt{MAPPING_THREADS} -c 100 -M -R \"\@RG\tID:$RG_ID\tSM:$RG_SM\tPL:$RG_PL\tLB:$RG_LB\" $opt{GENOME} $R1 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa_log | $opt{SAMBAMBA_PATH}/sambamba view -t $opt{MAPPING_THREADS} --format=bam -S -o $coreName.bam /dev/stdin\n";
    }

    print BWA_SH "echo \"mapping finished\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$mappingJobId.host\n\n";
    close BWA_SH;
    
    print $QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{MAPPING_PROJECT} -pe threaded $opt{MAPPING_THREADS} -o $opt{OUTPUT_DIR}/$sampleName/logs -e $opt{OUTPUT_DIR}/$sampleName/logs -N $mappingJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$mappingJobId.sh\n";
    ###################################
    
    ###############FLAGSTAT AFTER MAPPING JOB###############
    open FS1_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$mappingFSJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$mappingFSJobId.sh\n";
    print FS1_SH "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
    print FS1_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
    print FS1_SH "uname -n > $opt{OUTPUT_DIR}/$sampleName/logs/$mappingFSJobId.host\n";
    print FS1_SH "echo \"starting flagstat $coreName.bam\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$mappingFSJobId.host\n";
    print FS1_SH "$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MAPPING_THREADS} $coreName.bam > $coreName.flagstat\n";
    print FS1_SH "echo \"flagstat $coreName.bam finished\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$mappingFSJobId.host\n\n";
    close FS1_SH;
    
    print $QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{MAPPING_PROJECT} -pe threaded $opt{MAPPING_THREADS} -o $opt{OUTPUT_DIR}/$sampleName/logs -e $opt{OUTPUT_DIR}/$sampleName/logs -N $mappingFSJobId -hold_jid $mappingJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$mappingFSJobId.sh\n";
    ##################################
    
    ###############SORT JOB###############
    
    open SORT_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$sortJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$sortJobId.sh\n";
    print SORT_SH "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
    print SORT_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
    print SORT_SH "uname -n > $opt{OUTPUT_DIR}/$sampleName/logs/$sortJobId.host\n";
    print SORT_SH "echo \"starting sorting $coreName.bam\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sortJobId.host\n";
    print SORT_SH "$opt{SAMBAMBA_PATH}/sambamba sort -m $opt{MAPPING_MEM}GB -t $opt{MAPPING_THREADS} -o $coreName\_sorted.bam $coreName.bam 1>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_sort.log 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_sort.err\n";
    print SORT_SH "echo \"sorting finished\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sortJobId.host\n";
    close SORT_SH;
    
    print $QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{MAPPING_PROJECT} -pe threaded $opt{MAPPING_THREADS} -o $opt{OUTPUT_DIR}/$sampleName/logs -e $opt{OUTPUT_DIR}/$sampleName/logs -N $sortJobId -hold_jid $mappingJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$sortJobId.sh\n";
    ##################################
    
    ###############FLAGSTAT AFTER SORT JOB###############
    
    open FS2_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$sortFSJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$sortFSJobId.sh\n";
    print FS2_SH "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
    print FS2_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
    print FS2_SH "uname -n > $opt{OUTPUT_DIR}/$sampleName/logs/$sortJobId.host\n";
    print FS2_SH "echo \"starting flagstat $coreName\_sorted.bam\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sortFSJobId.host\n";
    print FS2_SH "$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MAPPING_THREADS} $coreName\_sorted.bam > $coreName\_sorted.flagstat\n";
    print FS2_SH "echo \"flagstat $coreName\_sorted.bam finished\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$sortFSJobId.host\n\n";
    close FS2_SH;
    
    print $QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{MAPPING_PROJECT} -pe threaded $opt{MAPPING_THREADS} -o $opt{OUTPUT_DIR}/$sampleName/logs -e $opt{OUTPUT_DIR}/$sampleName/logs -N $sortFSJobId -hold_jid $sortJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$sortFSJobId.sh\n";    
    #################################
    
    ###############INDEX JOB###############
    
    open INDEX_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$indexJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$indexJobId.sh\n";
    print INDEX_SH "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
    print INDEX_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
    print INDEX_SH "uname -n > $opt{OUTPUT_DIR}/$sampleName/logs/$indexJobId.host\n";
    print INDEX_SH "echo \"starting indexing $coreName\_sorted.bam\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$indexJobId.host\n";
    print INDEX_SH "$opt{SAMBAMBA_PATH}/sambamba index -t $opt{MAPPING_THREADS} $coreName\_sorted.bam $coreName\_sorted.bai 1>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_index.log 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_index.err\n";
    print INDEX_SH "echo \"indexing finished\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$indexJobId.host\n";
    close INDEX_SH;
    
    print $QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{MAPPING_PROJECT} -pe threaded $opt{MAPPING_THREADS} -o $opt{OUTPUT_DIR}/$sampleName/logs -e $opt{OUTPUT_DIR}/$sampleName/logs -N $indexJobId -hold_jid $sortJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$indexJobId.sh\n";    
    ################################

    ###############MARKDUP JOB###############
    open MARKDUP_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$markdupJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$markdupJobId.sh\n";
    print MARKDUP_SH "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
    print MARKDUP_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
    print MARKDUP_SH "uname -n > $opt{OUTPUT_DIR}/$sampleName/logs/$markdupJobId.host\n";
    print MARKDUP_SH "echo \"starting markdup $coreName\_sorted.bam\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$markdupJobId.host\n";
    print MARKDUP_SH "$opt{SAMBAMBA_PATH}/sambamba markdup -t $opt{MAPPING_THREADS} $coreName\_sorted.bam $coreName\_sorted_dedup.bam 1>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_dedup.log 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_dedup.err\n";
    print MARKDUP_SH "echo \"markdup finished\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$markdupJobId.host\n";
    close MARKDUP_SH;
    
    print $QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{MAPPING_PROJECT} -pe threaded $opt{MAPPING_THREADS} -o $opt{OUTPUT_DIR}/$sampleName/logs -e $opt{OUTPUT_DIR}/$sampleName/logs -N $markdupJobId -hold_jid $indexJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$markdupJobId.sh\n";
    ################################
    
    ###############FLAGSTAT AFTER MARKDUP JOB###############    
    open FS3_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$markdupFSJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$markdupFSJobId.sh\n";
    print FS3_SH "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
    print FS3_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
    print FS3_SH "uname -n > $opt{OUTPUT_DIR}/$sampleName/logs/$markdupFSJobId.host\n";
    print FS3_SH "echo \"starting flagstat $coreName\_sorted_dedup.bam\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$markdupFSJobId.host\n";
    print FS3_SH "$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MAPPING_THREADS} $coreName\_sorted_dedup.bam > $coreName\_sorted_dedup.flagstat\n";
    print FS3_SH "echo \"flagstat $coreName\_sorted_dedup.bam finished\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$markdupFSJobId.host\n";
    close FS3_SH;
    
    print $QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{MAPPING_PROJECT} -pe threaded $opt{MAPPING_THREADS} -o $opt{OUTPUT_DIR}/$sampleName/logs -e $opt{OUTPUT_DIR}/$sampleName/logs -N $markdupFSJobId -hold_jid $markdupJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$markdupFSJobId.sh\n";
    ###############################
    
    ###############CLEANUP JOB###############
    
    open CLEAN_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$cleanupJobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$cleanupJobId.sh\n";
    print CLEAN_SH "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
    print CLEAN_SH "cd $opt{OUTPUT_DIR}/$sampleName/mapping \n";
    print CLEAN_SH "uname -n > $opt{OUTPUT_DIR}/$sampleName/logs/$cleanupJobId.host\n";

    
    print CLEAN_SH "if [ -s $coreName.flagstat ] && [ -s $coreName\_sorted.flagstat ]\n";
    print CLEAN_SH "then\n";
    print CLEAN_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" $coreName.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
    print CLEAN_SH "\tFS2=\`grep -m 1 -P \"\\d+ \" $coreName\_sorted.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
    print CLEAN_SH "\tif [ \$FS1 -eq \$FS2 ]\n";
    print CLEAN_SH "\tthen\n";
    print CLEAN_SH "\t\trm $coreName.bam\n";
    print CLEAN_SH "\telse\n";
    print CLEAN_SH "\t\techo \"ERROR: $coreName.flagstat and $coreName\_sorted.flagstat do not have the same read counts\" >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n";
    print CLEAN_SH "\tfi\n";
    print CLEAN_SH "else\n";
    print CLEAN_SH "\techo \"ERROR: Either $coreName.flagstat or $coreName\_sorted.flagstat is empty.\" >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n";
    print CLEAN_SH "fi\n\n";
    
    print CLEAN_SH "if [ -s $coreName\_sorted.flagstat ] && [ -s $coreName\_sorted_dedup.flagstat ]\n";
    print CLEAN_SH "then\n";
    print CLEAN_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" $coreName\_sorted.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
    print CLEAN_SH "\tFS2=\`grep -m 1 -P \"\\d+ \" $coreName\_sorted_dedup.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
    print CLEAN_SH "\tif [ \$FS1 -eq \$FS2 ]\n";
    print CLEAN_SH "\tthen\n";
    print CLEAN_SH "\t\trm $coreName\_sorted.bam\n";
    print CLEAN_SH "\t\trm $coreName\_sorted.bai\n";
    print CLEAN_SH "\telse\n";
    print CLEAN_SH "\t\techo \"ERROR: $coreName\_sorted.flagstat and $coreName\_sorted_dedup.flagstat do not have the same read counts\" >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n";
    print CLEAN_SH "\tfi\n";
    print CLEAN_SH "else\n";
    print CLEAN_SH "\techo \"ERROR: Either $coreName\_sorted.flagstat or $coreName\_sorted_dedup.flagstat is empty.\" >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err\n";
    print CLEAN_SH "fi\n\n";

    print CLEAN_SH "if [ ! -s $opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_cleanup.err ]\n";
    print CLEAN_SH "then\n";
    print CLEAN_SH "\ttouch $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.done\n";
    print CLEAN_SH "fi\n\n";

    close CLEAN_SH;

    print $QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{MAPPING_PROJECT} -pe threaded $opt{MAPPING_THREADS} -o $opt{OUTPUT_DIR}/$sampleName/logs -e $opt{OUTPUT_DIR}/$sampleName/logs -N $cleanupJobId -hold_jid $mappingFSJobId,$sortFSJobId,$markdupFSJobId $opt{OUTPUT_DIR}/$sampleName/jobs/$cleanupJobId.sh\n\n";

    push(@{$samples->{$sampleName}}, {'jobId'=>$cleanupJobId, 'file'=>"$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.bam"});




}




sub readConfiguration {
    my $configuration = shift;
    
    my %opt = (
	
	'BWA_PATH'      	=> undef,
	'SAMBAMBA_PATH'		=> undef,
	'CLUSTER_PATH'  	=> undef,
	'CLUSTER_THREADS'	=> undef,
	'CLUSTER_MEM'		=> undef,
	'CLUSTER_TMP'		=> undef,
	'GENOME'		=> undef,
	'OUTPUT_DIR'		=> undef,
	'FASTQ'			=> []
    );
    
    foreach my $key (keys %{$configuration}){
	$opt{$key} = $configuration->{$key};
    }
    
    if(! $opt{BWA_PATH}){ die "ERROR: No BWA_PATH found in .conf file\n" }
    if(! $opt{SAMBAMBA_PATH}){ die "ERROR: No SAMBAMBA_PATH found in .conf file\n" }
    if(! $opt{MAPPING_THREADS}){ die "ERROR: No MAPPING_THREADS found in .ini file\n" }
    if(! $opt{MAPPING_MEM}){ die "ERROR: No MAPPING_MEM found in .ini file\n" }
    if(! $opt{MAPPING_QUEUE}){ die "ERROR: No MAPPING_QUEUE found in .ini file\n" }
    if(! $opt{MAPPING_PROJECT}){ die "ERROR: No MAPPING_PROJECT found in .ini file\n" }
    if(! $opt{MAPPING_MODE}){ die "ERROR: No MAPPING_MODE found in .conf file\n" }
    if(! $opt{CLUSTER_PATH}){ die "ERROR: No CLUSTER_PATH found in .conf file\n" }
    if(! $opt{CLUSTER_TMP}){ die "ERROR: No CLUSTER_TMP found in .conf file\n" }
    if(! $opt{GENOME}){ die "ERROR: No GENOME found in .conf file\n" }
    if(! $opt{OUTPUT_DIR}){ die "ERROR: No OUTPUT_DIR found in .conf file\n" }
    if(! $opt{FASTQ}){ die "ERROR: No FASTQ files specified\n" }
    
    
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