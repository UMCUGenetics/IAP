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

    open (QSUB,">$mainJobID") or die "ERROR: Couldn't create $mainJobID\n";
    print QSUB "\#!/bin/sh\n\n. $opt{CLUSTER_PATH}/settings.sh\n\n";

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
	
	    #print "R1 : $R1\n";
	    #print "R2 : $R2\n";
	
	
	}elsif(scalar(@files) == 1){
	    print "Switching to fragment mode!\n";
	    $R1 = $files[0];

	    if($R1 !~ m/fastq.gz$/){
	        die "ERROR: Invalid input file:\n\t$R1\n";
	    }
	
	    #print "R1 : $R1\n";

	}else{
	    die "ERROR: Invalid input pair: $input\n";
	}

	$coreName = (split("/", $R1))[-1];
	my ($sampleName, $flowcellID, $index, $lane, $tag) =  split("_", $coreName);
	$coreName =~ s/\.fastq.gz//;
	$coreName =~ s/\_R1//;
	$coreName =~ s/\_R2//;
    
	my ($RG_PL, $RG_ID, $RG_LB, $RG_SM) = ('ILLUMINA_HISEQ', $coreName, $sampleName, $sampleName);
	
	my $jobId = "MAPPING_$coreName\_".get_job_id();
	push(@{$samples->{$sampleName}}, {'jobId'=>$jobId, 'file'=>"$opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.bam"});
	
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

	print QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{MAPPING_PROJECT} -pe threaded $opt{MAPPING_THREADS} -o $opt{OUTPUT_DIR}/$sampleName/logs -e $opt{OUTPUT_DIR}/$sampleName/logs -N $jobId $opt{OUTPUT_DIR}/$sampleName/jobs/$jobId.sh\n\n";
    
	open BWA_SH,">$opt{OUTPUT_DIR}/$sampleName/jobs/$jobId.sh" or die "Couldn't create $opt{OUTPUT_DIR}/$sampleName/jobs/$jobId.sh\n";
	print BWA_SH "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
	print BWA_SH "mkdir $opt{CLUSTER_TMP}/$jobId/\n";
	print BWA_SH "cd $opt{CLUSTER_TMP}/$jobId/ \n\n";
	print BWA_SH "uname -n > $opt{OUTPUT_DIR}/$sampleName/logs/$jobId.host\n";
	print BWA_SH "echo \"Mapping pair\t\" `date` >> logs/$jobId.host\n";
    
	if($R2){
	    print "\t$R1\n\t$R2\n";
	    print BWA_SH "$opt{BWA_PATH}/bwa mem -t $opt{MAPPING_THREADS} -c 100 -M -R \"\@RG\tID:$RG_ID\tSM:$RG_SM\tPL:$RG_PL\tLB:$RG_LB\" $opt{GENOME} $R1 $R2 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa_log | $opt{SAMBAMBA_PATH}/sambamba view -t $opt{MAPPING_THREADS} --format=bam -S -o $coreName.bam /dev/stdin\n\n";
	}else{
	    print "\t$R1\n";
	    print BWA_SH "$opt{BWA_PATH}/bwa mem -t $opt{MAPPING_THREADS} -c 100 -M -R \"\@RG\tID:$RG_ID\tSM:$RG_SM\tPL:$RG_PL\tLB:$RG_LB\" $opt{GENOME} $R1 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_bwa_log | $opt{SAMBAMBA_PATH}/sambamba view -t $opt{MAPPING_THREADS} --format=bam -S -o $coreName.bam /dev/stdin\n\n";
	}
    
	print BWA_SH "$opt{SAMBAMBA_PATH}/sambamba sort -m $opt{MAPPING_MEM} -t $opt{MAPPING_THREADS} -o $coreName\_sorted.bam $coreName.bam 1>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_sort.log 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_sort.err\n\n";
	print BWA_SH "echo \"sortbam finished\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName.host\n";
	print BWA_SH "$opt{SAMBAMBA_PATH}/sambamba index -t $opt{MAPPING_MEM} $coreName\_sorted.bam $coreName\_sorted.bai 1>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_index.log 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_index.err\n\n";
	print BWA_SH "echo \"indexing finished\t\" `date` >> $opt{OUTPUT_DIR}/$sampleName/logs/$coreName.host\n";
	print BWA_SH "if [ -f $coreName\_sorted.bai ];then\n";
	print BWA_SH "\trm $coreName.bam\n";
	print BWA_SH "fi\n";
        
	print BWA_SH "$opt{SAMBAMBA_PATH}/sambamba markdup -t $opt{MAPPING_THREADS} $coreName\_sorted.bam $coreName\_sorted_dedup.bam 1>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_dedup.log 2>>$opt{OUTPUT_DIR}/$sampleName/logs/$coreName\_dedup.err\n\n";
	print BWA_SH "$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MAPPING_THREADS} $coreName\_sorted_dedup.bam > $coreName\_sorted_dedup.flagstat\n\n";
	print BWA_SH "if [ -f $coreName\_sorted_dedup.bam ];then\n";
	print BWA_SH "\trm $coreName\_sorted.bam\n";
	print BWA_SH "\trm $coreName\_sorted.bai\n";
	print BWA_SH "fi\n";

	print BWA_SH "mv -f $opt{CLUSTER_TMP}/$jobId/* $opt{OUTPUT_DIR}/$sampleName/mapping/\n";
	print BWA_SH "if [ -f $opt{OUTPUT_DIR}/$sampleName/mapping/$coreName\_sorted_dedup.bam ];then\n";
	print BWA_SH "\trm -r $opt{CLUSTER_TMP}/$jobId/\n";
	print BWA_SH "fi\n";

	close BWA_SH;
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
	print MERGE_SH ". $opt{CLUSTER_PATH}/settings.sh\n\n";
	print MERGE_SH "cd $opt{OUTPUT_DIR}/$sample\n";
	print MERGE_SH "uname -n > logs/$jobId.host\n";
	print MERGE_SH "$opt{SAMBAMBA_PATH}/sambamba merge -t $opt{MAPPING_THREADS} mapping/$sample\_dedup.bam ".join(" ",@bamList)." 1>>logs/merge.log 2>>logs/merge.err\n\n";
	print MERGE_SH "$opt{SAMBAMBA_PATH}/sambamba index -t $opt{MAPPING_THREADS} mapping/$sample\_dedup.bam\n\n";
	print MERGE_SH "$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{MAPPING_THREADS} mapping/$sample\_dedup.bam > mapping/$sample\_dedup.flagstat\n\n";
	close MERGE_SH;
    
	print QSUB "qsub -q $opt{MAPPING_QUEUE} -P $opt{MAPPING_PROJECT} -pe threaded $opt{MAPPING_THREADS} -o $opt{OUTPUT_DIR}/$sample/logs -e $opt{OUTPUT_DIR}/$sample/logs -N $jobId -hold_jid ".join(",",@jobIds)." $opt{OUTPUT_DIR}/$sample/jobs/$jobId.sh\n\n";
    
    }

    
    close QSUB;
    

    system("sh $mainJobID");
    
    return $mergeJobs;
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