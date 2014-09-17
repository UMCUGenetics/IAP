#!/usr/bin/perl -w

##################################################################################################################################################
###This script is designed to run GATK HC and UG variant callers using GATK Queue
###
###
###Author: R.F.Ernst
###Latest change: Added UG
###TODO:
##################################################################################################################################################

package illumina_calling;

use strict;
use POSIX qw(tmpnam);

sub runVariantCalling {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my @sampleBams;
    my @runningJobs;
    my $jobID = "VC_".get_job_id();

    ### Skip variant calling if .raw_variants.vcf already exists
    if (-e "$opt{OUTPUT_DIR}/logs/VariantCaller.done"){
	warn "WARNING: $opt{OUTPUT_DIR}/logs/VariantCaller.done exists, skipping \n";
	return $jobID;
    }
    
    ### Build Queue command
    my $javaMem = $opt{CALLING_MASTERTHREADS} * $opt{CALLING_MEM};
    my $command = "java -Xmx".$javaMem."G -Xms".$opt{CALLING_MEM}."G -jar $opt{QUEUE_PATH}/Queue.jar ";
    $command .= "-jobQueue $opt{CALLING_QUEUE} -jobNative \"-pe threaded $opt{CALLING_THREADS}\" -jobRunner GridEngine -jobReport $opt{OUTPUT_DIR}/logs/VariantCaller.jobReport.txt "; #Queue options

    ### Add caller and UG specific settings
    $command .= "-S $opt{CALLING_SCALA} ";
    if ($opt{CALLING_UGMODE}) {
	$command .= " -glm $opt{CALLING_UGMODE} ";
    }

    ### Common settings
    $command .= "-R $opt{GENOME} -O $runName -mem $opt{CALLING_MEM} -nct $opt{CALLING_THREADS} -nsc $opt{CALLING_SCATTER} -stand_call_conf $opt{CALLING_STANDCALLCONF} -stand_emit_conf $opt{CALLING_STANDEMITCONF} ";

    ### Add all bams
    foreach my $sample (@{$opt{SAMPLES}}){
	my $sampleBam;
        if ($opt{INDELREALIGNMENT} eq "yes" and $opt{BASEQUALITYRECAL} eq "yes") {
	    $sampleBam = "$sample/mapping/".$sample."_dedup_realigned_recalibrated.bam";
	} elsif ($opt{INDELREALIGNMENT} eq "yes" and $opt{BASEQUALITYRECAL} eq "no") {
	    $sampleBam = "$sample/mapping/".$sample."_dedup_realigned.bam";
	} elsif ($opt{INDELREALIGNMENT} eq "no" and $opt{BASEQUALITYRECAL} eq "yes") {
	    $sampleBam = "$sample/mapping/".$sample."_dedup_recalibrated.bam";
	} elsif ($opt{INDELREALIGNMENT} eq "no" and $opt{BASEQUALITYRECAL} eq "no") {
	    $sampleBam = "$sample/mapping/".$sample."_dedup.bam";
	}
	$command .= "-I $opt{OUTPUT_DIR}/$sampleBam ";
	push( @sampleBams, "$opt{OUTPUT_DIR}/$sampleBam");
	## Running jobs
	if ( $opt{RUNNING_JOBS}->{$sample} ){
	    push( @runningJobs, @{$opt{RUNNING_JOBS}->{$sample}} );
	}
    }

    ### Optional settings
    if ( $opt{CALLING_DBSNP} ) {
	$command .= "-D $opt{CALLING_DBSNP} ";
    }
    if ( $opt{CALLING_TARGETS} ) {
	$command .= "-L $opt{CALLING_TARGETS} ";
	if ( $opt{CALLING_INTERVALPADDING} ) {
	    $command .= "-ip $opt{CALLING_INTERVALPADDING} ";
	}
    }
    ### retry option
    if($opt{QUEUE_RETRY} eq 'yes'){
	$command  .= "-retry 1 ";
    }
    $command .= "-run";

    #Create main bash script
    my $bashFile = $opt{OUTPUT_DIR}."/jobs/VariantCalling_".$jobID.".sh";
    my $logDir = $opt{OUTPUT_DIR}."/logs";

    open CALLING_SH, ">$bashFile" or die "cannot open file $bashFile \n";
    print CALLING_SH "#!/bin/bash\n\n";
    print CALLING_SH "bash $opt{CLUSTER_PATH}/settings.sh\n\n";
    print CALLING_SH "cd $opt{OUTPUT_DIR}/tmp/\n";
    print CALLING_SH "echo \"Start variant caller\t\" `date` \"\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n\n";
    
    print CALLING_SH "if [ -f ".shift(@sampleBams)." ";
    foreach my $sampleBam (@sampleBams){
	print CALLING_SH "-a -f $sampleBam ";
    }
    print CALLING_SH "]\n";
    print CALLING_SH "then\n";
    print CALLING_SH "\t$command\n";
    print CALLING_SH "else\n";
    print CALLING_SH "\techo \"ERROR: One or more input bam files do not exist.\" >&2\n";
    print CALLING_SH "fi\n\n";
    
    print CALLING_SH "if [ -f $opt{OUTPUT_DIR}/tmp/.$runName\.raw_variants.vcf.done ]\n";
    print CALLING_SH "then\n";
    print CALLING_SH "\tmv $opt{OUTPUT_DIR}/tmp/$runName\.raw_variants.vcf $opt{OUTPUT_DIR}/\n";
    print CALLING_SH "\tmv $opt{OUTPUT_DIR}/tmp/$runName\.raw_variants.vcf.idx $opt{OUTPUT_DIR}/\n";
    print CALLING_SH "\ttouch $opt{OUTPUT_DIR}/logs/VariantCaller.done\n";
    print CALLING_SH "fi\n\n";
    print CALLING_SH "echo \"Finished variant caller\t\" `date` \"\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n";
    #Start main bash script
    if (@runningJobs){
	system "qsub -q $opt{CALLING_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{CALLING_MASTERTHREADS} -o $logDir/VariantCaller_$runName.out -e $logDir/VariantCaller_$runName.err -N $jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
    } else {
	system "qsub -q $opt{CALLING_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{CALLING_MASTERTHREADS} -o $logDir/VariantCaller_$runName.out -e $logDir/VariantCaller_$runName.err -N $jobID $bashFile";
    }
    
    return $jobID;
}

sub runVcfPrep {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];

    symlink($opt{VCF},"$opt{OUTPUT_DIR}/$runName.raw_variants.vcf");
    
    return $runName;
}

sub readConfiguration{
    my $configuration = shift;
    my %opt;

    foreach my $key (keys %{$configuration}){
	$opt{$key} = $configuration->{$key};
    }

    if(! $opt{QUEUE_PATH}){ die "ERROR: No PICARD_PATH found in .ini file\n" }
    if(! $opt{CALLING_MASTERQUEUE}){ die "ERROR: No CALLING_MASTERQUEUE found in .ini file\n" }
    if(! $opt{CALLING_MASTERTHREADS}){ die "ERROR: No CALLING_MASTERTHREADS found in .ini file\n" }
    if(! $opt{CALLING_QUEUE}){ die "ERROR: No CALLING_QUEUE found in .ini file\n" }
    if(! $opt{CALLING_THREADS}){ die "ERROR: No CALLING_THREADS found in .ini file\n" }
    if(! $opt{CALLING_MEM}){ die "ERROR: No CALLING_QUEUE found in .ini file\n" }
    if(! $opt{CALLING_SCATTER}){ die "ERROR: No CALLING_SCATTER found in .ini file\n" }
    if(! $opt{CALLING_SCALA}){ die "ERROR: No CALLING_SCALA found in .ini file\n" }
    if($opt{CALLING_UGMODE}){ 
	if($opt{CALLING_UGMODE} ne "SNP" and $opt{CALLING_UGMODE} ne "INDEL" and $opt{CALLING_UGMODE} ne "BOTH"){ die "ERROR: UGMODE: $opt{CALLING_UGMODE} does not exist use SNP, INDEL or BOTH\n"}
    }
    if(! $opt{CALLING_STANDCALLCONF}){ die "ERROR: No CALLING_STANDCALLCONF found in .ini file\n" }
    if(! $opt{CALLING_STANDEMITCONF}){ die "ERROR: No CALLING_STANDEMITCONF found in .ini file\n" }
    if(! $opt{QUEUE_RETRY}){ die "ERROR: No QUEUE_RETRY found in .ini file\n" }
    if(! $opt{GENOME}){ die "ERROR: No GENOME found in .ini file\n" }
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