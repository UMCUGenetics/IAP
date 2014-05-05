#!/usr/bin/perl -w

##################################################################################################################################################
###This script is designed to run GATK baseRecalibration using GATK Queue
###
###
###Author: R.F.Ernst
###Latest change: Created skeleton
###
###
##################################################################################################################################################

package illumina_baseRecal;

use strict;
use POSIX qw(tmpnam);

sub runBaseRecalibration {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
    my %baseRecalJobs;
    
    print "Running base recalibration for the following BAM-files:\n";
    
    foreach my $sample (@{$opt{SAMPLES}}){
	my $jobID = "BR_".$sample."_".get_job_id();
	
	### Check input .bam files
	my $inBam;
	my $outBam;
	if($opt{INDELREALIGNMENT} eq 'yes'){ 
	    $inBam = $opt{OUTPUT_DIR}."/".$sample."/mapping/".$sample."\_dedup_realigned.bam"; # expect $SAMPLE_dedup_realigned.bam
	    $outBam = $sample."\_dedup_realigned_recalibrated";
	}
	else {
	    $inBam = $opt{OUTPUT_DIR}."/".$sample."/mapping/".$sample."\_dedup.bam"; #expect $sample_dedup.bam
	    $outBam = $sample."\_dedup_recalibrated";
	}
	print "\t".$inBam."\n";
	
	### Check output .bam files
	if (-e "$opt{OUTPUT_DIR}/$sample/logs/BaseRecalibration_$sample.done"){
	    warn "\t WARNING: $opt{OUTPUT_DIR}/$sample/logs/BaseRecalibration_$sample.done exists, skipping \n";
	    next;
	}
	
	### Build Queue command
	my $javaMem = $opt{BASERECALIBRATION_MASTERTHREADS} * $opt{BASERECALIBRATION_MEM};
	my $command = "java -Xmx".$javaMem."G -Xms".$opt{BASERECALIBRATION_MEM}."G -jar $opt{QUEUE_PATH}/Queue.jar ";
	# cluster options
	$command .= "-jobQueue $opt{BASERECALIBRATION_QUEUE} -jobNative \"-pe threaded $opt{BASERECALIBRATION_THREADS}\" -jobRunner GridEngine -jobReport $opt{OUTPUT_DIR}/$sample/logs/BaseRecalibration.jobReport.txt "; #Queue options
	# baseRecalibration options
	$command .= "-S $opt{BASERECALIBRATION_SCALA} -R $opt{GENOME} -I $inBam -mem $opt{BASERECALIBRATION_MEM} -nct $opt{BASERECALIBRATION_THREADS} -nsc $opt{BASERECALIBRATION_SCATTER} ";
	
	### Parsing known files and add them to $command.
	my @knownFiles;
	if($opt{BASERECALIBRATION_KNOWN}) {
	    @knownFiles = split('\t', $opt{BASERECALIBRATION_KNOWN});
	    foreach my $knownFile (@knownFiles) { $command .= "-knownSites $knownFile "; }
	}
	
	$command .= "-run";

	### Create bash script
	my $bashFile = $opt{OUTPUT_DIR}."/".$sample."/jobs/".$jobID.".sh";
	my $logDir = $opt{OUTPUT_DIR}."/".$sample."/logs";

	open BASERECAL_SH, ">$bashFile" or die "cannot open file $bashFile \n";
	print BASERECAL_SH "#!/bin/bash\n\n";
	print BASERECAL_SH "bash $opt{CLUSTER_PATH}/settings.sh\n\n";
	print BASERECAL_SH "cd $opt{OUTPUT_DIR}/$sample/tmp/\n";
	print BASERECAL_SH "echo \"Start base recalibration\t\" `date` \"\t$inBam\t\" `uname -n` >> ../logs/$sample.log\n";
	print BASERECAL_SH "$command\n\n";
	
	### Generate FlagStats if gatk .done file present
	print BASERECAL_SH "if [ -f $opt{OUTPUT_DIR}/$sample/tmp/.$outBam\.bam.done ]\n";
	print BASERECAL_SH "then\n";
	print BASERECAL_SH "\t$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{BASERECALIBRATION_THREADS} $opt{OUTPUT_DIR}/$sample/tmp/$outBam\.bam > $opt{OUTPUT_DIR}/$sample/mapping/$outBam\.flagstat\n";
	print BASERECAL_SH "fi\n";
	
	### Check FlagStats and move files if correct else print error
	print BASERECAL_SH "if [ -s $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat ] && [ -s $opt{OUTPUT_DIR}/$sample/mapping/$outBam\.flagstat ]\n";
	print BASERECAL_SH "then\n";
	print BASERECAL_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	print BASERECAL_SH "\tFS2=\`grep -m 1 -P \"\\d+ \" $opt{OUTPUT_DIR}/$sample/mapping/$outBam\.flagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	print BASERECAL_SH "\tif [ \$FS1 -eq \$FS2 ]\n";
	print BASERECAL_SH "\tthen\n";
	print BASERECAL_SH "\t\tmv $opt{OUTPUT_DIR}/$sample/tmp/$outBam\.bam $opt{OUTPUT_DIR}/$sample/mapping/\n";
	print BASERECAL_SH "\t\tmv $opt{OUTPUT_DIR}/$sample/tmp/$outBam\.bai $opt{OUTPUT_DIR}/$sample/mapping/\n";
	print BASERECAL_SH "\t\ttouch $opt{OUTPUT_DIR}/$sample/logs/BaseRecalibration_$sample.done\n";
	print BASERECAL_SH "\telse\n";
	print BASERECAL_SH "\t\techo \"ERROR: $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat and $opt{OUTPUT_DIR}/$sample/mapping/$outBam\.flagstat do not have the same read counts\" >>../logs/BaseRecalibration_$sample.err\n";
	print BASERECAL_SH "\tfi\n";
	print BASERECAL_SH "else\n";
	print BASERECAL_SH "\techo \"ERROR: Either $opt{OUTPUT_DIR}/$sample/mapping/$sample\_dedup.flagstat or $opt{OUTPUT_DIR}/$sample/mapping/$outBam\.flagstat is empty.\" >> ../logs/BaseRecalibration_$sample.err\n";
	print BASERECAL_SH "fi\n";
	print BASERECAL_SH "echo \"End base recalibration\t\" `date` \"\t$inBam\t\" `uname -n` >> ../logs/$sample.log\n";
	close BASERECAL_SH;
	
	### Submit bash script
	if ( $opt{RUNNING_JOBS}->{$sample} ){
	    system "qsub -q $opt{BASERECALIBRATION_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{BASERECALIBRATION_MASTERTHREADS} -o $logDir/BaseRecalibration_$sample.out -e $logDir/BaseRecalibration_$sample.err -N $jobID -hold_jid ".join(",",@{$opt{RUNNING_JOBS}->{$sample}})." $bashFile";
	} else {
	    system "qsub -q $opt{BASERECALIBRATION_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{BASERECALIBRATION_MASTERTHREADS} -o $logDir/BaseRecalibration_$sample.out -e $logDir/BaseRecalibration_$sample.err -N $jobID $bashFile";
	}
	$baseRecalJobs{$sample} = $jobID;
    }
    return %baseRecalJobs;
}


sub readConfiguration{
    my $configuration = shift;
    my %opt;

    foreach my $key (keys %{$configuration}){
	$opt{$key} = $configuration->{$key};
    }

    if(! $opt{SAMBAMBA_PATH}){ die "ERROR: No SAMBAMBA_PATH found in .ini file\n" }
    if(! $opt{BASERECALIBRATION_MASTERQUEUE}){ die "ERROR: No BASERECALIBRATION_QUEUE found in .ini file\n" }
    if(! $opt{BASERECALIBRATION_MASTERTHREADS}){ die "ERROR: No BASERECALIBRATION_THREADS found in .ini file\n" }
    if(! $opt{BASERECALIBRATION_QUEUE}){ die "ERROR: No BASERECALIBRATION_QUEUE found in .ini file\n" }
    if(! $opt{BASERECALIBRATION_THREADS}){ die "ERROR: No BASERECALIBRATION_THREADS found in .ini file\n" }
    if(! $opt{BASERECALIBRATION_MEM}){ die "ERROR: No BASERECALIBRATION_MEM found in .ini file\n" }
    if(! $opt{BASERECALIBRATION_SCALA}){ die "ERROR: No BASERECALIBRATION_SCALA found in .ini file\n" }
    if(! $opt{BASERECALIBRATION_SCATTER}){ die "ERROR: No BASERECALIBRATION_SCATTER found in .ini file\n" }
    if(! $opt{CLUSTER_PATH}){ die "ERROR: No CLUSTER_PATH found in .ini file\n" }
    if(! $opt{GENOME}){ die "ERROR: No GENOME found in .ini file\n" }
    if(! $opt{OUTPUT_DIR}){ die "ERROR: No OUTPUT_DIR found in .conf file\n" }
    if(! $opt{SAMPLES}){ die "ERROR: No SAMPLES found\n" }
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