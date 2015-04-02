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
    my %opt = %{$configuration};

    print "Running base recalibration for the following BAM-files:\n";
    
    foreach my $sample (@{$opt{SAMPLES}}){
	my $jobID = "BR_".$sample."_".get_job_id();
	
	### Check input .bam files
	my $inBam = $opt{BAM_FILES}->{$sample};
	(my $inFlagstat = $inBam) =~ s/bam/flagstat/;
	(my $outBam = $inBam) =~ s/bam/recalibrated.bam/;
	(my $outBai = $inBam) =~ s/bam/recalibrated.bai/;
	(my $outFlagstat = $inBam) =~ s/bam/recalibrated.flagstat/;
	
	print "\t$opt{OUTPUT_DIR}/$sample/mapping/$inBam\n";
	$opt{BAM_FILES}->{$sample} = $outBam;
	
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
	$command .= "-S $opt{BASERECALIBRATION_SCALA} -R $opt{GENOME} -I $opt{OUTPUT_DIR}/$sample/mapping/$inBam -mem $opt{BASERECALIBRATION_MEM} -nct $opt{BASERECALIBRATION_THREADS} -nsc $opt{BASERECALIBRATION_SCATTER} ";
	
	### Parsing known files and add them to $command.
	my @knownFiles;
	if($opt{BASERECALIBRATION_KNOWN}) {
	    @knownFiles = split('\t', $opt{BASERECALIBRATION_KNOWN});
	    foreach my $knownFile (@knownFiles) {
		if(! -e $knownFile){ die"ERROR: $knownFile does not exist\n" }
		else { $command .= "-knownSites $knownFile " }
	    }
	}
	### retry option
	if($opt{QUEUE_RETRY} eq 'yes'){
	    $command  .= "-retry 1 ";
	}
	$command .= "-run";

	### Create bash script
	my $bashFile = $opt{OUTPUT_DIR}."/".$sample."/jobs/".$jobID.".sh";
	my $logDir = $opt{OUTPUT_DIR}."/".$sample."/logs";

	open BASERECAL_SH, ">$bashFile" or die "cannot open file $bashFile \n";
	print BASERECAL_SH "#!/bin/bash\n\n";
	print BASERECAL_SH "bash $opt{CLUSTER_PATH}/settings.sh\n\n";
	print BASERECAL_SH "cd $opt{OUTPUT_DIR}/$sample/tmp/\n";
	print BASERECAL_SH "echo \"Start base recalibration\t\" `date` \"\t$inBam\t\" `uname -n` >> ../logs/$sample.log\n\n";
	
	print BASERECAL_SH "if [ -f $opt{OUTPUT_DIR}/$sample/mapping/$inBam ]\n";
	print BASERECAL_SH "then\n";
	print BASERECAL_SH "\t$command\n";
	print BASERECAL_SH "else\n";
	print BASERECAL_SH "\techo \"ERROR: $opt{OUTPUT_DIR}/$sample/mapping/$inBam does not exist.\" >&2\n";
	print BASERECAL_SH "fi\n\n";

	### Generate FlagStats if gatk .done file present
	print BASERECAL_SH "if [ -f $opt{OUTPUT_DIR}/$sample/tmp/.$outBam.done ]\n";
	print BASERECAL_SH "then\n";
	print BASERECAL_SH "\t$opt{SAMBAMBA_PATH}/sambamba flagstat -t $opt{BASERECALIBRATION_THREADS} $opt{OUTPUT_DIR}/$sample/tmp/$outBam > $opt{OUTPUT_DIR}/$sample/mapping/$outFlagstat\n";
	print BASERECAL_SH "fi\n\n";

	### Check FlagStats and move files if correct else print error
	print BASERECAL_SH "if [ -s $opt{OUTPUT_DIR}/$sample/mapping/$inFlagstat ] && [ -s $opt{OUTPUT_DIR}/$sample/mapping/$outFlagstat ]\n";
	print BASERECAL_SH "then\n";
	print BASERECAL_SH "\tFS1=\`grep -m 1 -P \"\\d+ \" $opt{OUTPUT_DIR}/$sample/mapping/$inFlagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	print BASERECAL_SH "\tFS2=\`grep -m 1 -P \"\\d+ \" $opt{OUTPUT_DIR}/$sample/mapping/$outFlagstat | awk '{{split(\$0,columns , \"+\")} print columns[1]}'\`\n";
	print BASERECAL_SH "\tif [ \$FS1 -eq \$FS2 ]\n";
	print BASERECAL_SH "\tthen\n";
	print BASERECAL_SH "\t\tmv $opt{OUTPUT_DIR}/$sample/tmp/$outBam $opt{OUTPUT_DIR}/$sample/mapping/\n";
	print BASERECAL_SH "\t\tmv $opt{OUTPUT_DIR}/$sample/tmp/$outBai $opt{OUTPUT_DIR}/$sample/mapping/\n";
	print BASERECAL_SH "\t\ttouch $opt{OUTPUT_DIR}/$sample/logs/BaseRecalibration_$sample.done\n";
	print BASERECAL_SH "\telse\n";
	print BASERECAL_SH "\t\techo \"ERROR: $opt{OUTPUT_DIR}/$sample/mapping/$inFlagstat and $opt{OUTPUT_DIR}/$sample/mapping/$outFlagstat do not have the same read counts\" >>../logs/BaseRecalibration_$sample.err\n";
	print BASERECAL_SH "\tfi\n";
	print BASERECAL_SH "else\n";
	print BASERECAL_SH "\techo \"ERROR: Either $opt{OUTPUT_DIR}/$sample/mapping/$inFlagstat or $opt{OUTPUT_DIR}/$sample/mapping/$outFlagstat is empty.\" >> ../logs/BaseRecalibration_$sample.err\n";
	print BASERECAL_SH "fi\n\n";
	print BASERECAL_SH "echo \"End base recalibration\t\" `date` \"\t$inBam\t\" `uname -n` >> ../logs/$sample.log\n";
	close BASERECAL_SH;
	
	### Submit bash script
	if ( @{$opt{RUNNING_JOBS}->{$sample}} ){
	    system "qsub -q $opt{BASERECALIBRATION_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{BASERECALIBRATION_MASTERTHREADS} -o $logDir/BaseRecalibration_$sample.out -e $logDir/BaseRecalibration_$sample.err -N $jobID -hold_jid ".join(",",@{$opt{RUNNING_JOBS}->{$sample}})." $bashFile";
	} else {
	    system "qsub -q $opt{BASERECALIBRATION_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{BASERECALIBRATION_MASTERTHREADS} -o $logDir/BaseRecalibration_$sample.out -e $logDir/BaseRecalibration_$sample.err -N $jobID $bashFile";
	}

	push(@{$opt{RUNNING_JOBS}->{$sample}}, $jobID);
    }
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