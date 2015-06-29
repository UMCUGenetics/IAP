#!/usr/bin/perl -w

##################################################################################################################################################
###This script is designed to run SnpEff annotation on .vcf files
###
###
###Author: R.F.Ernst
###Latest change: 
###TODO:
##################################################################################################################################################

package illumina_annotateVariants;

use strict;
use POSIX qw(tmpnam);

sub runAnnotateVariants {
    my $configuration = shift;
    my %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my @runningJobs;
    my $command;
    my $jobID = "AV_".get_job_id();

    ### Skip variant annotation if .done file exists.
    if (-e "$opt{OUTPUT_DIR}/logs/VariantAnnotation.done"){
	print "WARNING: $opt{OUTPUT_DIR}/logs/VariantAnnotation.done exists, skipping \n";
	return $jobID;
    }

    ### vcf file
    my $invcf;
    my $outvcf;
    if ( $opt{FILTER_VARIANTS} eq "yes" ) {
	if ( $opt{FILTER_MODE} eq "BOTH" ) { $invcf = $runName.".filtered_variants.vcf"; }
	if ( $opt{FILTER_MODE} eq "SNP" ) { $invcf = $runName.".filtered_snps.vcf"; }
	if ( $opt{FILTER_MODE} eq "INDEL" ) { $invcf = $runName.".filtered_indels.vcf"; }
    } elsif ($opt{FILTER_VARIANTS} eq "no") { $invcf = $runName.".raw_variants.vcf"; }

    ### Create main bash script
    my $bashFile = $opt{OUTPUT_DIR}."/jobs/AnnotateVariants_".$jobID.".sh";
    my $logDir = $opt{OUTPUT_DIR}."/logs";

    open ANNOTATE_SH, ">$bashFile" or die "cannot open file $bashFile \n";
    print ANNOTATE_SH "#!/bin/bash\n\n";
    print ANNOTATE_SH "bash $opt{CLUSTER_PATH}/settings.sh\n\n";
    print ANNOTATE_SH "cd $opt{OUTPUT_DIR}/\n\n";
    print ANNOTATE_SH "echo \"Start variant annotation\t\" `date` \"\t$invcf\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n\n";

    ### basic eff prediction and annotation
    if($opt{ANNOTATE_SNPEFF} eq "yes"){
	$outvcf = $invcf;
	$outvcf =~ s/.vcf/_snpEff.vcf/;
	#$command = "java -Xmx".$opt{ANNOTATE_MEM}."g -jar $opt{SNPEFF_PATH}/snpEff.jar -c $opt{SNPEFF_PATH}/snpEff.config $opt{ANNOTATE_DB} -v $invcf -o gatk $opt{ANNOTATE_FLAGS} > $outvcf\n";
	$command = "java -Xmx".$opt{ANNOTATE_MEM}."g -jar $opt{SNPEFF_PATH}/snpEff.jar -c $opt{SNPEFF_PATH}/snpEff.config $opt{ANNOTATE_DB} -v $invcf $opt{ANNOTATE_FLAGS} > $outvcf\n";
	$command .= "\t$opt{IGVTOOLS_PATH}/igvtools index $outvcf\n";
	$command .= "\trm igv.log";
	print ANNOTATE_SH "if [ -f $invcf ]\n";
	print ANNOTATE_SH "then\n";
	print ANNOTATE_SH "\t$command\n";
	print ANNOTATE_SH "else\n";
	print ANNOTATE_SH "\techo \"ERROR: $invcf does not exist.\" >&2\n";
	print ANNOTATE_SH "fi\n\n";
	$invcf = $outvcf;
    }

    #run detailed annotation (sift, polyphen gerp, phyloP)
    if($opt{ANNOTATE_SNPSIFT} eq "yes"){
	$outvcf = $invcf;
	$outvcf =~ s/.vcf/_snpSift.vcf/;
	$command = "java -Xmx".$opt{ANNOTATE_MEM}."g -jar $opt{SNPEFF_PATH}/SnpSift.jar dbnsfp -v -f $opt{ANNOTATE_FIELDS} -db $opt{ANNOTATE_DBNSFP} $invcf > $outvcf\n";
	$command .= "\t$opt{IGVTOOLS_PATH}/igvtools index $outvcf\n";
	$command .= "\trm igv.log";
	print ANNOTATE_SH "if [ -f $invcf ]\n";
	print ANNOTATE_SH "then\n";
	print ANNOTATE_SH "\t$command\n";
	print ANNOTATE_SH "else\n";
	print ANNOTATE_SH "\techo \"ERROR: $invcf does not exist.\" >&2\n";
	print ANNOTATE_SH "fi\n\n";
	if($opt{ANNOTATE_SNPEFF} eq "yes"){
	    print ANNOTATE_SH "if [ -f $outvcf ]\nthen\n\trm $invcf $invcf.idx \nfi\n\n";
	}
	$invcf = $outvcf;
    }
    
    # Add ID from database vcf, for example Cosmic
    if($opt{ANNOTATE_IDFIELD} eq "yes"){
	$outvcf = $invcf;
	my $suffix = "_$opt{ANNOTATE_IDNAME}.vcf";
	$outvcf =~ s/.vcf/$suffix/;
	$command = "java -Xmx".$opt{ANNOTATE_MEM}."g -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T VariantAnnotator -nt $opt{ANNOTATE_THREADS} -R $opt{GENOME} -o $outvcf --variant $invcf --dbsnp $opt{ANNOTATE_IDDB} --alwaysAppendDbsnpId";
	print ANNOTATE_SH "if [ -f $invcf ]\n";
	print ANNOTATE_SH "then\n";
	print ANNOTATE_SH "\t$command\n";
	print ANNOTATE_SH "else\n";
	print ANNOTATE_SH "\techo \"ERROR: $invcf does not exist.\" >&2\n";
	print ANNOTATE_SH "fi\n\n";
	if($opt{ANNOTATE_SNPSIFT} eq "yes" || $opt{ANNOTATE_SNPEFF} eq "yes"){
	    print ANNOTATE_SH "if [ -f $outvcf ]\nthen\n\trm $invcf $invcf.idx \nfi\n\n";
	}
	$invcf = $outvcf;
    }
    
    #Add GoNL annotation
    if($opt{ANNOTATE_FREQUENCIES} eq "yes"){
	$outvcf = $invcf;
	my $suffix = "_$opt{ANNOTATE_FREQNAME}.vcf";
	$outvcf =~ s/.vcf/$suffix/;
	$command = "java -Xmx".$opt{ANNOTATE_MEM}."g -jar $opt{SNPEFF_PATH}/SnpSift.jar annotate -tabix -name $opt{ANNOTATE_FREQNAME}_ -info $opt{ANNOTATE_FREQINFO} $opt{ANNOTATE_FREQDB} $invcf > $outvcf \n";
	$command .= "\t$opt{IGVTOOLS_PATH}/igvtools index $outvcf\n";
	$command .= "\trm igv.log";
	print ANNOTATE_SH "if [ -f $invcf ]\n";
	print ANNOTATE_SH "then\n";
	print ANNOTATE_SH "\t$command\n";
	print ANNOTATE_SH "else\n";
	print ANNOTATE_SH "\techo \"ERROR: $invcf does not exist.\" >&2\n";
	print ANNOTATE_SH "fi\n\n";
	if($opt{ANNOTATE_SNPSIFT} eq "yes" || $opt{ANNOTATE_SNPEFF} eq "yes" || $opt{ANNOTATE_IDFIELD} eq "yes"){
	    print ANNOTATE_SH "if [ -f $outvcf ]\nthen\n\trm $invcf $invcf.idx \nfi\n\n";
	}
	$invcf = $outvcf;
    }
    
    print ANNOTATE_SH "if [ -s $outvcf ]\nthen\n\ttouch $opt{OUTPUT_DIR}/logs/VariantAnnotation.done\nfi\n\n"; ### check whether annotated vcf is not empty
    print ANNOTATE_SH "echo \"End variant annotation\t\" `date` \"\t$invcf\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n";
    
    ### Process runningjobs
    foreach my $sample (@{$opt{SAMPLES}}){
	if( exists $opt{RUNNING_JOBS}->{$sample} && @{$opt{RUNNING_JOBS}->{$sample}} ) {
	    push(@runningJobs, join(",",@{$opt{RUNNING_JOBS}->{$sample}}));
	}
    }

    ### Start main bash script
    if (@runningJobs){
	system "qsub -q $opt{ANNOTATE_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{ANNOTATE_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $logDir/VariantAnnotation_$runName.out -e $logDir/VariantAnnotation_$runName.err -N $jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
    } else {
	system "qsub -q $opt{ANNOTATE_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{ANNOTATE_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $logDir/VariantAnnotation_$runName.out -e $logDir/VariantAnnotation_$runName.err -N $jobID $bashFile";
    }

    return $jobID;
}

############
sub get_job_id {
    my $id = tmpnam();
    $id=~s/\/tmp\/file//;
    return $id;
}
############
1;