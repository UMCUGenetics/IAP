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
    my %opt = %{readConfiguration($configuration)};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my @runningJobs;
    my $command;
    my $jobID = "AV_".get_job_id();
    
    ### Skip variant annotation if .done file exists.
    if (-e "$opt{OUTPUT_DIR}/logs/VariantAnnotation.done"){
	warn "WARNING: $opt{OUTPUT_DIR}/logs/VariantAnnotation.done exists, skipping \n";
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
	$command = "java -Xmx".$opt{ANNOTATE_MEM}."g -jar $opt{SNPEFF_PATH}/snpEff.jar -c $opt{SNPEFF_PATH}/snpEff.config $opt{ANNOTATE_DB} -v $invcf -o gatk $opt{ANNOTATE_FLAGS} > $outvcf";
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
	$command = "java -Xmx".$opt{ANNOTATE_MEM}."g -jar $opt{SNPEFF_PATH}/SnpSift.jar dbnsfp -f $opt{ANNOTATE_FIELDS} -v $opt{ANNOTATE_DBNSFP} $invcf > $outvcf";
	print ANNOTATE_SH "if [ -f $invcf ]\n";
	print ANNOTATE_SH "then\n";
	print ANNOTATE_SH "\t$command\n";
	print ANNOTATE_SH "else\n";
	print ANNOTATE_SH "\techo \"ERROR: $invcf does not exist.\" >&2\n";
	print ANNOTATE_SH "fi\n\n";
	if($opt{ANNOTATE_SNPEFF} eq "yes"){
	    print ANNOTATE_SH "if [ -f $outvcf ]\nthen\n\trm $invcf\nfi\n\n";
	}
	$invcf = $outvcf;
    }
    
    #Add GoNL annotation
    if($opt{ANNOTATE_VCFTOOLS} eq "yes"){
	$outvcf = $invcf;
	$outvcf =~ s/.vcf/_GoNL.vcf/;
	$command = "$opt{VCFTOOLS_PATH}/vcf-annotate -a $opt{ANNOTATE_FREQ} -c $opt{ANNOTATE_COLUMNS} -d $opt{ANNOTATE_DESCR} $invcf > $outvcf";
	print ANNOTATE_SH "if [ -f $invcf ]\n";
	print ANNOTATE_SH "then\n";
	print ANNOTATE_SH "\t$command\n";
	print ANNOTATE_SH "else\n";
	print ANNOTATE_SH "\techo \"ERROR: $invcf does not exist.\" >&2\n";
	print ANNOTATE_SH "fi\n\n";
	if($opt{ANNOTATE_SNPSIFT} eq "yes" || $opt{ANNOTATE_SNPEFF} eq "yes"){
	    print ANNOTATE_SH "if [ -f $outvcf ]\nthen\n\trm $invcf\nfi\n\n";
	}
    }
    print ANNOTATE_SH "if [ -s $outvcf ]\nthen\n\ttouch $opt{OUTPUT_DIR}/logs/VariantAnnotation.done\nfi\n\n"; ### check whether annotated vcf is not empty
    print ANNOTATE_SH "echo \"End variant annotation\t\" `date` \"\t$invcf\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n";
    
    ### Process runningjobs
    foreach my $sample (@{$opt{SAMPLES}}){
	if(exists $opt{RUNNING_JOBS}->{$sample}) {
	    push(@runningJobs, join(",",@{$opt{RUNNING_JOBS}->{$sample}}));
	}
    }

    ### Start main bash script
    if (@runningJobs){
	system "qsub -q $opt{ANNOTATE_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{ANNOTATE_THREADS} -o $logDir/VariantAnnotation_$runName.out -e $logDir/VariantAnnotation_$runName.err -N $jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
    } else {
	system "qsub -q $opt{ANNOTATE_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{ANNOTATE_THREADS} -o $logDir/VariantAnnotation_$runName.out -e $logDir/VariantAnnotation_$runName.err -N $jobID $bashFile";
    }

    return $jobID;
}

sub readConfiguration{
    my $configuration = shift;

    my %opt;

    foreach my $key (keys %{$configuration}){
	$opt{$key} = $configuration->{$key};
    }
    if(! $opt{SNPEFF_PATH}){ die "ERROR: No SNPEFF_PATH found in .ini file\n" }
    if(! $opt{ANNOTATE_QUEUE}){ die "ERROR: No ANNOTATE_QUEUE found in .ini file\n" }
    if(! $opt{ANNOTATE_THREADS}){ die "ERROR: No ANNOTATE_THREADS found in .ini file\n" }
    if(! $opt{ANNOTATE_MEM}){ die "ERROR: No ANNOTATE_MEM found in .ini file\n" }
    if(! $opt{ANNOTATE_SNPEFF}){ die "ERROR: No ANNOTATE_SNPEFF found in .ini file\n" }
    if($opt{ANNOTATE_SNPEFF} eq "yes"){
	if(! $opt{ANNOTATE_DB}){ die "ERROR: No ANNOTATE_DB found in .ini file\n" }
	if(! $opt{ANNOTATE_FLAGS}){ die "ERROR: No ANNOTATE_FLAGS found in .ini file\n" }
    }
    if(! $opt{ANNOTATE_SNPSIFT}){ die "ERROR: No ANNOTATE_SNPSIFT found in .ini file\n" }
    if($opt{ANNOTATE_SNPSIFT} eq "yes"){
	if(! $opt{ANNOTATE_DBNSFP}){ die "ERROR: No ANNOTATE_DBNSFP found in .ini file\n" }
	elsif( $opt{ANNOTATE_DBNSFP} && ! -e $opt{ANNOTATE_DBNSFP}) { die"ERROR: $opt{ANNOTATE_DBNSFP} does not exist\n" }
	if(! $opt{ANNOTATE_FIELDS}){ die "ERROR: No ANNOTATE_FIELDS found in .ini file\n" }
    }
    if(! $opt{ANNOTATE_VCFTOOLS}){ die "ERROR: No ANNOTATE_VCFTOOLS found in .ini file\n" }
    if($opt{ANNOTATE_VCFTOOLS} eq "yes"){
	if(! $opt{ANNOTATE_FREQ}){ die "ERROR: No ANNOTATE_FREQ found in .ini file\n" }
	elsif( $opt{ANNOTATE_FREQ} && ! -e $opt{ANNOTATE_FREQ}) { die"ERROR: $opt{ANNOTATE_FREQ} does not exist\n" }
	if(! $opt{ANNOTATE_DESCR}){ die "ERROR: No ANNOTATE_DESCR found in .ini file\n" }
	elsif( $opt{ANNOTATE_DESCR} && ! -e $opt{ANNOTATE_DESCR}) { die"ERROR: $opt{ANNOTATE_DESCR} does not exist\n" }
	if(! $opt{ANNOTATE_COLUMNS}){ die "ERROR: No ANNOTATE_COLUMNS found in .ini file\n" }
    }
    if(! $opt{OUTPUT_DIR}){ die "ERROR: No OUTPUT_DIR found in .conf file\n" }

    if(! $opt{MAIL}){die "ERROR: No MAIL address found in .ini file \n" }

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