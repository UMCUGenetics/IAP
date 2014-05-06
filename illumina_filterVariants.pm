#!/usr/bin/perl -w

##################################################################################################################################################
###This script is designed to run GATK HardFilters on .vcf files
###
###
###Author: R.F.Ernst
###Latest change: 
###TODO:
##################################################################################################################################################

package illumina_filterVariants;

use strict;
use POSIX qw(tmpnam);


sub runFilterVariants {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my @runningJobs;
    my $jobID = "FV_".get_job_id();

    ### Skip variant calling if .raw_variants.vcf already exists
    ### add .done file checking?
    if (-e "$opt{OUTPUT_DIR}/logs/VariantFilter.done"){
	warn "WARNING: $opt{OUTPUT_DIR}/logs/VariantFilter.done exists, skipping \n";
	return $jobID;
    }

    ### Build Queue command
    my $javaMem = $opt{FILTER_MASTERTHREADS} * $opt{FILTER_MEM};
    my $command = "java -Xmx".$javaMem."G -Xms".$opt{FILTER_MEM}."G -jar $opt{QUEUE_PATH}/Queue.jar "; ### Change memory allocation here!!!!!
    $command .= "-jobQueue $opt{FILTER_QUEUE} -jobNative \"-pe threaded $opt{FILTER_THREADS}\" -jobRunner GridEngine -jobReport $opt{OUTPUT_DIR}/logs/VariantFilter.jobReport.txt "; #Queue options

    ### Common settings
    $command .= "-S $opt{FILTER_SCALA} -R $opt{GENOME} -V $opt{OUTPUT_DIR}/$runName\.raw_variants.vcf -O $runName -mem $opt{FILTER_MEM} -nsc $opt{FILTER_SCATTER} -mode $opt{FILTER_MODE} ";
    
    ### Mode dependent settings
    if ($opt{FILTER_MODE} eq "SNP" || $opt{FILTER_MODE} eq "BOTH") {
	my @SNPFilterNames = split("\t",$opt{FILTER_SNPNAME});
	my @SNPFilterExprs = split("\t",$opt{FILTER_SNPEXPR});

	if (scalar(@SNPFilterNames) ne scalar(@SNPFilterExprs)) {
	    die "FILTER_SNPNAME and FILTER_SNPEXPR do not have the same length.";
	}
	foreach my $i (0 .. scalar(@SNPFilterNames)-1 ){
	    $command .= "-snpFilterName $SNPFilterNames[$i] -snpFilterExpression \"$SNPFilterExprs[$i]\" ";
	}
	if ($opt{FILTER_CLUSTERSIZE} and $opt{FILTER_CLUSTERWINDOWSIZE}){
	    $command .= "-cluster $opt{FILTER_CLUSTERSIZE} -window $opt{FILTER_CLUSTERWINDOWSIZE} ";
	}
    }

    if ($opt{FILTER_MODE} eq "INDEL" || $opt{FILTER_MODE} eq "BOTH") {
	my @INDELFilterNames = split("\t",$opt{FILTER_INDELNAME});
	my @INDELFilterExprs = split("\t",$opt{FILTER_INDELEXPR});

	if (scalar(@INDELFilterNames) ne scalar(@INDELFilterExprs)) {
	    die "FILTER_INDELNAME and FILTER_INDELEXPR do not have the same length.";
	}

	foreach my $i (0 .. scalar(@INDELFilterNames)-1 ){
	    $command .= "-indelFilterName $INDELFilterNames[$i] -indelFilterExpression \"$INDELFilterExprs[$i]\" ";
	}
    }

    $command .= "-run";

    ### Create main bash script
    my $bashFile = $opt{OUTPUT_DIR}."/jobs/FilterVariants_".$jobID.".sh";
    my $logDir = $opt{OUTPUT_DIR}."/logs";

    open FILTER_SH, ">$bashFile" or die "cannot open file $bashFile \n";
    print FILTER_SH "#!/bin/bash\n\n";
    print FILTER_SH "bash $opt{CLUSTER_PATH}/settings.sh\n\n";
    print FILTER_SH "cd $opt{OUTPUT_DIR}/tmp/\n";
    print FILTER_SH "echo \"Start variant filter\t\" `date` \"\t$runName.raw_variants.vcf\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n\n";
    
    print FILTER_SH "if [ -f $opt{OUTPUT_DIR}/$runName\.raw_variants.vcf ]\n";
    print FILTER_SH "then\n";
    print FILTER_SH "\t$command\n";
    print FILTER_SH "else\n";
    print FILTER_SH "\techo \"ERROR: $runName\.raw_variants.vcf does not exist.\" >&2\n";
    print FILTER_SH "fi\n\n";

    if ($opt{FILTER_MODE} eq "SNP"){
	print FILTER_SH "if [ -f $opt{OUTPUT_DIR}/tmp/.$runName\.filtered_snps.vcf.done ]\n";
	print FILTER_SH "then\n";
	print FILTER_SH "\tmv $opt{OUTPUT_DIR}/tmp/$runName\.filtered_snps.vcf $opt{OUTPUT_DIR}/\n";
	print FILTER_SH "\tmv $opt{OUTPUT_DIR}/tmp/$runName\.filtered_snps.vcf.idx $opt{OUTPUT_DIR}/\n";
	print FILTER_SH "\ttouch $opt{OUTPUT_DIR}/logs/VariantFilter.done\n";
	print FILTER_SH "fi\n\n";
    } elsif ($opt{FILTER_MODE} eq "INDEL"){
	print FILTER_SH "if [ -f $opt{OUTPUT_DIR}/tmp/.$runName\.filtered_indels.vcf.done ]\n";
	print FILTER_SH "then\n";
	print FILTER_SH "\tmv $opt{OUTPUT_DIR}/tmp/$runName\.filtered_indels.vcf $opt{OUTPUT_DIR}/\n";
	print FILTER_SH "\tmv $opt{OUTPUT_DIR}/tmp/$runName\.filtered_indels.vcf.idx $opt{OUTPUT_DIR}/\n";
	print FILTER_SH "\ttouch $opt{OUTPUT_DIR}/logs/VariantFilter.done\n";
	print FILTER_SH "fi\n\n";
    } elsif ($opt{FILTER_MODE} eq "BOTH"){
	print FILTER_SH "if [ -f $opt{OUTPUT_DIR}/tmp/.$runName\.filtered_variants.vcf.done ]\n";
	print FILTER_SH "then\n";
	print FILTER_SH "\tmv $opt{OUTPUT_DIR}/tmp/$runName\.filtered_variants.vcf $opt{OUTPUT_DIR}/\n";
	print FILTER_SH "\tmv $opt{OUTPUT_DIR}/tmp/$runName\.filtered_variants.vcf.idx $opt{OUTPUT_DIR}/\n";
	print FILTER_SH "\ttouch $opt{OUTPUT_DIR}/logs/VariantFilter.done\n";
	print FILTER_SH "fi\n\n";
    }
    print FILTER_SH "echo \"End variant filter\t\" `date` \"\t$runName.raw_variants.vcf\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n";
    ### Process runningjobs
    foreach my $sample (keys %{$opt{RUNNING_JOBS}}){
	push(@runningJobs, join(",",@{$opt{RUNNING_JOBS}->{$sample}}));
    }

    ### Start main bash script
    if (@runningJobs){
	system "qsub -q $opt{FILTER_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{FILTER_MASTERTHREADS} -o $logDir/VariantFilter_$runName.out -e $logDir/VariantFilter_$runName.err -N $jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
    } else {
	system "qsub -q $opt{FILTER_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{FILTER_MASTERTHREADS} -o $logDir/VariantFilter_$runName.out -e $logDir/VariantFilter_$runName.err -N $jobID $bashFile";
    }

    return $jobID;
}

sub readConfiguration{
    my $configuration = shift;
    my %opt;

    foreach my $key (keys %{$configuration}){
	$opt{$key} = $configuration->{$key};
    }

    if(! $opt{QUEUE_PATH}){ die "ERROR: No PICARD_PATH found in .ini file\n" }
    if(! $opt{FILTER_MASTERQUEUE}){ die "ERROR: No FILTER_MASTERQUEUE found in .ini file\n" }
    if(! $opt{FILTER_MASTERTHREADS}){ die "ERROR: No FILTER_MASTERTHREADS found in .ini file\n" }    
    if(! $opt{FILTER_QUEUE}){ die "ERROR: No FILTER_QUEUE found in .ini file\n" }
    if(! $opt{FILTER_THREADS}){ die "ERROR: No FILTER_THREADS found in .ini file\n" }
    if(! $opt{FILTER_MEM}){ die "ERROR: No FILTER_QUEUE found in .ini file\n" }
    if(! $opt{FILTER_SCATTER}){ die "ERROR: No FILTER_SCATTER found in .ini file\n" }
    if(! $opt{FILTER_SCALA}){ die "ERROR: No FILTER_SCALA found in .ini file\n" }
    if(! $opt{FILTER_MODE}){ die "ERROR: No FILTER_MODE  found in .ini file\n" }
    if($opt{FILTER_MODE} ne "SNP" and $opt{FILTER_MODE} ne "INDEL" and $opt{FILTER_MODE} ne "BOTH"){ die "ERROR: FILTER_MODE $opt{FILTER_MODE} does not exist use SNP, INDEL or BOTH\n"}
    if ($opt{FILTER_MODE} eq "SNP" || $opt{FILTER_MODE} eq "BOTH") {
	if(! $opt{FILTER_SNPNAME}){ die "ERROR: No FILTER_SNPNAME found in .ini file\n" }
	if(! $opt{FILTER_SNPEXPR}){ die "ERROR: No FILTER_SNPEXPR  found in .ini file\n" }
    }
    if ($opt{FILTER_MODE} eq "INDEL" || $opt{FILTER_MODE} eq "BOTH") {
	if(! $opt{FILTER_INDELNAME}){ die "ERROR: No FILTER_INDELNAME found in .ini file\n" }
	if(! $opt{FILTER_INDELEXPR}){ die "ERROR: No FILTER_INDELEXPR found in .ini file\n" }
    }
    if(! $opt{GENOME}){ die "ERROR: No GENOME found in .ini file\n" }
    if(! $opt{OUTPUT_DIR}){ die "ERROR: No OUTPUT_DIR found in .conf file\n" }
    if(! $opt{MAIL}){die "ERROR: No MAIL address specified in .conf file\n"}

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