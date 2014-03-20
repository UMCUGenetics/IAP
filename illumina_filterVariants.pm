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
    if ($opt{FILTER_MODE} eq "SNP" && -e "$opt{OUTPUT_DIR}/$runName.filtered_snps.vcf"){
	warn "WARNING: $opt{OUTPUT_DIR}/$runName.filtered_snps.vcf already exists, skipping \n";
	return $jobID;
    } elsif ($opt{FILTER_MODE} eq "INDEL" && -e "$opt{OUTPUT_DIR}/$runName.filtered_indels.vcf"){
	warn "WARNING: $opt{OUTPUT_DIR}/$runName.filtered_indels.vcf already exists, skipping \n";
	return $jobID;
    } elsif ($opt{FILTER_MODE} eq "BOTH" && -e "$opt{OUTPUT_DIR}/$runName.filtered_variants.vcf"){
	warn "WARNING: $opt{OUTPUT_DIR}/$runName.filtered_variants.vcf already exists, skipping \n";
	return $jobID;
    }

    ### Build Queue command
    my $javaMem = $opt{FILTER_THREADS} * $opt{FILTER_MEM};
    my $command = "java -Xmx".$javaMem."G -Xms".$opt{FILTER_MEM}."G -jar $opt{QUEUE_PATH}/Queue.jar "; ### Change memory allocation here!!!!!
    $command .= "-jobQueue $opt{FILTER_QUEUE} -jobEnv \"threaded $opt{FILTER_THREADS}\" -jobRunner GridEngine -jobReport $opt{OUTPUT_DIR}/logs/filterVariants.jobReport.txt "; #Queue options

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
    print FILTER_SH "$command\n\n";

    if ($opt{FILTER_MODE} eq "SNP"){
	print FILTER_SH "if [ -f $opt{OUTPUT_DIR}/tmp/.$runName\.filtered_snps.vcf.done ]\n";
	print FILTER_SH "then\n";
	print FILTER_SH "\tmv $opt{OUTPUT_DIR}/tmp/$runName\.filtered_snps.vcf $opt{OUTPUT_DIR}/\n";
	print FILTER_SH "\tmv $opt{OUTPUT_DIR}/tmp/$runName\.filtered_snps.vcf.idx $opt{OUTPUT_DIR}/\n";
	print FILTER_SH "fi\n\n";
    } elsif ($opt{FILTER_MODE} eq "INDEL"){
	print FILTER_SH "if [ -f $opt{OUTPUT_DIR}/tmp/.$runName\.filtered_indels.vcf.done ]\n";
	print FILTER_SH "then\n";
	print FILTER_SH "\tmv $opt{OUTPUT_DIR}/tmp/$runName\.filtered_indels.vcf $opt{OUTPUT_DIR}/\n";
	print FILTER_SH "\tmv $opt{OUTPUT_DIR}/tmp/$runName\.filtered_indels.vcf.idx $opt{OUTPUT_DIR}/\n";
	print FILTER_SH "fi\n\n";
    } elsif ($opt{FILTER_MODE} eq "BOTH"){
	print FILTER_SH "if [ -f $opt{OUTPUT_DIR}/tmp/.$runName\.filtered_variants.vcf.done ]\n";
	print FILTER_SH "then\n";
	print FILTER_SH "\tmv $opt{OUTPUT_DIR}/tmp/$runName\.filtered_variants.vcf $opt{OUTPUT_DIR}/\n";
	print FILTER_SH "\tmv $opt{OUTPUT_DIR}/tmp/$runName\.filtered_variants.vcf.idx $opt{OUTPUT_DIR}/\n";
	print FILTER_SH "fi\n\n";
    }
    print FILTER_SH "touch filterVariants.done \n";
    
    ### Process runningjobs
    foreach my $sample (keys %{$opt{RUNNING_JOBS}}){
	push(@runningJobs, join(",",@{$opt{RUNNING_JOBS}->{$sample}}));
    }

    ### Start main bash script
    if (@runningJobs){
	system "qsub -q $opt{FILTER_QUEUE} -pe threaded $opt{FILTER_THREADS} -o $logDir -e $logDir -N $jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
    } else {
	system "qsub -q $opt{FILTER_QUEUE} -pe threaded $opt{FILTER_THREADS} -o $logDir -e $logDir -N $jobID $bashFile";
    }

    return $jobID;
}

sub readConfiguration{
    my $configuration = shift;
    my %opt;

    foreach my $key (keys %{$configuration}){
	$opt{$key} = $configuration->{$key};
    }

    if(! $opt{QUEUE_PATH}){ die "ERROR: No PICARD_PATH found in .conf file\n" }
    if(! $opt{FILTER_QUEUE}){ die "ERROR: No FILTER_QUEUE found in .conf file\n" }
    if(! $opt{FILTER_THREADS}){ die "ERROR: No FILTER_THREADS found in .conf file\n" }
    if(! $opt{FILTER_MEM}){ die "ERROR: No FILTER_QUEUE found in .conf file\n" }
    if(! $opt{FILTER_SCATTER}){ die "ERROR: No FILTER_SCATTER found in .conf file\n" }
    if(! $opt{FILTER_SCALA}){ die "ERROR: No FILTER_SCALA found in .conf file\n" }
    if(! $opt{FILTER_MODE}){ die "ERROR: No FILTER_MODE  found in .conf file\n" }
    if($opt{FILTER_MODE} ne "SNP" and $opt{FILTER_MODE} ne "INDEL" and $opt{FILTER_MODE} ne "BOTH"){ die "ERROR: FILTER_MODE $opt{FILTER_MODE} does not exist use SNP, INDEL or BOTH\n"}
    if ($opt{FILTER_MODE} eq "SNP" || $opt{FILTER_MODE} eq "BOTH") {
	if(! $opt{FILTER_SNPNAME}){ die "ERROR: No FILTER_SNPNAME found in .conf file\n" }
	if(! $opt{FILTER_SNPEXPR}){ die "ERROR: No FILTER_SNPEXPR  found in .conf file\n" }
    }
    if ($opt{FILTER_MODE} eq "INDEL" || $opt{FILTER_MODE} eq "BOTH") {
	if(! $opt{FILTER_INDELNAME}){ die "ERROR: No FILTER_INDELNAME found in .conf file\n" }
	if(! $opt{FILTER_INDELEXPR}){ die "ERROR: No FILTER_INDELEXPR found in .conf file\n" }
    }
    if(! $opt{GENOME}){ die "ERROR: No GENOME found in .conf file\n" }
    if(! $opt{OUTPUT_DIR}){ die "ERROR: No OUTPUT_DIR found in .conf file\n" }

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