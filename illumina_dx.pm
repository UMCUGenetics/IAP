#!/usr/bin/perl -w

##################################################################################################################################################
###This script is designed to run Diagnostics (Dx) specific exome sequencing pipeline steps.
###
###Author: R.F.Ernst
##################################################################################################################################################

package illumina_dx;

use strict;
use POSIX qw(tmpnam);

sub runDX {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my @runningJobs;
    my $jobID = "DX_".get_job_id();

    if (-e "$opt{OUTPUT_DIR}/logs/DX.done"){
	warn "WARNING: $opt{OUTPUT_DIR}/logs/DX.done exists, skipping \n";
	return $jobID;
    }
    
    ### vcf file
    my $vcf;
    if($opt{FILTER_VARIANTS} eq "yes"){
	if($opt{FILTER_MODE} eq "BOTH"){ $vcf = $runName.".filtered_variants.vcf";}
	if($opt{FILTER_MODE} eq "SNP"){ $vcf = $runName.".filtered_snps.vcf";}
	if($opt{FILTER_MODE} eq "INDEL"){ $vcf = $runName.".filtered_indels.vcf";}
    } elsif ($opt{FILTER_VARIANTS} eq "no"){ 
	$vcf = $runName.".raw_variants.vcf";
    }
    
    ### Create main bash script
    my $bashFile = $opt{OUTPUT_DIR}."/jobs/".$jobID.".sh";
    my $logDir = $opt{OUTPUT_DIR}."/logs";
    
    open DX_SH, ">$bashFile" or die "cannot open file $bashFile \n";
    print DX_SH "#!/bin/bash\n\n";
    print DX_SH "bash $opt{CLUSTER_PATH}/settings.sh\n\n";
    print DX_SH "cd $opt{OUTPUT_DIR}/\n\n";
    print DX_SH "echo \"Start DX\t\" `date` \"\t$vcf\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n\n";
    
    ### Run kinship analyses
    if ( $opt{DX_KINSHIP} eq "yes" ) {
	if (-e "$opt{OUTPUT_DIR}/logs/Kinship.done"){
	    warn "WARNING: $opt{OUTPUT_DIR}/logs/Kinship.done exists, skipping \n";
	} else {
	    print DX_SH "cd $opt{OUTPUT_DIR}/tmp/\n";
	    print DX_SH "$opt{VCFTOOLS_PATH}/vcftools --vcf $opt{OUTPUT_DIR}/$vcf --plink\n";
	    print DX_SH "$opt{PLINK_PATH}/plink --file out --make-bed --noweb\n";
	    print DX_SH "$opt{PLINK_PATH}/king -b plink.bed --kinship\n\n";
	    print DX_SH "cp king.kin0 $opt{OUTPUT_DIR}/$runName.kinship\n";
	    print DX_SH "mv $opt{OUTPUT_DIR}/tmp/plink.log $opt{OUTPUT_DIR}/logs/\n";
	    print DX_SH "mv $opt{OUTPUT_DIR}/tmp/out.log $opt{OUTPUT_DIR}/logs/\n";
	    print DX_SH "touch $opt{OUTPUT_DIR}/logs/Kinship.done\n\n";
	}
    }

    ### Phase by transmission
    if ( $opt{DX_PHASE} eq "yes" ) {
	if (-e "$opt{OUTPUT_DIR}/logs/Phase.done"){
	    warn "WARNING: $opt{OUTPUT_DIR}/logs/Phase.done exists, skipping \n";
	} else {
	    print DX_SH "cd $opt{OUTPUT_DIR}/tmp/\n";
	    print DX_SH "java -Xmx8G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T PhaseByTransmission -R $opt{GENOME} -V $opt{OUTPUT_DIR}/$vcf -ped $opt{DX_PED} -o $runName.phased.vcf.gz --MendelianViolationsFile $runName.MendelViol\n";
	    print DX_SH "mv $runName.phased.vcf.gz $opt{OUTPUT_DIR}/\n";
	    print DX_SH "mv $runName.MendelViol $opt{OUTPUT_DIR}/\n";
	    print DX_SH "touch $opt{OUTPUT_DIR}/logs/Phase.done\n\n";
	}
    }

    print DX_SH "touch $opt{OUTPUT_DIR}/logs/DX.done\n";
    print DX_SH "echo \"End DX\t\" `date` \"\t$vcf\t\" `uname -n` >> $opt{OUTPUT_DIR}/logs/$runName.log\n";

    ### Process runningjobs
    foreach my $sample (@{$opt{SAMPLES}}){
	if( exists $opt{RUNNING_JOBS}->{$sample} && @{$opt{RUNNING_JOBS}->{$sample}} ) {
	    push(@runningJobs, join(",",@{$opt{RUNNING_JOBS}->{$sample}}));
	}
    }
    
    ### Run job
    if (@runningJobs){
	system "qsub -q $opt{DX_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{DX_THREADS} -o $logDir/DX_$runName.out -e $logDir/DX_$runName.err -N $jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
    } else {
	system "qsub -q $opt{DX_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{DX_THREADS} -o $logDir/DX_$runName.out -e $logDir/DX_$runName.err -N $jobID $bashFile";
    }
    
    return $jobID;
}

sub readConfiguration{
    my $configuration = shift;
    my %opt;

    foreach my $key (keys %{$configuration}){
	$opt{$key} = $configuration->{$key};
    }

    if(! $opt{DX_QUEUE}){ die "ERROR: No DX_QUEUE found in .ini file\n" }
    if(! $opt{DX_THREADS}){ die "ERROR: No DX_THREADS found in .ini file\n" }
    #if(! $opt{DX_SCATTER}){ die "ERROR: No DX_SCATTER found in .ini file\n" }
    if(! $opt{DX_MEM}){ die "ERROR: No DX_MEM found in .ini file\n" }
    if(! $opt{DX_KINSHIP}){ die "ERROR: No DX_KINSHIP found in .ini file\n" }
    if(! $opt{PLINK_PATH}){ die "ERROR: No PLINK_PATH found in .ini file\n" }
    if(! $opt{VCFTOOLS_PATH}){ die "ERROR: No VCFTOOLS_PATH found in .ini file\n" }
    if(! $opt{DX_PHASE}){ die "ERROR: No DX_PHASE found in .ini file\n" }
    if(! $opt{DX_PED}){ die "ERROR: No PED found in .conf file\n" }
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