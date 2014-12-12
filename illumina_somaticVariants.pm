#!/usr/bin/perl -w

##################################################################################################################################################
###This script is designed to run Somatic Variant callers and combine the result.
###
###
###Author: R.F.Ernst
###Latest change:
###TODO: CPCT sample name parsing, add callers, combine the result
##################################################################################################################################################

package illumina_somaticVariants;

use strict;
use POSIX qw(tmpnam);
use File::Path qw(make_path);

### Parse sample names
# Expects CPCT samples (CPCT........T/R)
# Creates directory structure
sub parseSamples {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
    my %somatic_samples;

    foreach my $sample (@{$opt{SAMPLES}}){
	# Parse cpct samples based on expected naming
	my ($cpct_name,$origin) = ($sample =~ /(CPCT\d{8})([TR].*)/);

	# Reference sample
	if ($origin =~ m/R.*/){
	    if ($somatic_samples{$cpct_name}{"ref"}){
		warn "\t WARNING: $cpct_name has multiple reference samples, using: $somatic_samples{$cpct_name}{'ref'} \n";
	    } else {
		$somatic_samples{$cpct_name}{"ref"} = $sample;
	    }
	}

	# Tumor samples
	elsif ($origin =~ m/T.*/){
	    push(@{$somatic_samples{$cpct_name}{"tumor"}},$sample);
	}
    }

    $opt{SOMATIC_SAMPLES} = {%somatic_samples};
    return \%opt;
}

### Somatic Variant Callers
sub runSomaticVariantCallers {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
    my @somvar_jobs;
    
    if($opt{SOMVAR_STRELKA} eq "yes"){
	print "\n###SCHEDULING STRELKA####\n";
	my @strelka_jobs = @{runStrelka(\%opt)};
	push(@somvar_jobs, @strelka_jobs);
    }
    if($opt{SOMVAR_VARSCAN} eq "yes"){
	print "\n###SCHEDULING VARSCAN####\n";
	my @varscan_jobs = @{runVarscan(\%opt)};
	push(@somvar_jobs, @varscan_jobs);
    }
    if($opt{SOMVAR_FREEBAYES} eq "yes"){
	print "\n###SCHEDULING FREEBAYES####\n";
	my @freebayes_jobs = @{runFreeBayes(\%opt)};
	push(@somvar_jobs, @freebayes_jobs);
    }
    print "\n\n";
    print @somvar_jobs;
    print "\n\n";
}


sub runStrelka {
    my %opt = %{$_[0]};
    my @strelka_jobs;
    my $job_dir = "$opt{OUTPUT_DIR}/somaticVariants/strelka/jobs";
    my $log_dir = "$opt{OUTPUT_DIR}/somaticVariants/strelka/logs";

    ### Create strelka output, job and log dirs
    if(! -e "$opt{OUTPUT_DIR}/somaticVariants/strelka"){
	make_path("$opt{OUTPUT_DIR}/somaticVariants/strelka") or die "Couldn't create directory: $opt{OUTPUT_DIR}/somaticVariants/strelka\n";
    }
    if(! -e $job_dir){
	make_path($job_dir) or die "Couldn't create directory: $job_dir\n";
    }
    if(! -e $log_dir){
	make_path($log_dir) or die "Couldn't create directory: $log_dir\n";
    }
    
    ### Run Strelka per tumor sample
    foreach my $sample (keys(%{$opt{SOMATIC_SAMPLES}})){
	foreach my $sample_tumor (@{$opt{SOMATIC_SAMPLES}{$sample}{'tumor'}}){

	    ## Lookup bams and running jobs
	    my $sample_tumor_bam = "$opt{OUTPUT_DIR}/$sample_tumor/mapping/$opt{BAM_FILES}->{$sample_tumor}";
	    my @running_jobs;
	    if ( @{$opt{RUNNING_JOBS}->{$sample_tumor}} ){
		push(@running_jobs, @{$opt{RUNNING_JOBS}->{$sample_tumor}});
	    }

	    my $sample_ref_bam = "$opt{OUTPUT_DIR}/$opt{SOMATIC_SAMPLES}{$sample}{'ref'}/mapping/$opt{BAM_FILES}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}";
	    if ( @{$opt{RUNNING_JOBS}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}} ){
		push(@running_jobs, @{$opt{RUNNING_JOBS}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}});
	    }

	    # Print sample and bam info
	    print "$sample \t $sample_ref_bam \t $sample_tumor_bam \n";

	    ### Skip Strelka if .done file exist
	    if (-e "$log_dir/strelka.done"){
		warn "WARNING: $log_dir/$sample_tumor.done, skipping \n";
		next;
	    }

	    ## Create strelka bash script
	    my $job_id = "STR_".$sample_tumor."_".get_job_id();
	    my $bash_file = $job_dir."/".$job_id.".sh";
	    my $sample_tumor_ref_out = "$opt{OUTPUT_DIR}/somaticVariants/strelka/$opt{SOMATIC_SAMPLES}{$sample}{'ref'}\_$sample_tumor";

	    open STRELKA_SH, ">$bash_file" or die "cannot open file $bash_file \n";
	    print STRELKA_SH "#!/bin/bash\n\n";
	    print STRELKA_SH "echo \"Start Strelka\t\" `date` \"\t$sample \t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/$sample_tumor.log\n\n";

	    # Run Strelka
	    print STRELKA_SH "$opt{STRELKA_PATH}/bin/configureStrelkaWorkflow.pl --tumor $sample_tumor_bam --normal $sample_ref_bam --ref $opt{GENOME} --config $opt{STRELKA_INI} --output-dir $sample_tumor_ref_out\n\n";

	    print STRELKA_SH "cd $sample_tumor_ref_out\n";
	    print STRELKA_SH "make -j 8\n\n";

	    # Check strelka completed
	    print STRELKA_SH "if [ -f $sample_tumor_ref_out/task.complete ]\n";
	    print STRELKA_SH "then\n";
	    print STRELKA_SH "\ttouch $log_dir/$sample_tumor.done\n";
	    print STRELKA_SH "fi\n\n";
	    print STRELKA_SH "echo \"End Strelka\t\" `date` \"\t$sample \t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/$sample_tumor.log\n\n";
	    close STRELKA_SH;

	    ## Run job
	    if ( @running_jobs ){
		system "qsub -q $opt{STRELKA_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{STRELKA_THREADS} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
		push(@strelka_jobs, $job_id);
	    } else {
		system "qsub -q $opt{STRELKA_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{STRELKA_THREADS} -o $log_dir -e $log_dir -N $job_id $bash_file";
		push(@strelka_jobs, $job_id);
	    }
	}
    }
    return \@strelka_jobs;
}

sub runVarscan {
    my %opt = %{$_[0]};
    my @varscan_jobs;
    my $job_dir = "$opt{OUTPUT_DIR}/somaticVariants/varscan/jobs";
    my $log_dir = "$opt{OUTPUT_DIR}/somaticVariants/varscan/logs";

    ### Create varscan output, job and log dirs
    if(! -e $job_dir){
	make_path($job_dir) or die "Couldn't create directory: $job_dir\n";
    }
    if(! -e $log_dir){
	make_path($log_dir) or die "Couldn't create directory: $log_dir\n";
    }
    
    ### Run Varscan for tumor sample
    foreach my $sample (keys(%{$opt{SOMATIC_SAMPLES}})){
	foreach my $sample_tumor (@{$opt{SOMATIC_SAMPLES}{$sample}{'tumor'}}){

	    ## Lookup bams and running jobs
	    my $sample_tumor_bam = "$opt{OUTPUT_DIR}/$sample_tumor/mapping/$opt{BAM_FILES}->{$sample_tumor}";
	    my @running_jobs;
	    if ( @{$opt{RUNNING_JOBS}->{$sample_tumor}} ){
		push(@running_jobs, @{$opt{RUNNING_JOBS}->{$sample_tumor}});
	    }

	    my $sample_ref_bam = "$opt{OUTPUT_DIR}/$opt{SOMATIC_SAMPLES}{$sample}{'ref'}/mapping/$opt{BAM_FILES}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}";
	    if ( @{$opt{RUNNING_JOBS}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}} ){
		push(@running_jobs, @{$opt{RUNNING_JOBS}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}});
	    }

	    #  Print sample and bam info
	    print "$sample \t $sample_ref_bam \t $sample_tumor_bam \n";

	    ### Skip VarScan if .done file exist
	    if (-e "$log_dir/varscan.done"){
		warn "WARNING: $log_dir/varscan.done, skipping \n";
		next;
	    }

	    ## Setup varscan bash script and output dir
	    my $job_id = "VS_".$sample_tumor."_".get_job_id();
	    my $bash_file = $job_dir."/".$job_id.".sh";

	    # Create output dir
	    my $output_name = "$opt{SOMATIC_SAMPLES}{$sample}{'ref'}\_$sample_tumor";
	    my $sample_tumor_ref_out = "$opt{OUTPUT_DIR}/somaticVariants/varscan/$output_name";
	    if(! -e $sample_tumor_ref_out){
		make_path($sample_tumor_ref_out) or die "Couldn't create directory: $sample_tumor_ref_out\n";
	    }

	    # Create bash script
	    open VARSCAN_SH, ">$bash_file" or die "cannot open file $bash_file \n";
	    print VARSCAN_SH "#!/bin/bash\n\n";
	    print VARSCAN_SH "cd $sample_tumor_ref_out\n";
	    print VARSCAN_SH "echo \"Start Pileup\t\" `date` \"\t$sample \t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/$sample_tumor.log\n";
	    # run pileups
	    if (!-e "$sample_tumor_bam.pileup") {
		print VARSCAN_SH "samtools mpileup -q 1 -f $opt{GENOME} $sample_tumor_bam > $sample_tumor_bam.pileup\n";
	    }
	    if (!-e "$sample_ref_bam.pileup") {    
		print VARSCAN_SH "samtools mpileup -q 1 -f $opt{GENOME} $sample_ref_bam > $sample_ref_bam.pileup\n";
	    }
	    print VARSCAN_SH "echo \"End Pileup\t\" `date` \"\t$sample \t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/$sample_tumor.log\n\n";

	    # run varscan
	    print VARSCAN_SH "echo \"Start Varscan\t\" `date` \"\t$sample \t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/$sample_tumor.log\n";
	    print VARSCAN_SH "java -Xmx12g -jar $opt{VARSCAN_PATH} somatic $sample_ref_bam.pileup $sample_tumor_bam.pileup $output_name $opt{VARSCAN_SETTINGS} --output-vcf 1\n\n";

	    # postprocessing
	    print VARSCAN_SH "java -Xmx12g -jar $opt{VARSCAN_PATH} processSomatic $output_name.indel.vcf $opt{VARSCAN_POSTSETTINGS}\n";
	    print VARSCAN_SH "java -Xmx12g -jar $opt{VARSCAN_PATH} processSomatic $output_name.snp.vcf $opt{VARSCAN_POSTSETTINGS}\n\n";
	    
	    # Check varscan completed
	    print VARSCAN_SH "touch $log_dir/$sample_tumor.done\n\n"; ## Check on complete output!!!
	    print VARSCAN_SH "echo \"End Varscan\t\" `date` \"\t$sample \t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/$sample_tumor.log\n";
	    
	    close VARSCAN_SH;

	    # Run job
	    if ( @running_jobs ){
		system "qsub -q $opt{VARSCAN_QUEUE} -pe threaded $opt{VARSCAN_THREADS} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
		push(@varscan_jobs, $job_id);
	    } else {
		system "qsub -q $opt{VARSCAN_QUEUE} -pe threaded $opt{VARSCAN_THREADS} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id $bash_file";
		push(@varscan_jobs, $job_id);
	    }
	}
    }
    return \@varscan_jobs;
}

sub runFreeBayes {
    my %opt = %{$_[0]};
    my @freebayes_jobs;
    my $job_dir = "$opt{OUTPUT_DIR}/somaticVariants/freebayes/jobs";
    my $log_dir = "$opt{OUTPUT_DIR}/somaticVariants/freebayes/logs";
    
    ### Create freebayes output, job and log dirs
    if(! -e $job_dir){
	make_path($job_dir) or die "Couldn't create directory: $job_dir\n";
    }
    if(! -e $log_dir){
	make_path($log_dir) or die "Couldn't create directory: $log_dir\n";
    }

    foreach my $sample (keys(%{$opt{SOMATIC_SAMPLES}})){
	foreach my $sample_tumor (@{$opt{SOMATIC_SAMPLES}{$sample}{'tumor'}}){

	    # Lookup bams and running jobs
	    my $sample_tumor_bam = "$opt{OUTPUT_DIR}/$sample_tumor/mapping/$opt{BAM_FILES}->{$sample_tumor}";
	    my @running_jobs;
	    if ( @{$opt{RUNNING_JOBS}->{$sample_tumor}} ){
		push(@running_jobs, @{$opt{RUNNING_JOBS}->{$sample_tumor}});
	    }
	    my $sample_ref_bam = "$opt{OUTPUT_DIR}/$opt{SOMATIC_SAMPLES}{$sample}{'ref'}/mapping/$opt{BAM_FILES}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}";
	    if ( @{$opt{RUNNING_JOBS}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}} ){
		push(@running_jobs, @{$opt{RUNNING_JOBS}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}});
	    }

	    #Print sample and bam info
	    print "$sample \t $sample_ref_bam \t $sample_tumor_bam \n";

	    ### Skip FreeBayes if .done file exist
	    if (-e "$log_dir/freebayes.done"){
		warn "WARNING: $log_dir/freebayes.done, skipping \n";
		next;
	    }

	    ## Setup varscan bash script and output dir
	    my $job_id = "FB_".$sample_tumor."_".get_job_id();
	    my $bash_file = $job_dir."/".$job_id.".sh";

	    # Create output dir
	    my $output_name = "$opt{SOMATIC_SAMPLES}{$sample}{'ref'}\_$sample_tumor";
	    my $sample_tumor_ref_out = "$opt{OUTPUT_DIR}/somaticVariants/freebayes/$output_name";
	    if(! -e $sample_tumor_ref_out){
		make_path($sample_tumor_ref_out) or die "Couldn't create directory: $sample_tumor_ref_out\n";
	    }

	    # Create bash script
	    open FREEBAYES_SH, ">$bash_file" or die "cannot open file $bash_file \n";
	    print FREEBAYES_SH "#!/bin/bash\n\n";
	    print FREEBAYES_SH "echo \"Start Freebayes\t\" `date` \"\t$sample \t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/$sample_tumor.log\n";
	    
	    # Run freebayes
	    print FREEBAYES_SH "$opt{FREEBAYES_PATH}/freebayes -f $opt{GENOME} -t $opt{FREEBAYES_TARGETS} $opt{FREEBAYES_SETTINGS} $sample_ref_bam $sample_tumor_bam > $sample_tumor_ref_out/$output_name.vcf\n\n";

	    # Uniqify freebayes output 
	    print FREEBAYES_SH "uniq $sample_tumor_ref_out/$output_name.vcf > $sample_tumor_ref_out/$output_name.uniq.vcf\n";
	    print FREEBAYES_SH "mv $sample_tumor_ref_out/$output_name.uniq.vcf > $sample_tumor_ref_out/$output_name.vcf\n\n";

	    # get sample ids
	    print FREEBAYES_SH "sample_R=`grep -P \"^#CHROM\" $sample_tumor_ref_out/$output_name.vcf | cut -f 10`\n";
	    print FREEBAYES_SH "sample_T=`grep -P \"^#CHROM\" $sample_tumor_ref_out/$output_name.vcf | cut -f 11`\n\n";

	    # annotate somatic and germline scores
	    print FREEBAYES_SH "/hpc/local/CentOS6/cog_bioinf/vcflib/bin/vcfsamplediff VT \$sample_R \$sample_T $sample_tumor_ref_out/$output_name.vcf> $sample_tumor_ref_out/$output_name\_VTannot.vcf\n";
	    print FREEBAYES_SH "grep -P \"^#\" $sample_tumor_ref_out/$output_name\_VTannot.vcf > $sample_tumor_ref_out/$output_name\_germline.vcf\n";
	    print FREEBAYES_SH "grep -P \"^#\" $sample_tumor_ref_out/$output_name\_VTannot.vcf > $sample_tumor_ref_out/$output_name\_somatic.vcf\n";
	    print FREEBAYES_SH "grep -i \"VT=germline\" $sample_tumor_ref_out/$output_name\_VTannot.vcf >> $sample_tumor_ref_out/$output_name\_germline.vcf\n";
	    print FREEBAYES_SH "grep -i \"VT=somatic\" $sample_tumor_ref_out/$output_name\_VTannot.vcf >> $sample_tumor_ref_out/$output_name\_somatic.vcf\n";
	    print FREEBAYES_SH "rm $sample_tumor_ref_out/$output_name\_VTannot.vcf\n\n";
	    
	    print FREEBAYES_SH "touch $log_dir/$sample_tumor.done\n\n"; ## Check on complete output!!!
	    print FREEBAYES_SH "echo \"End Freebayes\t\" `date` \"\t$sample \t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/$sample_tumor.log\n";

	    close FREEBAYES_SH;

	    # Run job
	    if ( @running_jobs ){
		system "qsub -q $opt{FREEBAYES_QUEUE} -pe threaded $opt{FREEBAYES_THREADS} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
		push(@freebayes_jobs, $job_id);
	    } else {
		system "qsub -q $opt{FREEBAYES_QUEUE} -pe threaded $opt{FREEBAYES_THREADS} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id $bash_file";
		push(@freebayes_jobs, $job_id);
	    }
	}
    }
    return \@freebayes_jobs;
}

### Merge results
sub mergeSomaticVariants {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
}

### Read config and check input
sub readConfiguration{
    my $configuration = shift;
    my %opt;

    foreach my $key (keys %{$configuration}){
	$opt{$key} = $configuration->{$key};
    }


    if(! $opt{GENOME}){ die "ERROR: No GENOME found in .ini file\n" }
    elsif(! -e $opt{GENOME}){ die"ERROR: $opt{GENOME} does not exist\n"}
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
