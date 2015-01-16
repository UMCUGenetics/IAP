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

### Run and merge
sub runSomaticVariantCallers {
    my $configuration = shift;
    my %opt = %{readConfiguration($configuration)};
    my @all_somvar_jobs;
    ### Loop over tumor samples
    foreach my $sample (keys(%{$opt{SOMATIC_SAMPLES}})){
	foreach my $sample_tumor (@{$opt{SOMATIC_SAMPLES}{$sample}{'tumor'}}){
	    my @somvar_jobs;
	    ## Create output, log and job directories
	    my $sample_tumor_name = "$opt{SOMATIC_SAMPLES}{$sample}{'ref'}\_$sample_tumor";
	    my $sample_tumor_out_dir = "$opt{OUTPUT_DIR}/somaticVariants/$sample_tumor_name";
	    my $sample_tumor_log_dir = "$sample_tumor_out_dir/logs/";
	    my $sample_tumor_job_dir = "$sample_tumor_out_dir/jobs/";

	    if(! -e $sample_tumor_out_dir){
		make_path($sample_tumor_out_dir) or die "Couldn't create directory:  $sample_tumor_out_dir\n";
	    }
	    if(! -e $sample_tumor_job_dir){
		make_path($sample_tumor_job_dir) or die "Couldn't create directory: $sample_tumor_job_dir\n";
	    }
	    if(! -e $sample_tumor_log_dir){
		make_path($sample_tumor_log_dir) or die "Couldn't create directory: $sample_tumor_log_dir\n";
	    }

	    ## Lookup running jobs and bams
	    my $sample_tumor_bam = "$opt{OUTPUT_DIR}/$sample_tumor/mapping/$opt{BAM_FILES}->{$sample_tumor}";
	    my @running_jobs;
	    if ( @{$opt{RUNNING_JOBS}->{$sample_tumor}} ){
		push(@running_jobs, @{$opt{RUNNING_JOBS}->{$sample_tumor}});
	    }
	    my $sample_ref_bam = "$opt{OUTPUT_DIR}/$opt{SOMATIC_SAMPLES}{$sample}{'ref'}/mapping/$opt{BAM_FILES}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}";
	    if ( @{$opt{RUNNING_JOBS}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}} ){
		push(@running_jobs, @{$opt{RUNNING_JOBS}->{$opt{SOMATIC_SAMPLES}{$sample}{'ref'}}});
	    }

	    ## Print sample and bam info
	    print "$sample \t $sample_ref_bam \t $sample_tumor_bam \n";

	    ## Skip Somatic Callers if .done file exist
	    if (-e "$sample_tumor_log_dir/$sample_tumor_name.done"){
		warn "WARNING: $sample_tumor_log_dir/$sample_tumor_name.done, skipping \n";
		next;
	    }

	    ## Run somatic callers
	    if($opt{SOMVAR_STRELKA} eq "yes"){
		print "\n###SCHEDULING STRELKA####\n";
		my $strelka_job = runStrelka($sample_tumor, $sample_tumor_out_dir, $sample_tumor_job_dir, $sample_tumor_log_dir, $sample_tumor_bam, $sample_ref_bam, \@running_jobs, \%opt);
		if($strelka_job){push(@somvar_jobs, $strelka_job)};
	    }
	    if($opt{SOMVAR_VARSCAN} eq "yes"){
		print "\n###SCHEDULING VARSCAN####\n";
		my $varscan_job = runVarscan($sample_tumor, $sample_tumor_name, $sample_tumor_out_dir, $sample_tumor_job_dir, $sample_tumor_log_dir, $sample_tumor_bam, $sample_ref_bam, \@running_jobs, \%opt);
		if($varscan_job){push(@somvar_jobs, $varscan_job)};
	    }
	    if($opt{SOMVAR_FREEBAYES} eq "yes"){
		print "\n###SCHEDULING FREEBAYES####\n";
		my $freebayes_job = runFreeBayes($sample_tumor, $sample_tumor_name, $sample_tumor_out_dir, $sample_tumor_job_dir, $sample_tumor_log_dir, $sample_tumor_bam, $sample_ref_bam, \@running_jobs, \%opt);
		if($freebayes_job){push(@somvar_jobs, $freebayes_job)};
	    }

	    ## Merge somatic vcfs
	    print "\n###SCHEDULING MERGE SOMATIC VCFS####\n";

	    my $job_id = "MERGE_".$sample_tumor."_".get_job_id();
	    my $bash_file = $sample_tumor_job_dir."/".$job_id.".sh";

	    open MERGE_SH, ">$bash_file" or die "cannot open file $bash_file \n";
	    print MERGE_SH "#!/bin/bash\n\n";
	    print MERGE_SH "echo \"Start Merge\t\" `date` `uname -n` >> $sample_tumor_log_dir/merge.log\n\n";

	    print MERGE_SH "java -Xmx6G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T CombineVariants -R $opt{GENOME} -o $sample_tumor_out_dir/$sample_tumor_name\_merged_somatics.vcf --genotypemergeoption uniquify --assumeIdenticalSamples ";
	    if($opt{SOMVAR_STRELKA} eq "yes"){ print MERGE_SH "-V $sample_tumor_out_dir/strelka/passed.somatic.merged.vcf "; }
	    if($opt{SOMVAR_VARSCAN} eq "yes"){ print MERGE_SH "-V $sample_tumor_out_dir/varscan/$sample_tumor_name.merged.Somatic.hc.vcf "; }
	    if($opt{SOMVAR_FREEBAYES} eq "yes"){ print MERGE_SH "-V $sample_tumor_out_dir/freebayes/$sample_tumor_name\_somatic_filtered.vcf "; }

	    print MERGE_SH "\n\ntouch $sample_tumor_log_dir/$sample_tumor_name.done\n\n";

	    print MERGE_SH "echo \"END Merge\t\" `date` `uname -n` >> $sample_tumor_log_dir/merge.log\n\n";
	    close MERGE_SH;

	    # Run job
	    if ( @somvar_jobs ){
		system "qsub -q $opt{SOMVARMERGE_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{SOMVARMERGE_THREADS} -o $sample_tumor_log_dir -e $sample_tumor_log_dir -N $job_id -hold_jid ".join(",",@somvar_jobs)." $bash_file";
	    } else {
		system "qsub -q $opt{SOMVARMERGE_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{SOMVARMERGE_THREADS} -o $sample_tumor_log_dir -e $sample_tumor_log_dir -N $job_id $bash_file";
	    }
	    push(@all_somvar_jobs, @somvar_jobs);
	}
    }
    return \@all_somvar_jobs;
}

### Somatic Variant Callers
sub runStrelka {
    my ($sample_tumor, $out_dir, $job_dir, $log_dir, $sample_tumor_bam, $sample_ref_bam, $running_jobs, $opt) = (@_);
    my @running_jobs = @{$running_jobs};
    my %opt = %{$opt};
    my $strelka_out_dir = "$out_dir/strelka";

    ## Skip Strelka if .done file exist
    if (-e "$log_dir/strelka.done"){
	warn "WARNING: $log_dir/strelka.done, skipping \n";
	return;
    }

    ## Create strelka bash script
    my $job_id = "STR_".$sample_tumor."_".get_job_id();
    my $bash_file = $job_dir."/".$job_id.".sh";

    open STRELKA_SH, ">$bash_file" or die "cannot open file $bash_file \n";
    print STRELKA_SH "#!/bin/bash\n\n";
    print STRELKA_SH "echo \"Start Strelka\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/strelka.log\n\n";

    # Run Strelka
    print STRELKA_SH "$opt{STRELKA_PATH}/bin/configureStrelkaWorkflow.pl --tumor $sample_tumor_bam --normal $sample_ref_bam --ref $opt{GENOME} --config $opt{STRELKA_INI} --output-dir $strelka_out_dir\n\n";

    print STRELKA_SH "cd $strelka_out_dir\n";
    print STRELKA_SH "make -j 8\n\n";

    # Check strelka completed
    print STRELKA_SH "if [ -f $strelka_out_dir/task.complete ]\n";
    print STRELKA_SH "then\n";
    print STRELKA_SH "\tjava -Xmx6G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T CombineVariants -R $opt{GENOME} -o passed.somatic.merged.vcf -V results/passed.somatic.snvs.vcf -V results/passed.somatic.indels.vcf \n\n";
    print STRELKA_SH "\ttouch $log_dir/strelka.done\n";
    print STRELKA_SH "fi\n\n";
    print STRELKA_SH "echo \"End Strelka\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/strelka.log\n\n";
    close STRELKA_SH;

    ## Run job
    if ( @running_jobs ){
	system "qsub -q $opt{STRELKA_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{STRELKA_THREADS} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
    } else {
	system "qsub -q $opt{STRELKA_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{STRELKA_THREADS} -o $log_dir -e $log_dir -N $job_id $bash_file";
    }

    return $job_id;
}

sub runVarscan {
    my ($sample_tumor, $sample_tumor_name, $out_dir, $job_dir, $log_dir, $sample_tumor_bam, $sample_ref_bam, $running_jobs, $opt) = (@_);
    my @running_jobs = @{$running_jobs};
    my %opt = %{$opt};
    my $varscan_out_dir = "$out_dir/varscan";

    ## Create output dir
    if( ! -e $varscan_out_dir ){
	make_path($varscan_out_dir) or die "Couldn't create directory: $varscan_out_dir\n";
    }

    ## Skip varscan if .done file exist
    if ( -e "$log_dir/varscan.done" ){
	warn "WARNING: $log_dir/varscan.done, skipping \n";
	return;
    }

    ## Setup varscan bash script and output dir
    my $job_id = "VS_".$sample_tumor."_".get_job_id();
    my $bash_file = $job_dir."/".$job_id.".sh";

    ## Create bash script
    open VARSCAN_SH, ">$bash_file" or die "cannot open file $bash_file \n";
    print VARSCAN_SH "#!/bin/bash\n\n";
    print VARSCAN_SH "cd $varscan_out_dir\n";
    print VARSCAN_SH "echo \"Start Pileup\t\" `date` \"$sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/varscan.log\n";

    # run pileups
    if ((!-e "$sample_tumor_bam.pileup")||(-z "$sample_tumor_bam.pileup")) {
	print VARSCAN_SH "samtools mpileup -q 1 -f $opt{GENOME} $sample_tumor_bam > $sample_tumor_bam.pileup\n";
    }
    if ((!-e "$sample_ref_bam.pileup")||(-z "$sample_ref_bam.pileup")) {
	print VARSCAN_SH "samtools mpileup -q 1 -f $opt{GENOME} $sample_ref_bam > $sample_ref_bam.pileup\n";
    }
    print VARSCAN_SH "echo \"End Pileup\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/varscan.log\n\n";

    # run varscan
    print VARSCAN_SH "echo \"Start Varscan\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/varscan.log\n";
    print VARSCAN_SH "java -Xmx12g -jar $opt{VARSCAN_PATH} somatic $sample_ref_bam.pileup $sample_tumor_bam.pileup $sample_tumor_name $opt{VARSCAN_SETTINGS} --output-vcf 1\n\n";

    # postprocessing
    print VARSCAN_SH "java -Xmx12g -jar $opt{VARSCAN_PATH} processSomatic $sample_tumor_name.indel.vcf $opt{VARSCAN_POSTSETTINGS}\n";
    print VARSCAN_SH "java -Xmx12g -jar $opt{VARSCAN_PATH} processSomatic $sample_tumor_name.snp.vcf $opt{VARSCAN_POSTSETTINGS}\n\n";
    
    # merge varscan hc snps and indels
    print VARSCAN_SH "java -Xmx6G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T CombineVariants -R $opt{GENOME} -o $sample_tumor_name.merged.Somatic.hc.vcf -V $sample_tumor_name.snp.Somatic.hc.vcf -V $sample_tumor_name.indel.Somatic.hc.vcf\n";
    print VARSCAN_SH "sed -i 's/SSC/VS_SSC/' $sample_tumor_name.merged.Somatic.hc.vcf\n\n"; # to resolve merge conflict with FB vcfs
    
    # Check varscan completed
    print VARSCAN_SH "touch $log_dir/varscan.done\n\n"; ## Check on complete output!!!
    print VARSCAN_SH "echo \"End Varscan\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/varscan.log\n";

    close VARSCAN_SH;

    ## Run job
    if ( @running_jobs ){
        system "qsub -q $opt{VARSCAN_QUEUE} -pe threaded $opt{VARSCAN_THREADS} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
    } else {
        system "qsub -q $opt{VARSCAN_QUEUE} -pe threaded $opt{VARSCAN_THREADS} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id $bash_file";
    }

    return $job_id;
}

sub runFreeBayes {
    my ($sample_tumor, $sample_tumor_name, $out_dir, $job_dir, $log_dir, $sample_tumor_bam, $sample_ref_bam, $running_jobs, $opt) = (@_);
    my @running_jobs = @{$running_jobs};
    my %opt = %{$opt};
    my $freebayes_out_dir = "$out_dir/freebayes";

    ## Create output dir
    if(! -e $freebayes_out_dir){
	make_path($freebayes_out_dir) or die "Couldn't create directory: $freebayes_out_dir\n";
    }

    ## Skip varscan if .done file exist
    if (-e "$log_dir/freebayes.done"){
	warn "WARNING: $log_dir/freebayes.done, skipping \n";
	return;
    }

    ## Setup freebayes bash script and output dir
    my $job_id = "FB_".$sample_tumor."_".get_job_id();
    my $bash_file = $job_dir."/".$job_id.".sh";

    # Create bash script
    open FREEBAYES_SH, ">$bash_file" or die "cannot open file $bash_file \n";
    print FREEBAYES_SH "#!/bin/bash\n\n";
    print FREEBAYES_SH "echo \"Start Freebayes\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/freebayes.log\n";

    # Run freebayes
    print FREEBAYES_SH "$opt{FREEBAYES_PATH}/freebayes -f $opt{GENOME} -t $opt{FREEBAYES_TARGETS} $opt{FREEBAYES_SETTINGS} $sample_ref_bam $sample_tumor_bam > $freebayes_out_dir/$sample_tumor_name.vcf\n\n";

    # Uniqify freebayes output 
    print FREEBAYES_SH "uniq $freebayes_out_dir/$sample_tumor_name.vcf > $freebayes_out_dir/$sample_tumor_name.uniq.vcf\n";
    print FREEBAYES_SH "mv $freebayes_out_dir/$sample_tumor_name.uniq.vcf $freebayes_out_dir/$sample_tumor_name.vcf\n\n";

    # get sample ids
    print FREEBAYES_SH "sample_R=`grep -P \"^#CHROM\" $freebayes_out_dir/$sample_tumor_name.vcf | cut -f 10`\n";
    print FREEBAYES_SH "sample_T=`grep -P \"^#CHROM\" $freebayes_out_dir/$sample_tumor_name.vcf | cut -f 11`\n\n";

    # annotate somatic and germline scores
    print FREEBAYES_SH "$opt{VCFSAMPLEDIFF_PATH}/vcfsamplediff VT \$sample_R \$sample_T $freebayes_out_dir/$sample_tumor_name.vcf > $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf\n";
    print FREEBAYES_SH "sed -i 's/SSC/FB_SSC/' $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf\n"; # to resolve merge conflicts with varscan vcfs
    print FREEBAYES_SH "grep -P \"^#\" $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf > $freebayes_out_dir/$sample_tumor_name\_germline.vcf\n";
    print FREEBAYES_SH "grep -P \"^#\" $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf > $freebayes_out_dir/$sample_tumor_name\_somatic.vcf\n";
    print FREEBAYES_SH "grep -i \"VT=germline\" $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf >> $freebayes_out_dir/$sample_tumor_name\_germline.vcf\n";
    print FREEBAYES_SH "grep -i \"VT=somatic\" $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf >> $freebayes_out_dir/$sample_tumor_name\_somatic.vcf\n";
    print FREEBAYES_SH "rm $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf\n\n";

    # Filter
    print FREEBAYES_SH "cat $freebayes_out_dir/$sample_tumor_name\_somatic.vcf | $opt{BIOVCF_PATH}/bio-vcf $opt{FREEBAYES_SOMATICFILTER} > $freebayes_out_dir/$sample_tumor_name\_somatic_filtered.vcf\n";
    print FREEBAYES_SH "cat $freebayes_out_dir/$sample_tumor_name\_germline.vcf | $opt{BIOVCF_PATH}/bio-vcf $opt{FREEBAYES_GERMLINEFILTER} > $freebayes_out_dir/$sample_tumor_name\_germline_filtered.vcf\n\n";

    print FREEBAYES_SH "touch $log_dir/freebayes.done\n\n"; ## Check on complete output!!!
    print FREEBAYES_SH "echo \"End Freebayes\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/freebayes.log\n";

    close FREEBAYES_SH;

    # Run job
    if ( @running_jobs ){
	system "qsub -q $opt{FREEBAYES_QUEUE} -pe threaded $opt{FREEBAYES_THREADS} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
    } else {
	system "qsub -q $opt{FREEBAYES_QUEUE} -pe threaded $opt{FREEBAYES_THREADS} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id $bash_file";
    }
    return $job_id;
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
    if(! $opt{SOMVAR_STRELKA}){ die "ERROR: NO SOMVAR_STRELKA found in .ini file\n" }
    if($opt{SOMVAR_STRELKA} eq "yes"){
	if(! $opt{STRELKA_PATH}){ die "ERROR: NO STRELKA_PATH found in .ini file\n" }
	if(! $opt{STRELKA_INI}){ die "ERROR: NO STRELKA_INI found in .ini file\n" }
	if(! $opt{STRELKA_QUEUE}){ die "ERROR: NO STRELKA_QUEUE found in .ini file\n" }
	if(! $opt{STRELKA_THREADS}){ die "ERROR: NO STRELKA_THREADS found in .ini file\n" }
    }
    if(! $opt{SOMVAR_VARSCAN}){ die "ERROR: NO SOMVAR_VARSCAN found in .ini file\n" }
    if($opt{SOMVAR_VARSCAN} eq "yes"){
	if(! $opt{VARSCAN_PATH}){ die "ERROR: NO VARSCAN_PATH found in .ini file\n" }
	if(! $opt{VARSCAN_QUEUE}){ die "ERROR: NO VARSCAN_QUEUE found in .ini file\n" }
	if(! $opt{VARSCAN_THREADS}){ die "ERROR: NO VARSCAN_THREADS found in .ini file\n" }
	if(! $opt{VARSCAN_SETTINGS}){ die "ERROR: NO VARSCAN_SETTINGS found in .ini file\n" }
	if(! $opt{VARSCAN_POSTSETTINGS}){ die "ERROR: NO VARSCAN_POSTSETTINGS found in .ini file\n" }
    }
    if(! $opt{SOMVAR_FREEBAYES}){ die "ERROR: NO SOMVAR_FREEBAYES found in .ini file\n" }
    if($opt{SOMVAR_FREEBAYES} eq "yes"){
	if(! $opt{FREEBAYES_PATH}){ die "ERROR: NO found in .ini file\n" }
	if(! $opt{VCFSAMPLEDIFF_PATH}){ die "ERROR: NO found in .ini file\n" }
	if(! $opt{BIOVCF_PATH}){ die "ERROR: NO found in .ini file\n" }
	if(! $opt{FREEBAYES_QUEUE}){ die "ERROR: NO found in .ini file\n" }
	if(! $opt{FREEBAYES_THREADS}){ die "ERROR: NO found in .ini file\n" }
	if(! $opt{FREEBAYES_TARGETS}){ die "ERROR: NO found in .ini file\n" }
	if(! $opt{FREEBAYES_SETTINGS}){ die "ERROR: NO found in .ini file\n" }
	if(! $opt{FREEBAYES_SOMATICFILTER}){ die "ERROR: NO found in .ini file\n" }
	if(! $opt{FREEBAYES_GERMLINEFILTER}){ die "ERROR: NO found in .ini file\n" }
    }
    if(! $opt{SOMVARMERGE_QUEUE}){ die "ERROR: NO SOMVARMERGE_QUEUE found in .ini file\n" }
    if(! $opt{SOMVARMERGE_THREADS}){ die "ERROR: NO found in .ini file\n" }
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
