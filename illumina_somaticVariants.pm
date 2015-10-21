#!/usr/bin/perl -w

#########################################################
### illumina_somaticVariants.pm
### - Run somatic variant callers
###   - Varscan, Strelka, FreeBayes
### - Merge and annotate somatic high confidence calls.
###
### Author: R.F.Ernst
#########################################################

package illumina_somaticVariants;

use strict;
use POSIX qw(tmpnam);
use File::Path qw(make_path);


sub parseSamples {
    ###
    # Parse sample names
    # Expects CPCT samples (CPCT........T/R)
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my %somatic_samples;

    foreach my $sample (@{$opt{SAMPLES}}){
	# Parse cpct samples based on regular expression defining two groups, sample name and sample origin.
	my ($cpct_name,$origin) = ($sample =~ /$opt{SOMVAR_REGEX}/);

	if ( (! $cpct_name) || (! $origin) ){
	    print "WARNING: $sample is not passing somatic samplename parsing, skipping \n\n";
	    next;
	}
	
	# Reference sample
	if ($origin =~ m/R.*/){
	    if ($somatic_samples{$cpct_name}{"ref"}){
		print "\t WARNING: $cpct_name has multiple reference samples, using: $somatic_samples{$cpct_name}{'ref'} \n\n";
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
    ###
    # Run somatic variant callers
    # Merge and annotate high confidence calls
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my @merge_somvar_jobs;
    ### Loop over tumor samples
    foreach my $sample (keys(%{$opt{SOMATIC_SAMPLES}})){

	# Check correct sample ref
	if (! $opt{SOMATIC_SAMPLES}{$sample}{'ref'}){
	    print "WARNING: No ref sample for $sample, skipping \n";
	    next;
	}

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
	    print "\n$sample \t $sample_ref_bam \t $sample_tumor_bam \n";

	    ## Skip Somatic Callers if .done file exist
	    if (-e "$sample_tumor_log_dir/$sample_tumor_name.done"){
		print "WARNING: $sample_tumor_log_dir/$sample_tumor_name.done exists, skipping \n";
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
	    if($opt{SOMVAR_MUTECT} eq "yes"){
		print "\n###SCHEDULING MUTECT####\n";
		my $mutect_job = runMutect($sample_tumor, $sample_tumor_name, $sample_tumor_out_dir, $sample_tumor_job_dir, $sample_tumor_log_dir, $sample_tumor_bam, $sample_ref_bam, \@running_jobs, \%opt);
		if($mutect_job){push(@somvar_jobs, $mutect_job)};
	    }
	    ## Merge somatic vcfs
	    print "\n###SCHEDULING MERGE SOMATIC VCFS####\n";

	    my $job_id = "MERGE_".$sample_tumor."_".get_job_id();
	    my $bash_file = $sample_tumor_job_dir."/".$job_id.".sh";

	    open MERGE_SH, ">$bash_file" or die "cannot open file $bash_file \n";
	    print MERGE_SH "#!/bin/bash\n\n";
	    print MERGE_SH "echo \"Start Merge\t\" `date` `uname -n` >> $sample_tumor_log_dir/merge.log\n\n";

	    # Merge vcfs
	    my $invcf;
	    my $outvcf = "$sample_tumor_out_dir/$sample_tumor_name\_merged_somatics.vcf";
	    print MERGE_SH "java -Xmx6G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T CombineVariants -R $opt{GENOME} -o $outvcf --genotypemergeoption uniquify ";
	    if($opt{SOMVAR_STRELKA} eq "yes"){ print MERGE_SH "-V:strelka $sample_tumor_out_dir/strelka/passed.somatic.merged.vcf "; }
	    if($opt{SOMVAR_VARSCAN} eq "yes"){ print MERGE_SH "-V:varscan $sample_tumor_out_dir/varscan/$sample_tumor_name.merged.Somatic.hc.vcf "; }
	    if($opt{SOMVAR_FREEBAYES} eq "yes"){ print MERGE_SH "-V:freebayes $sample_tumor_out_dir/freebayes/$sample_tumor_name\_somatic_filtered.vcf "; }
	    if($opt{SOMVAR_MUTECT} eq "yes"){ print MERGE_SH "-V:mutect $sample_tumor_out_dir/mutect/$sample_tumor_name\_mutect.vcf ";}

	    # Filter vcf on target
	    if($opt{SOMVAR_TARGETS}){
		$invcf = $outvcf;
		$outvcf = "$sample_tumor_out_dir/$sample_tumor_name\_filtered_merged_somatics.vcf";
		print MERGE_SH "\n\njava -Xmx6G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T SelectVariants -R $opt{GENOME} -L $opt{SOMVAR_TARGETS} -V $invcf -o $outvcf\n";
		print MERGE_SH "rm $invcf*";
	    }

	    # Annotate somatic vcf
	    if($opt{SOMVAR_ANNOTATE} eq "yes"){
		$invcf = $outvcf;
		my $preAnnotateVCF = $invcf;
		$outvcf =~ s/.vcf/_snpEff.vcf/;
		print MERGE_SH "\n\njava -Xmx6G -jar $opt{SNPEFF_PATH}/snpEff.jar -c $opt{SNPEFF_PATH}/snpEff.config $opt{ANNOTATE_DB} -v $invcf $opt{ANNOTATE_FLAGS} > $outvcf\n";
		
		## dbsnp
		$invcf = $outvcf;
		my $suffix = "_dbSNP.vcf";
		$outvcf =~ s/.vcf/$suffix/;
		print MERGE_SH "java -Xmx6G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T VariantAnnotator -nt $opt{SOMVARMERGE_THREADS} -R $opt{GENOME} -o $outvcf --variant $invcf --dbsnp $opt{CALLING_DBSNP} --alwaysAppendDbsnpId\n";
		print MERGE_SH "if [ -f $outvcf ]\nthen\n\trm $invcf $invcf.idx \nfi\n";
		
		## cosmic
		$invcf = $outvcf;
		$suffix = "_$opt{ANNOTATE_IDNAME}.vcf";
		$outvcf =~ s/.vcf/$suffix/;
		print MERGE_SH "java -Xmx6G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T VariantAnnotator -nt $opt{SOMVARMERGE_THREADS} -R $opt{GENOME} -o $outvcf --variant $invcf --dbsnp $opt{ANNOTATE_IDDB} --alwaysAppendDbsnpId\n";
		print MERGE_SH "if [ -f $outvcf ]\nthen\n\trm $invcf $invcf.idx \nfi\n";
		
		## Check annotated vcf using the last position
		print MERGE_SH "\nif [ \"\$(tail -n 1 $preAnnotateVCF | cut -f 1,2)\" = \"\$(tail -n 1 $outvcf | cut -f 1,2)\" ]\n";
		print MERGE_SH "then\n";
		print MERGE_SH "\ttouch $sample_tumor_log_dir/$sample_tumor_name.done\n";
		print MERGE_SH "fi\n";
		print MERGE_SH "echo \"END Merge\t\" `date` `uname -n` >> $sample_tumor_log_dir/merge.log\n\n";
		close MERGE_SH;
		
	    } else {
		print MERGE_SH "\nif [ -f $outvcf ]\n";
		print MERGE_SH "then\n";
		print MERGE_SH "\ttouch $sample_tumor_log_dir/$sample_tumor_name.done\n";
		print MERGE_SH "fi\n";
		print MERGE_SH "echo \"END Merge\t\" `date` `uname -n` >> $sample_tumor_log_dir/merge.log\n\n";
		close MERGE_SH;
	    }

	    # Run job
	    if ( @somvar_jobs ){
		system "qsub -q $opt{SOMVARMERGE_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{SOMVARMERGE_THREADS} -P $opt{CLUSTER_PROJECT} -o $sample_tumor_log_dir -e $sample_tumor_log_dir -N $job_id -hold_jid ".join(",",@somvar_jobs)." $bash_file";
	    } else {
		system "qsub -q $opt{SOMVARMERGE_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{SOMVARMERGE_THREADS} -P $opt{CLUSTER_PROJECT} -o $sample_tumor_log_dir -e $sample_tumor_log_dir -N $job_id $bash_file";
	    }
	    push(@merge_somvar_jobs, $job_id);
	}
    }
    return \@merge_somvar_jobs;
}

### Somatic Variant Callers
sub runStrelka {
    my ($sample_tumor, $out_dir, $job_dir, $log_dir, $sample_tumor_bam, $sample_ref_bam, $running_jobs, $opt) = (@_);
    my @running_jobs = @{$running_jobs};
    my %opt = %{$opt};
    my $strelka_out_dir = "$out_dir/strelka";

    ## Skip Strelka if .done file exist
    if (-e "$log_dir/strelka.done"){
	print "WARNING: $log_dir/strelka.done, skipping \n";
	return;
    }

    ## Create strelka bash script
    my $job_id = "STR_".$sample_tumor."_".get_job_id();
    my $bash_file = $job_dir."/".$job_id.".sh";

    open STRELKA_SH, ">$bash_file" or die "cannot open file $bash_file \n";
    print STRELKA_SH "#!/bin/bash\n\n";
    print STRELKA_SH "if [ -f $sample_tumor_bam -a -f $sample_ref_bam ]\n";
    print STRELKA_SH "then\n";
    print STRELKA_SH "\techo \"Start Strelka\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/strelka.log\n\n";

    # Run Strelka
    print STRELKA_SH "\t$opt{STRELKA_PATH}/bin/configureStrelkaWorkflow.pl --tumor $sample_tumor_bam --normal $sample_ref_bam --ref $opt{GENOME} --config $opt{STRELKA_INI} --output-dir $strelka_out_dir\n\n";

    print STRELKA_SH "\tcd $strelka_out_dir\n";
    print STRELKA_SH "\tmake -j 8\n\n";

    # Check strelka completed
    print STRELKA_SH "\tif [ -f $strelka_out_dir/task.complete ]\n";
    print STRELKA_SH "\tthen\n";
    print STRELKA_SH "\t\tjava -Xmx6G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T CombineVariants -R $opt{GENOME} --genotypemergeoption unsorted -o passed.somatic.merged.vcf -V results/passed.somatic.snvs.vcf -V results/passed.somatic.indels.vcf \n";
    print STRELKA_SH "\t\tperl -p -e 's/\\t([A-Z][A-Z]:)/\\tGT:\$1/g' passed.somatic.merged.vcf | perl -p -e 's/(:T[UO]R?)\\t/\$1\\t0\\/0:/g' | perl -p -e 's/(:\\d+,\\d+)\\t/\$1\\t0\\/1:/g' | perl -p -e 's/(#CHROM.*)/##StrelkaGATKCompatibility=Added GT fields to strelka calls for gatk compatibility.\\n\$1/g' > temp.vcf\n";
    print STRELKA_SH "\t\tmv temp.vcf passed.somatic.merged.vcf\n";
    print STRELKA_SH "\t\trm -r chromosomes/ \n";
    print STRELKA_SH "\t\ttouch $log_dir/strelka.done\n";
    print STRELKA_SH "\tfi\n\n";
    print STRELKA_SH "\techo \"End Strelka\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/strelka.log\n\n";
    
    print STRELKA_SH "else\n";
    print STRELKA_SH "\techo \"ERROR: $sample_tumor_bam or $sample_ref_bam does not exist.\" >&2\n";
    print STRELKA_SH "fi\n";
    
    close STRELKA_SH;

    ## Run job
    if ( @running_jobs ){
	system "qsub -q $opt{STRELKA_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{STRELKA_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
    } else {
	system "qsub -q $opt{STRELKA_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{STRELKA_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id $bash_file";
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
	print "WARNING: $log_dir/varscan.done exists, skipping \n";
	return;
    }

    ## Setup varscan bash script and output dir
    my $job_id = "VS_".$sample_tumor."_".get_job_id();
    my $bash_file = $job_dir."/".$job_id.".sh";

    ## Create bash script
    open VARSCAN_SH, ">$bash_file" or die "cannot open file $bash_file \n";
    print VARSCAN_SH "#!/bin/bash\n\n";
    print VARSCAN_SH "cd $varscan_out_dir\n";
    print VARSCAN_SH "if [ -f $sample_tumor_bam -a -f $sample_ref_bam ]\n";
    print VARSCAN_SH "then\n";
    
    # run pileups
    print VARSCAN_SH "\techo \"Start Pileup\t\" `date` \"$sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/varscan.log\n";
    print VARSCAN_SH "\tPATH=$opt{SAMTOOLS_PATH}:\$PATH\n";
    print VARSCAN_SH "\texport PATH\n";
    if ((!-e "$sample_tumor_bam.pileup")||(-z "$sample_tumor_bam.pileup")) {
	#print VARSCAN_SH "\t$opt{SAMTOOLS_PATH}/samtools mpileup -q 1 -f $opt{GENOME} -l $opt{SOMVAR_TARGETS} $sample_tumor_bam > $sample_tumor_bam.pileup_samtools\n";
	print VARSCAN_SH "\t$opt{SAMBAMBA_PATH}/sambamba mpileup -t $opt{VARSCAN_THREADS} --tmpdir=$opt{OUTPUT_DIR}/tmp/ -L $opt{SOMVAR_TARGETS} -o $sample_tumor_bam.pileup $sample_tumor_bam --samtools \"-q 1 -f $opt{GENOME}\"\n";
    }
    if ((!-e "$sample_ref_bam.pileup")||(-z "$sample_ref_bam.pileup")) {
	#print VARSCAN_SH "\t$opt{SAMTOOLS_PATH}/samtools mpileup -q 1 -f $opt{GENOME} -l $opt{SOMVAR_TARGETS} $sample_ref_bam > $sample_ref_bam.pileup_samtools\n";
	print VARSCAN_SH "\t$opt{SAMBAMBA_PATH}/sambamba mpileup -t $opt{VARSCAN_THREADS} --tmpdir=$opt{OUTPUT_DIR}/tmp/ -L $opt{SOMVAR_TARGETS} -o $sample_ref_bam.pileup $sample_ref_bam --samtools \"-q 1 -f $opt{GENOME}\"\n";
    }
    print VARSCAN_SH "\techo \"End Pileup\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/varscan.log\n\n";

    # run varscan
    print VARSCAN_SH "\techo \"Start Varscan\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/varscan.log\n";
    print VARSCAN_SH "\tjava -Xmx12g -jar $opt{VARSCAN_PATH} somatic $sample_ref_bam.pileup $sample_tumor_bam.pileup $sample_tumor_name $opt{VARSCAN_SETTINGS} --output-vcf 1\n\n";

    # postprocessing
    print VARSCAN_SH "\tjava -Xmx12g -jar $opt{VARSCAN_PATH} processSomatic $sample_tumor_name.indel.vcf $opt{VARSCAN_POSTSETTINGS}\n";
    print VARSCAN_SH "\tjava -Xmx12g -jar $opt{VARSCAN_PATH} processSomatic $sample_tumor_name.snp.vcf $opt{VARSCAN_POSTSETTINGS}\n\n";
    
    # merge varscan hc snps and indels
    print VARSCAN_SH "\tjava -Xmx6G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T CombineVariants -R $opt{GENOME} --genotypemergeoption unsorted -o $sample_tumor_name.merged.Somatic.hc.vcf -V $sample_tumor_name.snp.Somatic.hc.vcf -V $sample_tumor_name.indel.Somatic.hc.vcf\n";
    print VARSCAN_SH "\tsed -i 's/SSC/VS_SSC/' $sample_tumor_name.merged.Somatic.hc.vcf\n\n"; # to resolve merge conflict with FB vcfs
    
    # Check varscan completed
    print VARSCAN_SH "\tif [ -f $sample_tumor_name.merged.Somatic.hc.vcf ]\n";
    print VARSCAN_SH "\tthen\n";
    print VARSCAN_SH "\t\ttouch $log_dir/varscan.done\n";
    print VARSCAN_SH "\tfi\n\n";
    print VARSCAN_SH "\techo \"End Varscan\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/varscan.log\n";

    print VARSCAN_SH "else\n";
    print VARSCAN_SH "\techo \"ERROR: $sample_tumor_bam or $sample_ref_bam does not exist.\" >&2\n";
    print VARSCAN_SH "fi\n";
    close VARSCAN_SH;

    ## Run job
    if ( @running_jobs ){
        system "qsub -q $opt{VARSCAN_QUEUE} -pe threaded $opt{VARSCAN_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
    } else {
        system "qsub -q $opt{VARSCAN_QUEUE} -pe threaded $opt{VARSCAN_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id $bash_file";
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

    ## Skip freebayes if .done file exist
    if (-e "$log_dir/freebayes.done"){
	print "WARNING: $log_dir/freebayes.done exists, skipping \n";
	return;
    }

    ## Setup freebayes bash script and output dir
    my $job_id = "FB_".$sample_tumor."_".get_job_id();
    my $bash_file = $job_dir."/".$job_id.".sh";

    # Create bash script
    open FREEBAYES_SH, ">$bash_file" or die "cannot open file $bash_file \n";
    print FREEBAYES_SH "#!/bin/bash\n\n";
    print FREEBAYES_SH "if [ -f $sample_tumor_bam -a -f $sample_ref_bam ]\n";
    print FREEBAYES_SH "then\n";
    
    print FREEBAYES_SH "\techo \"Start Freebayes\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/freebayes.log\n";

    # Run freebayes
    print FREEBAYES_SH "\t$opt{FREEBAYES_PATH}/freebayes -f $opt{GENOME} -t $opt{SOMVAR_TARGETS} $opt{FREEBAYES_SETTINGS} $sample_ref_bam $sample_tumor_bam > $freebayes_out_dir/$sample_tumor_name.vcf\n\n";

    # Uniqify freebayes output 
    print FREEBAYES_SH "\tuniq $freebayes_out_dir/$sample_tumor_name.vcf > $freebayes_out_dir/$sample_tumor_name.uniq.vcf\n";
    print FREEBAYES_SH "\tmv $freebayes_out_dir/$sample_tumor_name.uniq.vcf $freebayes_out_dir/$sample_tumor_name.vcf\n\n";

    # get sample ids
    print FREEBAYES_SH "\tsample_R=`grep -P \"^#CHROM\" $freebayes_out_dir/$sample_tumor_name.vcf | cut -f 10`\n";
    print FREEBAYES_SH "\tsample_T=`grep -P \"^#CHROM\" $freebayes_out_dir/$sample_tumor_name.vcf | cut -f 11`\n\n";

    # annotate somatic and germline scores
    print FREEBAYES_SH "\t$opt{VCFSAMPLEDIFF_PATH}/vcfsamplediff VT \$sample_R \$sample_T $freebayes_out_dir/$sample_tumor_name.vcf > $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf\n";
    print FREEBAYES_SH "\tsed -i 's/SSC/FB_SSC/' $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf\n"; # to resolve merge conflicts with varscan vcfs
    print FREEBAYES_SH "\tgrep -P \"^#\" $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf > $freebayes_out_dir/$sample_tumor_name\_germline.vcf\n";
    print FREEBAYES_SH "\tgrep -P \"^#\" $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf > $freebayes_out_dir/$sample_tumor_name\_somatic.vcf\n";
    print FREEBAYES_SH "\tgrep -i \"VT=germline\" $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf >> $freebayes_out_dir/$sample_tumor_name\_germline.vcf\n";
    print FREEBAYES_SH "\tgrep -i \"VT=somatic\" $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf >> $freebayes_out_dir/$sample_tumor_name\_somatic.vcf\n";
    print FREEBAYES_SH "\trm $freebayes_out_dir/$sample_tumor_name\_VTannot.vcf\n\n";

    # Filter
    print FREEBAYES_SH "\tcat $freebayes_out_dir/$sample_tumor_name\_somatic.vcf | $opt{BIOVCF_PATH}/bio-vcf $opt{FREEBAYES_SOMATICFILTER} > $freebayes_out_dir/$sample_tumor_name\_somatic_filtered.vcf\n";
    print FREEBAYES_SH "\tcat $freebayes_out_dir/$sample_tumor_name\_germline.vcf | $opt{BIOVCF_PATH}/bio-vcf $opt{FREEBAYES_GERMLINEFILTER} > $freebayes_out_dir/$sample_tumor_name\_germline_filtered.vcf\n\n";
    
    #Check freebayes completed
    print FREEBAYES_SH "\tif [ -f $freebayes_out_dir/$sample_tumor_name\_somatic_filtered.vcf -a -f $freebayes_out_dir/$sample_tumor_name\_germline_filtered.vcf ]\n";
    print FREEBAYES_SH "\tthen\n";
    print FREEBAYES_SH "\t\ttouch $log_dir/freebayes.done\n\n"; ## Check on complete output!!!
    print FREEBAYES_SH "\tfi\n";
    
    print FREEBAYES_SH "\techo \"End Freebayes\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/freebayes.log\n";
    print FREEBAYES_SH "else\n";
    print FREEBAYES_SH "\techo \"ERROR: $sample_tumor_bam or $sample_ref_bam does not exist.\" >&2\n";
    print FREEBAYES_SH "fi\n";

    close FREEBAYES_SH;

    # Run job
    if ( @running_jobs ){
	system "qsub -q $opt{FREEBAYES_QUEUE} -pe threaded $opt{FREEBAYES_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
    } else {
	system "qsub -q $opt{FREEBAYES_QUEUE} -pe threaded $opt{FREEBAYES_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -m a -M $opt{MAIL} -o $log_dir -e $log_dir -N $job_id $bash_file";
    }
    return $job_id;
}

sub runMutect {
    my ($sample_tumor, $sample_tumor_name, $out_dir, $job_dir, $log_dir, $sample_tumor_bam, $sample_ref_bam, $running_jobs, $opt) = (@_);
    my @running_jobs = @{$running_jobs};
    my %opt = %{$opt};
    my $mutect_out_dir = "$out_dir/mutect";
    my $mutect_tmp_dir = "$mutect_out_dir/tmp";

    ## Create output and tmp dir
    if(! -e $mutect_out_dir){
	make_path($mutect_out_dir) or die "Couldn't create directory: $mutect_out_dir\n";
    }
    if(! -e $mutect_tmp_dir){
	make_path($mutect_tmp_dir) or die "Couldn't create directory: $mutect_tmp_dir\n";
    }
    
    ## Skip Mutect if .done file exist
    if (-e "$log_dir/mutect.done"){
	print "WARNING: $log_dir/mutect.done, skipping \n";
	return;
    }

    ## Build Queue command 
    ## wait for GATK-MuTect integration, see http://gatkforums.broadinstitute.org/discussion/comment/24614#Comment_24614
    #my $javaMem = $opt{MUTECT_MASTERTHREADS} * $opt{MUTECT_MEM};
    #my $javaJobMem = $opt{MUTECT_THREADS} * $opt{MUTECT_MEM};
    #my $command = "java -Xmx".$javaMem."G -Xms".$opt{MUTECT_MEM}."G -jar $opt{QUEUE_PATH}/Queue.jar ";
    #$command .= "-jobQueue $opt{MUTECT_QUEUE} -jobNative \"-pe threaded $opt{MUTECT_THREADS} -P $opt{CLUSTER_PROJECT}\" -jobRunner GridEngine -jobReport $log_dir/mutect.jobReport.txt -memLimit $javaJobMem "; #Queue options
    #$command .= "-S $opt{MUTECT_SCALA} ";
    #$command .= "-R $opt{GENOME} -O $sample_tumor_name -mem $opt{MUTECT_MEM} -nsc $opt{MUTECT_SCATTER} ";
    #$command .= "-tb $sample_tumor_bam -nb $sample_ref_bam ";
    #$command .= "-D $opt{CALLING_DBSNP} -C $opt{MUTECT_COSMIC} ";
    
    ### Optional settings
    #if ( $opt{CALLING_TARGETS} ) {
	#$command .= "-L $opt{CALLING_TARGETS} ";
	#if ( $opt{CALLING_INTERVALPADDING} ) {
	    #$command .= "-ip $opt{CALLING_INTERVALPADDING} ";
	#}
    #}
    #if($opt{QUEUE_RETRY} eq 'yes'){
	#$command  .= "-retry 1 ";
    #}
    
    ## Set run option
    #$command .= "-run";
    
    ### Mutect .jar command
    my $javaJobMem = $opt{MUTECT_THREADS} * $opt{MUTECT_MEM};
    my $command = "java -Xmx".$javaJobMem."G -jar $opt{MUTECT_PATH}/mutect.jar -T MuTect ";
    $command .= "-R $opt{GENOME} --cosmic $opt{MUTECT_COSMIC} --dbsnp $opt{CALLING_DBSNP} --intervals $opt{CALLING_TARGETS} ";
    $command .= "--input_file:normal $sample_ref_bam --input_file:tumor $sample_tumor_bam ";
    $command .= "--out call_stats.out --vcf $sample_tumor_name\_mutect.vcf";
    ## Create mutect bash script
    my $job_id = "MUT_".$sample_tumor."_".get_job_id();
    my $bash_file = $job_dir."/".$job_id.".sh";
    open MUTECT_SH, ">$bash_file" or die "cannot open file $bash_file \n";
    print MUTECT_SH "#!/bin/bash\n\n";
    print MUTECT_SH "if [ -f $sample_tumor_bam -a -f $sample_ref_bam ]\n";
    print MUTECT_SH "then\n";
    print MUTECT_SH "\techo \"Start Mutect\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/mutect.log\n\n";

    # Run Mutect
    print MUTECT_SH "\tcd $mutect_tmp_dir\n";
    print MUTECT_SH "\t$command\n";
    
    # Filter Mutect result
    $command = "cat $sample_tumor_name\_mutect.vcf | java -Xmx".$javaJobMem."G -jar $opt{SNPEFF_PATH}/SnpSift.jar filter \"( na FILTER ) | (FILTER = 'PASS')\" > $sample_tumor_name\_mutect_passed.vcf \n";
    print MUTECT_SH "\t$command\n\n";
    # Check Mutect completed
    print MUTECT_SH "\tif [ -f $sample_tumor_name\_mutect.vcf -a -f $sample_tumor_name\_mutect_passed.vcf ]\n";
    print MUTECT_SH "\tthen\n";
    print MUTECT_SH "\t\tmv $sample_tumor_name\_mutect.vcf $mutect_out_dir/\n";
    print MUTECT_SH "\t\tmv $sample_tumor_name\_mutect.vcf.idx $mutect_out_dir/\n";
    print MUTECT_SH "\t\tmv $sample_tumor_name\_mutect_passed.vcf $mutect_out_dir/\n";
    print MUTECT_SH "\t\tmv $sample_tumor_name\_mutect_passed.vcf.idx $mutect_out_dir/\n";
    print MUTECT_SH "\t\tcd $mutect_out_dir/\n";
    print MUTECT_SH "\t\trm -r tmp/\n";
    print MUTECT_SH "\t\ttouch $log_dir/mutect.done\n";
    print MUTECT_SH "\tfi\n\n";
    print MUTECT_SH "\techo \"End Mutect\t\" `date` \"\t $sample_ref_bam \t $sample_tumor_bam\t\" `uname -n` >> $log_dir/mutect.log\n\n";

    print MUTECT_SH "else\n";
    print MUTECT_SH "\techo \"ERROR: $sample_tumor_bam or $sample_ref_bam does not exist.\" >&2\n";
    print MUTECT_SH "fi\n";

    close MUTECT_SH;

    ## Run job
    if ( @running_jobs ){
	#system "qsub -q $opt{MUTECT_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{MUTECT_MASTERTHREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
	system "qsub -q $opt{MUTECT_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{MUTECT_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id -hold_jid ".join(",",@running_jobs)." $bash_file";
    } else {
	#system "qsub -q $opt{MUTECT_MASTERQUEUE} -m a -M $opt{MAIL} -pe threaded $opt{MUTECT_MASTERTHREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id $bash_file";
	system "qsub -q $opt{MUTECT_QUEUE} -m a -M $opt{MAIL} -pe threaded $opt{MUTECT_THREADS} -R $opt{CLUSTER_RESERVATION} -P $opt{CLUSTER_PROJECT} -o $log_dir -e $log_dir -N $job_id $bash_file";
    }

    return $job_id;
}

############
sub get_job_id {
    my $id = tmpnam();
    $id=~s/\/tmp\/file//;
    return $id;
}
############

1;
