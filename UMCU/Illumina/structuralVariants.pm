#!/usr/bin/perl -w

##################################################################
### illumina_structuralVariants.pm
### - Run structural variant callers Delly and Manta
###
### Author: R.F.Ernst , M. van Roosmalen ,H.H.D.Kerstens
##################################################################

package UMCU::Illumina::structuralVariants;

use strict;
use POSIX qw(tmpnam);
use File::Path qw(make_path);
use lib "$FindBin::Bin"; #locates pipeline directory
use UMCU::Illumina::sge;

sub runStructuralVariantCallers {
    ###
    # Run structural variant callers
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my @sv_jobs;
    my @sv_samples;

    if($opt{SV_DELLY} eq "yes"){
	my $delly_jobs = runDelly(\%opt);
	push(@sv_jobs, @{$delly_jobs});
    }
    
    if($opt{SV_MANTA} eq "yes"){
	my $manta_jobs = runManta(\%opt);
	push(@sv_jobs, @{$manta_jobs});
    }
    return(\@sv_jobs);
}

sub runManta {
    ###
    # Run structural variant caller Manta
    ###
    my $configuration = shift;
    my %opt = %{$configuration};
    my @manta_jobs;
    my %somatic_samples = %{$opt{SOMATIC_SAMPLES}};
    my @single_samples = @{$opt{SINGLE_SAMPLES}};
    
    ### Run single samples
    foreach my $sample (@single_samples){
	# Setup output, log and job directories.
	my $manta_out_dir = "$opt{OUTPUT_DIR}/structuralVariants/manta/$sample/";
	my $manta_log_dir = "$opt{OUTPUT_DIR}/structuralVariants/manta/logs/";
	my $manta_job_dir = "$opt{OUTPUT_DIR}/structuralVariants/manta/jobs/";

	if(! -e $manta_out_dir){make_path($manta_out_dir) or die "Couldn't create directory: $manta_out_dir\n"}
	if(! -e $manta_log_dir){make_path($manta_log_dir) or die "Couldn't create directory: $manta_log_dir\n"}
	if(! -e $manta_job_dir){make_path($manta_job_dir) or die "Couldn't create directory: $manta_job_dir\n"}
	
	# Skip Manta if .done file exists
	if (-e "$manta_log_dir/SV_MANTA_$sample.done"){
	    print "WARNING: $manta_log_dir/SV_MMANTA_$sample.done exists, skipping \n";
	} else {
	    # Store running jobs
	    my @running_jobs;
	    if (@{$opt{RUNNING_JOBS}->{$sample}}){ push(@running_jobs, @{$opt{RUNNING_JOBS}->{$sample}}) }

	    # Setup manta commands
	    my $sample_bam = "$opt{OUTPUT_DIR}/$sample/mapping/$opt{BAM_FILES}->{$sample}";
	    my $config_manta = "$opt{MANTA_PATH}/configManta.py --referenceFasta $opt{GENOME} --runDir $manta_out_dir --bam $sample_bam ";
	    my $run_manta = "$manta_out_dir/runWorkflow.py -m local -j $opt{MANTA_THREADS} ";
	    
	    # Create manta bash script
	    my $job_id = "SV_MANTA_$sample\_".get_job_id();
	    my $bashFile = $manta_job_dir.$job_id.".sh";
    
	    open MANTA_SH, ">$bashFile" or die "cannot open file $bashFile \n";
	    print MANTA_SH "#!/bin/bash\n\n";
	    print MANTA_SH "bash $opt{CLUSTER_PATH}/settings.sh\n\n";
	    print MANTA_SH "cd $opt{OUTPUT_DIR}\n\n";

	    print MANTA_SH "$config_manta\n";
	    print MANTA_SH "$run_manta\n\n";

	    print MANTA_SH "if [ -s $manta_out_dir/results/variants/diploidSV.vcf.gz.tbi ]\n";
	    print MANTA_SH "then\n";
	    print MANTA_SH "\ttouch $manta_log_dir/SV_MANTA_$sample.done\n";
	    print MANTA_SH "fi\n";
    
	    my $qsub = &qsubTemplate(\%opt,"MANTA");
	    if (@running_jobs){
		system "$qsub -o $manta_log_dir/$job_id.out -e $manta_log_dir/$job_id.err -N $job_id -hold_jid ".join(",",@running_jobs)." $bashFile";
	    } else {
		system "$qsub -o $manta_log_dir/$job_id.out -e $manta_log_dir/$job_id.err -N $job_id $bashFile";
	    }
	    push(@manta_jobs, $job_id);
	}
    }
    
    ### Run somatic samples
    foreach my $sample (keys %somatic_samples){
	foreach my $sample_tumor (@{$somatic_samples{$sample}{'tumor'}}){
	    foreach my $sample_ref (@{$somatic_samples{$sample}{'ref'}}){
		# Setup output, log and job directories.
		my $sample_tumor_name = "$sample_ref\_$sample_tumor";
		my $manta_out_dir = "$opt{OUTPUT_DIR}/structuralVariants/manta/$sample_tumor_name/";
		my $manta_log_dir = "$opt{OUTPUT_DIR}/structuralVariants/manta/logs/";
		my $manta_job_dir = "$opt{OUTPUT_DIR}/structuralVariants/manta/jobs/";

		if(! -e $manta_out_dir){make_path($manta_out_dir) or die "Couldn't create directory: $manta_out_dir\n"}
		if(! -e $manta_log_dir){make_path($manta_log_dir) or die "Couldn't create directory: $manta_log_dir\n"}
		if(! -e $manta_job_dir){make_path($manta_job_dir) or die "Couldn't create directory: $manta_job_dir\n"}
		
		# Skip Manta if .done file exists
		if (-e "$manta_log_dir/SV_MANTA_$sample_tumor_name.done"){
		    print "WARNING: $manta_log_dir/SV_MANTA_$sample_tumor_name.done exists, skipping \n";
		} else {
		    #Store running jobs
		    my @running_jobs;
		    if (@{$opt{RUNNING_JOBS}->{$sample_tumor}}){ push(@running_jobs, @{$opt{RUNNING_JOBS}->{$sample_tumor}}) }
		    if (@{$opt{RUNNING_JOBS}->{$sample_ref}}){ push(@running_jobs, @{$opt{RUNNING_JOBS}->{$sample_ref}}) }
		
		    # Setup manta commands
		    my $ref_bam = "$opt{OUTPUT_DIR}/$sample_ref/mapping/$opt{BAM_FILES}->{$sample_ref}";
		    my $tumor_bam = "$opt{OUTPUT_DIR}/$sample_tumor/mapping/$opt{BAM_FILES}->{$sample_tumor}";
		    my $config_manta = "$opt{MANTA_PATH}/configManta.py --referenceFasta $opt{GENOME} --runDir $manta_out_dir --normalBam $ref_bam --tumorBam $tumor_bam ";
		    my $run_manta = "$manta_out_dir/runWorkflow.py -m local -j $opt{MANTA_THREADS} ";
		
		    # Create manta bash script
		    my $job_id = "SV_MANTA_$sample_tumor_name\_".get_job_id();
		    my $bashFile = $manta_job_dir.$job_id.".sh";

		    open MANTA_SH, ">$bashFile" or die "cannot open file $bashFile \n";
		    print MANTA_SH "#!/bin/bash\n\n";
		    print MANTA_SH "bash $opt{CLUSTER_PATH}/settings.sh\n\n";
		    print MANTA_SH "cd $opt{OUTPUT_DIR}\n\n";

		    print MANTA_SH "$config_manta\n";
		    print MANTA_SH "$run_manta\n\n";

		    print MANTA_SH "if [ -s $manta_out_dir/results/variants/diploidSV.vcf.gz.tbi -a -s $manta_out_dir/results/variants/somaticSV.vcf.gz.tbi ]\n";
		    print MANTA_SH "then\n";
		    print MANTA_SH "\ttouch $manta_log_dir/SV_MANTA_$sample_tumor_name.done\n";
		    print MANTA_SH "fi\n";
    
		    my $qsub = &qsubTemplate(\%opt,"MANTA");
		    if (@running_jobs){
			system "$qsub -o $manta_log_dir/$job_id.out -e $manta_log_dir/$job_id.err -N $job_id -hold_jid ".join(",",@running_jobs)." $bashFile";
		    } else {
			system "$qsub -o $manta_log_dir/$job_id.out -e $manta_log_dir/$job_id.err -N $job_id $bashFile";
		    }
		    push(@manta_jobs, $job_id);
		}
	    }
	}
    }
    return(\@manta_jobs);
}

sub runDelly {
    ## Former globals?? # Probably breaks code....
    my @runningJobs;
    my @sampleBams;
    my ($delly_out_dir, $delly_log_dir, $delly_job_dir, $delly_tmp_dir);
    my %opt;

    ###
    # Run structural variant caller Delly
    ###
    my $configuration = shift;
    %opt = %{$configuration};
    my $runName = (split("/", $opt{OUTPUT_DIR}))[-1];
    my $jobID = "SV_DELLY_".get_job_id();
    my $qsub = &qsubTemplate(\%opt,"DELLY_MERGE");


    # Skip Delly if .done file exists
    if (-e "$opt{OUTPUT_DIR}/logs/SV_DELLY.done"){
	print "WARNING: $opt{OUTPUT_DIR}/logs/SV_DELLY.done exists, skipping \n";
	return \%opt;
    }

    # Create output, log, job and tmp directories.
    $delly_out_dir = "$opt{OUTPUT_DIR}/structuralVariants/delly/";
    $delly_log_dir = "$delly_out_dir/logs/";
    $delly_job_dir = "$delly_out_dir/jobs/";
    $delly_tmp_dir = "$delly_out_dir/tmp/";
    if (! -e $delly_log_dir) {
	make_path($delly_log_dir) or die "Couldn't create directory: $delly_log_dir\n $! \n";
    }
    if (! -e $delly_job_dir) {
	make_path($delly_job_dir) or die "Couldn't create directory: $delly_job_dir\n $! \n";
    }
    if (! -e $delly_tmp_dir) {
	make_path($delly_tmp_dir) or die "Couldn't create directory: $delly_tmp_dir\n $! \n";
    }

    # Get sample bam files and store running jobs
    foreach my $sample (@{$opt{SAMPLES}}){
	my $sampleBam = "$opt{OUTPUT_DIR}/$sample/mapping/$opt{BAM_FILES}->{$sample}";
	push @sampleBams, $sampleBam;
	if ( @{$opt{RUNNING_JOBS}->{$sample}} ){
    	    push( @runningJobs, @{$opt{RUNNING_JOBS}->{$sample}} );
        }
    }

    # Get SV types and split settings configured in ini/config
    my @svTypes = split/\t/, $opt{DELLY_SVTYPE};
    my @svSplit = split/\t/, $opt{DELLY_SPLIT};

    my @jobIDs_concat;
    my %chrs;

    for(my $i = 0; $i <= $#svTypes; $i++) {
	my $type = $svTypes[$i];

	# Skip SV calling for type if .done exist
	if (-e "$delly_log_dir/DELLY_$type.done"){
	    print "WARNING: $delly_log_dir/DELLY_$type.done exists, skipping \n";
	    next;
	}

	# Split per chromosome
	if ($svSplit[$i] eq "yes") {
	    get_chrs_from_dict(\%chrs,\%opt) unless scalar(keys %chrs);
	    my ( $jobIDs_chunks, $logFiles );
	    ( $jobIDs_chunks, $logFiles ) = create_interchromosomal_chunks(\@sampleBams, \%chrs, $type, $delly_tmp_dir, $delly_job_dir, $delly_log_dir,\%opt, \@runningJobs) if $type eq "TRA";
	    ( $jobIDs_chunks, $logFiles ) = create_intrachromosomal_chunks(\@sampleBams, \%chrs, $type, $delly_tmp_dir, $delly_job_dir, $delly_log_dir,\%opt, \@runningJobs) if $type =~ /DEL|DUP|INV/;

	    # Translocation jobs
	    if ($type eq "TRA") {
		my $jobID = "CONVERT_".get_job_id();
		my $convert_file = "$delly_job_dir/$type\_".$jobID.".sh";
		open CONVERT, ">$convert_file";
		print CONVERT "#!/bin/bash\n\n";
		print CONVERT "FILES=(".join(" ", @$logFiles).")\n";
		print CONVERT "FINISHED=()\n";
		print CONVERT "FAILED=()\n";
		print CONVERT "for FILE in \${FILES[\@]}\n";
		print CONVERT "do\n";
		print CONVERT "\tTAIL=`tail -n 1 \$FILE`\n";
		print CONVERT "\tif [[ \$TAIL =~ Done.\$ ]] ; then\n";
		print CONVERT "\t\tDONE_FILE=\${FILE//DELLY_[a-zA-Z0-9]*.out/DELLY.done}\n";
		print CONVERT "\t\ttouch \$DONE_FILE\n";
		print CONVERT "\t\tVCF_FILE=\$FILE\n";
		print CONVERT "\t\tVCF_FILE=\${VCF_FILE//\/logs\//\/tmp\/}\n";
		print CONVERT "\t\tVCF_FILE=\${VCF_FILE//_DELLY_[a-zA-Z0-9]*.out/.vcf}\n";
		print CONVERT "\t\tif [ -f \"\$VCF_FILE\" ]\n";
		print CONVERT "\t\tthen\n";
		print CONVERT "\t\t\techo \$VCF_FILE >> $delly_tmp_dir/$type\_vcf_files.txt\n";
		print CONVERT "\t\tfi\n";
		print CONVERT "\t\tFINISHED+=(\$FILE)\n";
	        print CONVERT "\telse\n";
	        print CONVERT "\t\tFAILED+=(\$FILE)\n";
	        print CONVERT "\tfi\n";
	        print CONVERT "done\n\n";
	        print CONVERT "if [[ \${#FAILED[@]} > 0 ]] ; then\n";
	        print CONVERT "\t>&2 echo \"error\"\n";
	        print CONVERT "else\n";
	        print CONVERT "\t$opt{VCFTOOLS_PATH}/vcf-concat -f $delly_tmp_dir/$type\_vcf_files.txt | $opt{VCFTOOLS_PATH}/vcf-sort -c > $delly_tmp_dir/$runName\_$type.vcf\n";
	        print CONVERT "\t$opt{IAP_PATH}/scripts/delly_TRA_convert.pl $delly_tmp_dir/$runName\_$type.vcf\n";
	        print CONVERT "fi\n";
		close CONVERT;
		system "$qsub -o $delly_log_dir/$type\_CONVERT.out -e $delly_log_dir/$type\_CONVERT.err -N $jobID -hold_jid ".join(",",@$jobIDs_chunks). " $convert_file";
		push @$jobIDs_chunks, $jobID;
		
		my $jobID2 = "VCF_CONCAT_".get_job_id();
	        my $vcf_concat_file = "$delly_job_dir/$type\_".$jobID2.".sh";
		open VCF_CONCAT, ">$vcf_concat_file";
		print VCF_CONCAT "#!/bin/bash\n\n";
		print VCF_CONCAT "FILE=$delly_log_dir/$type\_CONVERT.out\n";
		print VCF_CONCAT "TAIL=`tail -n 1 \$FILE`\n";
		print VCF_CONCAT "if [[ \$TAIL =~ Done.\$ ]] ; then\n";
		print VCF_CONCAT "\t$opt{VCFTOOLS_PATH}/vcf-sort -c $delly_tmp_dir/$runName\_$type\_CONVERT.vcf > $delly_tmp_dir/$runName\_$type\_CONVERT_SORT.vcf\n";
		print VCF_CONCAT "\tmv $delly_tmp_dir/$runName\_$type\_CONVERT_SORT.vcf $delly_out_dir/$runName\_$type.vcf\n";
	        print VCF_CONCAT "\ttouch $delly_log_dir/DELLY_$type.done\n";
	        print VCF_CONCAT "fi\n\n";
		close VCF_CONCAT;
		
		system "$qsub -o $delly_log_dir/$type\_VCF_CONCAT.out -e $delly_log_dir/$type\_VCF_CONCAT.err -N $jobID2 -hold_jid ".join(",",@$jobIDs_chunks). " $vcf_concat_file";

		push @jobIDs_concat, $jobID2;
	    # Other sv types
	    } else {
		my $jobID = "VCF_CONCAT_".get_job_id();
	        my $vcf_concat_file = "$delly_job_dir/$type\_".$jobID.".sh";
		open VCF_CONCAT, ">$vcf_concat_file";
		print VCF_CONCAT "#!/bin/bash\n\n";
		print VCF_CONCAT "FILES=(".join(" ", @$logFiles).")\n";
		print VCF_CONCAT "FINISHED=()\n";
		print VCF_CONCAT "FAILED=()\n";
		print VCF_CONCAT "for FILE in \${FILES[\@]}\n";
		print VCF_CONCAT "do\n";
		print VCF_CONCAT "\tTAIL=`tail -n 1 \$FILE`\n";
		print VCF_CONCAT "\tif [[ \$TAIL =~ Done.\$ ]] ; then\n";
		print VCF_CONCAT "\t\tDONE_FILE=\${FILE//DELLY_[a-zA-Z0-9]*.out/DELLY.done}\n";
		print VCF_CONCAT "\t\ttouch \$DONE_FILE\n";
		print VCF_CONCAT "\t\tVCF_FILE=\$FILE\n";
		print VCF_CONCAT "\t\tVCF_FILE=\${VCF_FILE//\/logs\//\/tmp\/}\n";
		print VCF_CONCAT "\t\tVCF_FILE=\${VCF_FILE//_DELLY_[a-zA-Z0-9]*.out/.vcf}\n";
		print VCF_CONCAT "\t\tif [ -f \"\$VCF_FILE\" ]\n";
		print VCF_CONCAT "\t\tthen\n";
		print VCF_CONCAT "\t\t\techo \$VCF_FILE >> $delly_tmp_dir/$type\_vcf_files.txt\n";
		print VCF_CONCAT "\t\tfi\n";
		print VCF_CONCAT "\t\tFINISHED+=(\$FILE)\n";
	        print VCF_CONCAT "\telse\n";
	        print VCF_CONCAT "\t\tFAILED+=(\$FILE)\n";
	        print VCF_CONCAT "\tfi\n";
	        print VCF_CONCAT "done\n\n";
	        print VCF_CONCAT "if [[ \${#FAILED[@]} > 0 ]] ; then\n";
	        print VCF_CONCAT "\t>&2 echo \"error\"\n";
	        print VCF_CONCAT "else\n";
	        print VCF_CONCAT "\t$opt{VCFTOOLS_PATH}/vcf-concat -f $delly_tmp_dir/$type\_vcf_files.txt | $opt{VCFTOOLS_PATH}/vcf-sort -c > $delly_tmp_dir/$runName\_$type.vcf\n";
	        print VCF_CONCAT "\tmv $delly_tmp_dir/$runName\_$type.vcf $delly_out_dir/$runName\_$type.vcf\n";
	        print VCF_CONCAT "\ttouch $delly_log_dir/DELLY_$type.done\n";
	        print VCF_CONCAT "fi\n\n";
	        close VCF_CONCAT;
	        system "$qsub -o $delly_log_dir/$type\_VCF_CONCAT.out -e $delly_log_dir/$type\_VCF_CONCAT.err -N $jobID -hold_jid ".join(",",@$jobIDs_chunks). " $vcf_concat_file";
	        push @jobIDs_concat, $jobID;
	    }
	# Non split jobs
	} else {
	    my $jobID = "DELLY_".get_job_id();
	    my $dellyFile = "$delly_job_dir/$type\_".$jobID.".sh";
	    my $logFile = $dellyFile;
	    $logFile =~ s/jobs/logs/;
	    $logFile =~ s/.sh$/.out/;

	    submit_delly($dellyFile, $jobID, $type, "", "$delly_tmp_dir/$runName\_$type.vcf", \%opt, \@runningJobs, \@sampleBams);
	    # Translocation jobs
	    if ($type eq "TRA") {
		    my $jobID2 = "CONVERT_".get_job_id();
		    my $convert_file = "$delly_job_dir/$type\_".$jobID2.".sh";
		    open CONVERT, ">$convert_file";
		    print CONVERT "#!/bin/bash\n\n";
		    print CONVERT "TAIL=`tail -n 1 $logFile`\n";
    		    print CONVERT "if [[ \$TAIL =~ Done.\$ ]] ; then\n";
    		    print CONVERT "\t$opt{IAP_PATH}/scripts/delly_TRA_convert.pl $delly_tmp_dir/$runName\_$type.vcf\n";
    		    print CONVERT "fi\n";
    		    close CONVERT;
		    system "$qsub -o $delly_log_dir/$type\_CONVERT.out -e $delly_log_dir/$type\_CONVERT.err -N $jobID2 -hold_jid $jobID $convert_file";

		    my $jobID3 = "VCF_CONCAT_".get_job_id();
		    my $vcf_concat_file = "$delly_job_dir/$type\_".$jobID3.".sh";
		    open VCF_CONCAT, ">$vcf_concat_file";
		    print VCF_CONCAT "#!/bin/bash\n\n";
    		    print VCF_CONCAT "FILE=$delly_log_dir/$type\_CONVERT.out\n";
		    print VCF_CONCAT "TAIL=`tail -n 1 \$FILE`\n";
		    print VCF_CONCAT "if [[ \$TAIL =~ Done.\$ ]] ; then\n";
		    print VCF_CONCAT "\t$opt{VCFTOOLS_PATH}/vcf-sort -c $delly_tmp_dir/$runName\_$type\_CONVERT.vcf > $delly_tmp_dir/$runName\_$type\_CONVERT_SORT.vcf\n";
		    print VCF_CONCAT "\tmv $delly_tmp_dir/$runName\_$type\_CONVERT_SORT.vcf $delly_out_dir/$runName\_$type.vcf\n";
	    	    print VCF_CONCAT "\ttouch $delly_log_dir/DELLY_$type.done\n";
	    	    print VCF_CONCAT "fi\n\n";
		    close VCF_CONCAT;

		    system "$qsub -o $delly_log_dir/$type\_VCF_CONCAT.out -e $delly_log_dir/$type\_VCF_CONCAT.err -N $jobID3 -hold_jid $jobID2 $vcf_concat_file";
		
		    push @jobIDs_concat, $jobID3;
	    # Other sv types
	    } else {
		my $jobID2 = "VCF_CONCAT_".get_job_id();
		my $vcf_concat_file = "$delly_job_dir/$type\_".$jobID2.".sh";
		open VCF_CONCAT, ">$vcf_concat_file";
		print VCF_CONCAT "#!/bin/bash\n\n";
		print VCF_CONCAT "TAIL=`tail -n 1 $logFile`\n";
		print VCF_CONCAT "if [[ \$TAIL =~ Done.\$ ]] ; then\n";
		print VCF_CONCAT "\tmv $delly_tmp_dir/$runName\_$type.vcf $delly_out_dir/$runName\_$type.vcf\n";
		print VCF_CONCAT "\ttouch $delly_log_dir/DELLY_$type.done\n";
		print VCF_CONCAT "else\n";
		print VCF_CONCAT "\t>&2 echo \"error\"\n";
		print VCF_CONCAT "fi\n\n";
		close VCF_CONCAT;
		system "$qsub -o $delly_log_dir/$type\_VCF_CONCAT.out -e $delly_log_dir/$type\_VCF_CONCAT.err -N $jobID2 -hold_jid $jobID $vcf_concat_file";
		
		push @jobIDs_concat, $jobID2;
	    }
	}
    }
    return(\@jobIDs_concat);
}

### Submit delly jobs
sub submit_delly {
    my ($bashFile, $jobID, $type, $excludeFile, $vcfFile, $config, $jobs, $bams) = @_;
    my ($logFile, $errorFile) = ($bashFile, $bashFile);
    my %opt = %{$config};
    my @runningJobs = @{$jobs};
    my @sampleBams = @{$bams};
    $logFile =~ s/jobs/logs/;
    $logFile =~ s/.sh$/.out/;
    $errorFile =~ s/jobs/logs/;
    $errorFile =~ s/.sh$/.err/;
    my $omp_num_threads = $opt{DELLY_THREADS} * 2;
    open DELLY_SH , ">$bashFile" or die "cannot open file $bashFile\n $! \n";
    print DELLY_SH "#!/bin/bash\n\n";
    print DELLY_SH "export OMP_NUM_THREADS=".$omp_num_threads."\n";
    print DELLY_SH "$opt{DELLY_PATH}/delly";
    print DELLY_SH " -t " . $type;
    print DELLY_SH " -g " . $opt{GENOME};
    print DELLY_SH " -x " . $excludeFile if $excludeFile;
    print DELLY_SH " -q " . $opt{DELLY_MAPQUAL};
    print DELLY_SH " -s " . $opt{DELLY_MAD};
    print DELLY_SH " -m " . $opt{DELLY_FLANK};
    print DELLY_SH " -u " . $opt{DELLY_GENO_QUAL};
    print DELLY_SH " -v " . $opt{DELLY_VCF_GENO} if $opt{DELLY_VCF_GENO};
    print DELLY_SH " -o " . $vcfFile;
    print DELLY_SH " ".join(" ", @sampleBams);
    close DELLY_SH;
    my $qsub = &qsubTemplate(\%opt,"DELLY");
    if (@runningJobs) {
	system "$qsub -o $logFile -e $errorFile -N $jobID -hold_jid ".join(",",@runningJobs)." $bashFile";
    } else {
        system "$qsub -o $logFile -e $errorFile -N $jobID $bashFile";
    }
    return($logFile);
}

### Create inter chromosomal chunks
sub create_interchromosomal_chunks {
    my ($bams, $chrs, $type, $delly_tmp_dir, $delly_job_dir, $delly_log_dir, $config, $jobs) = @_;
    my @runningJobs = @{$jobs};
    my @jobIDs;
    my @logFiles;
    my %opt = %{$config};
    open VCF_FILES, ">$delly_tmp_dir/$type\_vcf_files.txt";
    foreach my $chr1 (keys %{$chrs}) {
	foreach my $chr2 (keys %{$chrs}) {
	    next unless $chr2 gt $chr1;
	    my $vcfFile = "$delly_tmp_dir/$type\_$chr1\_$chr2.vcf";
	    print VCF_FILES $vcfFile . "\n" if -e $vcfFile;
	    if (-e "$delly_log_dir/$type\_$chr1\_$chr2\_DELLY.done") {
		print "WARNING: $delly_log_dir/$type\_$chr1\_$chr2\_DELLY.done exists, skipping \n";
		next;
	    }
	    my $excludeFile = "$delly_tmp_dir/$type\_$chr1\_$chr2\_exclude.txt";
	    open EXC, ">$excludeFile";
	    foreach my $chrom (keys %{$chrs}) {
		#print EXC join("\t", $chrom, -1, $chrs->{$chrom}) . "\n" unless $chrom =~ /^$chr1$/ or $chrom =~ /^$chr2$/;
		print EXC $chrom . "\n" unless $chrom =~ /^$chr1$/ or $chrom =~ /^$chr2$/;
	    }
	    close EXC;
	    my $jobID = "DELLY_".get_job_id();
	    my $dellyFile = "$delly_job_dir/$type\_$chr1\_$chr2\_".$jobID.".sh";
	    push @jobIDs, $jobID;
	    my ( $logFile ) = submit_delly($dellyFile, $jobID, $type, $excludeFile, $vcfFile,\%opt, \@runningJobs, $bams);
	    push @logFiles, $logFile;
	}
    }
    close VCF_FILES;
    return(\@jobIDs, \@logFiles);
}

### Create intra chromosomal chunks
sub create_intrachromosomal_chunks {
    my ($bams, $chrs, $type, $delly_tmp_dir, $delly_job_dir, $delly_log_dir,$config, $jobs) = @_;
    my @runningJobs = @{$jobs};
    my @jobIDs;
    my @logFiles;
    my %opt = %{$config};
    open VCF_FILES, ">$delly_tmp_dir/$type\_vcf_files.txt";
    foreach my $chr (keys %{$chrs}) {
	my $vcfFile = "$delly_tmp_dir/$type\_$chr.vcf";
	print VCF_FILES $vcfFile . "\n" if -e $vcfFile;
	if ( -e "$delly_log_dir/$type\_$chr\_DELLY.done" ) {
	    print "WARNING: $delly_log_dir/$type\_$chr\_DELLY.done exists, skipping \n";
	    next;
	}
	my $excludeFile = "$delly_tmp_dir/$type\_$chr\_exclude.txt";
	open EXC, ">$excludeFile";
	foreach my $chrom (keys %{$chrs}) {
	    #print EXC join("\t", $chrom, -1, $chrs->{$chrom}) . "\n" unless $chrom =~ /^$chr$/;
	    print EXC $chrom . "\n" unless $chrom =~ /^$chr$/;
	}
	close EXC;
	
	my $jobID = "DELLY_".get_job_id();
	my $dellyFile = "$delly_job_dir/$type\_$chr\_".$jobID.".sh";
	push @jobIDs, $jobID;
	my ( $logFile ) = submit_delly($dellyFile, $jobID, $type, $excludeFile, $vcfFile,\%opt, \@runningJobs, $bams);
	push @logFiles, $logFile;
    }
    close VCF_FILES;
    return(\@jobIDs, \@logFiles);
}

### Get chromosomes from genome.dict
sub get_chrs_from_dict {
    my ($chrs, $config) = @_;
    my %chrs = %{$chrs};
    my %opt = %{$config};
    my $dictFile = $opt{GENOME};
    $dictFile =~ s/.fasta$/.dict/;
    open DICT, $dictFile;
    while(<DICT>) {
	chomp;
	my ($chr, $length) = ($1, $2) if $_ =~ /SN:(\w+)\s*LN:(\d+)/;
	$chrs->{$chr} = $length if $chr;
    }
    close DICT;
}

############
sub get_job_id {
    my $id = tmpnam();
    $id=~s/\/tmp\/file//;
    return $id;
}
############

1;
