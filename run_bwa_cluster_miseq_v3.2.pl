#!/usr/bin/perl -w
use strict;
use POSIX qw(tmpnam);
use Getopt::Long;
use Data::Dumper;

#################################################################################
#										#
#	Mi/Hiseq mapper pipeline	Ies Nijman	sept 2011		#
#										#
#	v3.2	uses sambamba for merging, indexing and flagstat		#
#										#
#################################################################################
my %opt;

BEGIN{ 
    %opt = (
		'help'          => undef, 
		'cluster_path'  => "/opt/sge/default/common/",
		'bwa_path'      => "/hpc/cog_bioinf/common_scripts/bwa-0.7.5a",
		'gatk_path'     => "/hpc/cog_bioinf/common_scripts/GenomeAnalysisTK-2.8-1/",
		'picard_path'   => "/hpc/cog_bioinf/common_scripts/picard-tools-1.98",
		'samtools_path' => "/hpc/cog_bioinf/common_scripts/samtools",
		'fastqc_path' 	=> "/hpc/cog_bioinf/common_scripts/FastQC",
		'sambamba_path' => "/hpc/cog_bioinf/common_scripts/sambamba/sambamba_v0.4.3",
		'qualimap_path' => "/hpc/cog_bioinf/common_scripts/qualimap-build-11-11-13",
		'species'       => 'HUMAN',
		'genome'	=> "/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta",
		'design'	=> "",
		'threads'	=> 10,
		'pre_stats'	=> "yes",
		'post_stats'	=> "no",
		'queue'		=> "veryshort",
		'email'		=> "i.nijman\@umcutrecht.nl",
		'mark_dups'	=> "yes",
    );
    
    
    sub usage{ 

	print <<END;

 Usage: 
 
  Required params:
   -genome		[s]	reference genome [$opt{genome}]
 
  Options:
   -cluster_path	[s]	Path to clusterfiles [$opt{cluster_path}]
   -bwa_path		[s]	path to bwa [$opt{bwa_path}]
   -gatk_path		[s]	path to gatk [$opt{gatk_path}]
   -picard_path		[s]	path to picard [$opt{picard_path}]
   -samtools_path	[s]	path to samtools [$opt{samtools_path}]
   -fastqc_path		[s]	path to FastQc [$opt{fastqc_path}]
   -sambamba_path	[s]	path to sambamba [$opt{sambamba_path}]
   -qualimap_path	[s]	path to Qualimap [$opt{qualimap_path}]
   -species		[s]	species [$opt{species}]
   -design		[s]	design for enrichment [$opt{design}]
   -mark_duplicates	[s]	mark duplicates? [$opt{mark_dups}]
   -threads		[i]	number of threads [$opt{threads}]
   -queue		[s]	cluster queue [$opt{queue}]
   -email		[s]	Email adres [$opt{email}]
   -pre_stats		[s]	read statistics FastQC [$opt{pre_stats}]
   -post_stats		[s]	mapping statistics [$opt{post_stats}]
 
   -h/help		[s]	help
  
  
  Description:

  Notes:

	  
END
	exit;
    }

    #die usage() if @ARGV == 0;
    GetOptions (
		'h|help'      => \$opt{help},
		'cluster_path=s'	=> \$opt{cluster_path},
		'bwa_path=s'		=> \$opt{bwa_path},
		'gatk_path=s'		=> \$opt{gatk_path},
		'species=s'		=> \$opt{species},
		'samtools_path=s'	=> \$opt{samtools_path},
		'picard_path=s'		=> \$opt{picard_path},
		'fastqc_path=s'		=> \$opt{fastqc_path},
		'sambamba_path=s'	=> \$opt{sambamba_path},
		'qualimap_path=s'	=> \$opt{qualimap_path},
		'genome=s'		=> \$opt{genome},
		'design=s'		=> \$opt{design},
		'mark_duplicates=s'	=> \$opt{mark_dups},
		'threads|t=i'		=> \$opt{threads},
		'queue|q=s'		=> \$opt{queue},
		'email=s'		=> \$opt{email},
		'pre_stats=s'		=> \$opt{pre_stats},
		'post_stats=s'		=> \$opt{post_stats},
    ) 
    or die usage();
    die usage() if $opt{help};
}

#my $GENOME = $ARGV[0];
my $FAI = "$opt{genome}\.fai";
die "genome $opt{genome} does not exists!!\t$!\n" if !-e "$opt{genome}.bwt";
die "fai file $FAI does not exists!!\n" if !-e $FAI;


my @to_run;
while (my $f=<reads/*_R1*.fastq.gz>) {
  #ssc11Bcell_S6_L001_R1_001.fastq.gz
  
  $f=~s/reads\///;
  push @to_run, $f;
}


if (scalar @to_run ==0) {
    print "Nothing to do!!\n\n";
    usage;
}

my $workdir = `pwd`; 
$workdir=~s/\n//;

open QSUB,">qsub2.sh";
print QSUB "\#!/bin/sh\n\n. $opt{cluster_path}/settings.sh\n\n";

mkdir "results" unless -e "results";

my (@jobs_to_wait, @bamfiles);
foreach my $p (@to_run) {
	my $coreName = $p;
	print $p, " found\n";
	$coreName =~ s/\.fastq.gz//;
	$coreName =~ s/\_R1//;
	$coreName =~ s/\_R2//;
	my $sam = "$coreName.sam";
	my $bam = "$coreName.bam";
	
	print "corename: $coreName\n";
	
	#create part ids
	my $S1 = $p;
	my $S2 = $p;
	$S2 =~ s/R1/R2/;
	
	print $S1," = S1\n";
	print $S2," = S2\n";
	print $coreName," = corename\n";
	
	#detect if run was paired of single fragment
	my $paired = 0;
	if (-e "reads/$S2") {
	    $paired = 1;
	    print "Switching to paired end mode!\n";
	}else{
	    print "single tag found, switching to fragment mode\n";
	}
	
	#ssc11Bcell_S6_L001_R1_001.fastq.gz
	#CONTROLP25_H7YRLADXX_ATCACG_L001_001
	##trying to parse readgroup parameters from readnames.
	$coreName =~ s/_\d{3}$//;
	my ($RG_PL, $RG_ID, $RG_LB, $RG_SM, $RG_PU) = ('ILLUMINA', $coreName, $coreName, $coreName, '');
	
	#get readgroup ID from readname
	my @fields = split(/\_/, $coreName);
	($RG_SM, $RG_PU, $RG_LB)  = ($fields[0], $fields[1], $fields[2]);
	#($RG_SM, $RG_LB, $RG_ID) = ($1,$2,$3) if ($coreName =~ /^([\w\d]+?)_([\w\d]+?)_(L\d+)_\d+/);

	print "readgroups deduced from $coreName: RG_PL: $RG_PL, RG_ID: $RG_ID, RG_LB: $RG_LB, RG_SM: $RG_SM, RG_PU: $RG_PU\n";

	if (-e "results/$coreName\_sorted.bam"){
		#warn "Skipping $p\n";
		#next;
	}else{
	    print "using $p\n";
	}
	
	
	print "$p >> RG $RG_ID / LIB $RG_LB / SAMPLE $RG_SM / PLATFORM $RG_PL\n";

	my $job_id = "BWA$p".get_job_id();
	my ($prestats1_job_id, $prestats2_job_id) = ("preStats1".get_job_id(), "preStats2".get_job_id()); 
	
	warn "writing submission script for $p\n";
	#my $nreads = `grep -c '>' reads/$p`;
	#$nreads=~s/\s//g;    
	open SH,">reads/$job_id.sh";
	print SH "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
	print SH "cd $workdir\n\n";
	print SH "uname -n > results/$coreName.host\n";
	
	if ($opt{pre_stats} eq "yes") {
	    open PS,">reads/$prestats1_job_id.sh";
	    print PS "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
	    print PS "cd $workdir\n\n";
	    print PS "uname -n > results/preStats_$coreName.host\n";
	    if ($paired == 1) {
		open PS2,">reads/$prestats2_job_id.sh";
		print PS2 "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
		print PS2 "cd $workdir\n\n";
		print PS2 "uname -n > results/preStats2_$coreName.host\n";
	    }
	}
	
	if ($paired == 1) {
	    if ($opt{pre_stats} eq "yes") {
		print PS "echo \"FastQC started\t\" `date` >> results/preStats_$coreName.host\n";
		print PS "$opt{fastqc_path}/fastqc reads/$S1 2>&1 >>results/prestats1_log\n";
		print PS "mv reads/*.zip ./\n";
		print PS "rm -r reads/*fastqc\n";
		print PS "echo \"FastQC finished\t\" `date` >> results/preStats_$coreName.host\n";
		print PS2 "echo \"FastQC started\t\" `date` >> results/preStats2_$coreName.host\n";
		print PS2 "$opt{fastqc_path}/fastqc reads/$S2 2>&1 >>results/prestats2_log\n";
		print PS2 "mv reads/*.zip ./\n";
		print PS2 "rm -r reads/*fastqc\n";
		print PS2 "echo \"FastQC finished\t\" `date` >> results/preStats2_$coreName.host\n";
		print QSUB "qsub -q $opt{queue} -o $workdir/reads -e $workdir/reads -N $prestats1_job_id $workdir/reads/$prestats1_job_id.sh\n";
		print QSUB "qsub -q $opt{queue} -o $workdir/reads -e $workdir/reads -N $prestats2_job_id $workdir/reads/$prestats2_job_id.sh\n";
	    }
	    print SH "echo \"mapping pair started\t\" `date` >> results/$coreName.host\n";
	    print SH "$opt{bwa_path}/bwa mem -t $opt{threads} -c 100 -M -R \"\@RG\tID:$RG_ID\tSM:$RG_SM\tPL:$RG_PL\tLB:$RG_LB\tPU:$RG_PU\" $opt{genome} reads/$S1 reads/$S2 2>>results/bwa_log | samtools view -b -S -o results/$coreName.bam -\n";
	    print SH "echo \"mapping pair finished\t\" `date` >> results/$coreName.host\n";
	}else{
	    if ($opt{pre_stats} eq "yes") {
		print PS "echo \"FastQC started\t\" `date` >> results/preStats_$coreName.host\n";
		print PS "$opt{fastqc_path}/fastqc reads/$S1 2>&1 >>results/prestats_log\n";
		print PS "mv reads/*.zip ./\n";
		print PS "rm -r reads/*fastqc\n";
		print PS "echo \"FastQC finished\t\" `date` >> results/preStats_$coreName.host\n";
		print QSUB "qsub -q $opt{queue} -o $workdir/reads -e $workdir/reads -N $prestats1_job_id $workdir/reads/$prestats1_job_id.sh\n";
	    }
	    print SH "echo \"mapping F3\t\" `date` >> results/$coreName.host\n";
	    print SH "$opt{bwa_path}\/bwa aln -t $opt{threads} $opt{genome} reads/$S1 -f results/$S1.sai 2> results/$S1.err\n";
	    print SH "echo \"starting samse\t\" `date` >> results/$coreName.host\n";
	    print SH "$opt{bwa_path}\/bwa samse -r \"\@RG\tID:$RG_ID\tSM:$RG_SM\tPL:$RG_PL\tLB:$RG_LB\tPU:$RG_PU\" $opt{genome} results/$S1.sai reads/$S1 | samtools view -b -S -o results/$coreName.bam -\n";
	    print SH "echo \"samse finished\t\" `date` >> results/$coreName.host\n";
	}
	
	print SH "echo \"Indexing\t\" `date` >> results/$coreName.host\n";
	print SH "$opt{sambamba_path} index -t $opt{threads} results/$coreName.bam\n";
	print SH "echo \"Indexing finished\t\" `date` >> results/$coreName.host\n";
	
	print SH "echo \"flagstat started\t\" `date` >> results/$coreName.host\n";
	print SH "$opt{sambamba_path} flagstat -t $opt{threads} results/$coreName.bam > results/$coreName.flagstat\n";
	print SH "echo \"flagstat finished\t\" `date` >> results/$coreName.host\n";
	
        
        ## take care: only samtools 0.19 supports this!!
        print SH "echo \"sortbam started\t\" `date` >> results/$coreName.host\n";
        print SH "$opt{samtools_path}-0.1.19/samtools sort -m 3G -@ $opt{threads} results/$coreName.bam results/$coreName\_sorted\n";		#bam ext automatically added by samtools
        print SH "$opt{samtools_path}-0.1.19/samtools index results/$coreName\_sorted.bam results/$coreName\_sorted.bai\n";
	print SH "echo \"sortbam finished\t\" `date` >> results/$coreName.host\n";
	
	#print SH "echo \"sambamba sort started\t\" `date` >> results/$coreName.host\n";
	#print SH "sleep 10\n"; # allow copying to complete
        #print SH "$opt{sambamba_path} sort -m 20GB -t $opt{threads} --tmpdir reads results/$coreName.bam -o results/$coreName\_sorted.bam\n";	
        #print SH "$opt{sambamba_path} index -t $opt{threads} results/$coreName.bam\n";
	#print SH "echo \"sambamba sort finished\t\" `date` >> results/$coreName.host\n";
	
	
	
	print SH "rm results/$coreName.bam\n";

	print SH "sleep 10\n"; # allow merging and copying to complete
	print SH "echo \"flagstat started\t\" `date` >> results/$coreName.host\n";
	print SH "$opt{sambamba_path} flagstat -t $opt{threads} results/$coreName\_sorted.bam > results/$coreName\_sorted.flagstat\n";
	print SH "echo \"flagstat finished\t\" `date` >> results/$coreName.host\n";
	print SH "date >> results/$coreName.host\n";
	print SH "mv results/*sorted.* ./\n";
	$bam = "$coreName\_sorted.bam";

	if ($opt{mark_dups} eq 'yes') {
	    #MarkDuplicates
	    #print SH "java -Xmx28G -jar $opt{picard_path}/MarkDuplicates.jar I=$coreName\_sorted.bam O=$coreName\_sorted_dedup.bam REMOVE_DUPLICATES=false M=$coreName\_sorted.bam\_picard_rmdup_metrics VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true CREATE_INDEX=true MAX_RECORDS_IN_RAM=1000000000 2>&1 >>results/dedup_log\n";
	    print SH "echo \"Markdups started\t\" `date` >> results/$coreName.host\n";
	    print SH "$opt{sambamba_path} markdup -t $opt{threads} $coreName\_sorted.bam $coreName\_sorted_dedup.bam\n";
	    print SH "echo \"Markdups finished\t\" `date` >> results/$coreName.host\n";
	    print SH "echo \"Indexing\t\" `date` >> results/$coreName.host\n";
	    print SH "$opt{sambamba_path} index -t $opt{threads} $coreName\_sorted_dedup.bam\n";
	    print SH "echo \"Indexing finished\t\" `date` >> results/$coreName.host\n";
	    print SH "sleep 30\n"; # allow merging and copying to complete
	    print SH "echo \"Flagstat started\t\" `date` >> results/$coreName.host\n";
	    print SH "$opt{sambamba_path} flagstat -t $opt{threads} $coreName\_sorted_dedup.bam > $coreName\_sorted_dedup.flagstat\n";
	    print SH "echo \"Flagstat finished\t\" `date` >> results/$coreName.host\n";
	    $bam = "$coreName\_sorted_dedup.bam";
	    push (@bamfiles, $bam);
	}




	if ($opt{post_stats} eq "yes") {
	    my ($poststats_job_id) = ("postStats".get_job_id()); 
	    open POS,">reads/$poststats_job_id.sh";
	    print POS "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
	    print POS "cd $workdir\n\n";
	    print POS "uname -n > results/postStats_$coreName.host\n";
	    my $mem = $opt{threads} * 3;
	    $mem .= 'G';
	    
	    
	    print POS "echo \"Qualimap started\t\" `date` >> results/$coreName.host\n";
	    print POS "$opt{qualimap_path}/qualimap bamqc -bam $bam -gd $opt{species} --outdir $workdir/$coreName\_qualimap --outfile $workdir/$coreName\_qualimap/$coreName\_report.pdf--java-mem-size=$mem -outformat PDF -gff $opt{design} -os -c -nt $opt{threads} 2>&1 >>results/poststats_log\n";
	    print POS "echo \"Qualimap finished\t\" `date` >> results/$coreName.host\n";
	    print SH "qsub -q veryshort -R y -M $opt{email} -m as -pe threaded $opt{threads} -o $workdir/reads -e $workdir/reads -N POST$poststats_job_id $workdir/reads/$poststats_job_id.sh\n";
	}


	close SH;
		
    print QSUB "qsub -q $opt{queue} -R y -M $opt{email} -m as -pe threaded $opt{threads} -o $workdir/reads -e $workdir/reads -N $job_id $workdir/reads/$job_id.sh\n";
    print QSUB "\n";
    push @jobs_to_wait, $job_id;
}


#submit mapping
close QSUB;
system "sh qsub2.sh";




#  merging job======================================================================================================================================
my $merge_jid = get_job_id();
open M, ">MERGE_$merge_jid.sh" or die "cannot open merge scriptfile\n";
print M "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
print M "cd $workdir\n\n";
my $outfile = $bamfiles[0];
$outfile =~ s/L\d{3}/MERGED/;


print M "uname -n > results/merge_$outfile.host\n";
print M "echo \"Merging started\t\" `date` >> results/merge.host\n";


#die "cannot guess correct outputfile name or file exists...\n\n" if -e $outfile;
print "merging ", scalar @bamfiles, " files into $outfile... \n";
my $infiles = join(" ",@bamfiles);
print M "$opt{sambamba_path} merge -t $opt{threads} -p $outfile $infiles\n";
print M "$opt{sambamba_path} index -t $opt{threads} -p $outfile > $outfile.bai\n";
my $flg = $outfile;
$flg =~s/bam/flagstat/;
print M "$opt{sambamba_path} flagstat -t $opt{threads} -p $outfile > $flg\n";
print M "echo \"Merging finished\t\" `date` >> results/merge.host\n";
foreach my $bam (@bamfiles) {
    $bam =~ s/bam/bam.bai/;
    $infiles .= " $bam ";
}
    
print M "if [ -f $outfile ]\nthen\n   rm *sorted.ba* *sorted.flagstat $infiles\nfi\n";
system "qsub -q $opt{queue} -R y -M $opt{email} -m as -pe threaded $opt{threads} -o $workdir/reads -e $workdir/reads -N MERGE_$merge_jid -hold_jid ". join(",", @jobs_to_wait)." $workdir/MERGE_$merge_jid.sh \n";
#=========================================================================================================================================================


# indelreaglinemtn======================================================================================================================================
my $scalaScript="/hpc/cog_bioinf/common_scripts/GATK_v2.7/ies/preProcessing_indelRealign.scala";
my $inputKnownIndelFile="/hpc/cog_bioinf/common_scripts/GATK_v2.7/bundle/Mills_and_1000G_gold_standard.indels.b37.vcf";
my $inputReference="/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta";
my $outputCoreName = $outfile;
$outputCoreName =~s/.bam//;
my $realign_jid = get_job_id();

open REA, ">reads/REALIGN_$realign_jid.sh" or die "cannot open realign.sh:$!\n";
print REA "#!/bin/bash\n\n";
print REA "cd $workdir\n";
print REA "echo \"Realign started\t\" `date` >> results/merge.host\n";
print REA "java -Xmx5G -Xms2G -jar /hpc/cog_bioinf/common_scripts/GATK_v2.7/Queue.jar -jobRunner GridEngine -jobQueue short -jobEnv \"threaded 1\" -S $scalaScript -I $workdir/$outfile -R $inputReference -known $inputKnownIndelFile -O $workdir/$outputCoreName -run\n";
print REA "echo \"Realign finished\t\" `date` >> results/merge.host\n";
close REA;
	
system "qsub -q short -o $workdir/reads/realign_out -e $workdir/reads/realign_err -N realign_$realign_jid -hold_jid MERGE_$merge_jid $workdir/reads/REALIGN_$realign_jid.sh";
#=========================================================================================================================================================


#=statsjob================================================================================================================================================
my $ps_jid = get_job_id();
open PS, ">reads/PICARD_$ps_jid.sh" or die "cannot open realign.sh:$!\n";
print PS "#!/bin/bash\n\n";
print PS "cd $workdir\n";
print PS "echo \"Picardstats started\t\" `date` >> results/merge.host\n";
print PS "/hpc/cog_bioinf/data/ies/scripts/illumina_proc/run_picard_metrics_single.pl $outfile \n";
print PS "echo \"Picardstats finished\t\" `date` >> results/merge.host\n";
close PS;

system "qsub -q veryshort -o $workdir/reads/picard_out -e $workdir/reads/picard_err -N picard_sub_$ps_jid  -hold_jid MERGE_$merge_jid $workdir/reads/PICARD_$ps_jid.sh";

#=========================================================================================================================================================


############
sub get_job_id {
   my $id = tmpnam(); 
      $id=~s/\/tmp\/file//;
   return $id;
}
############
