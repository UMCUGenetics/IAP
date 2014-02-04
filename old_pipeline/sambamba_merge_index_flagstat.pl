#!/usr/bin/perl -w
use strict;

#usage script OUTFILE infile1 infile2 ...
#my $outfile = shift(@ARGV);
#my $infiles = join(" ", @ARGV);

system "cd \$PWD";

my @infiles;
my $outfile;
my $merged = 0;

foreach my $file (<*_sorted_dedup.bam>) {
    if ($file =~ /MERGED/) {	#merging already done.
	print "Merging already done; outfile set to $file...\n";
	$merged=1;
	$outfile = $file;
	next;
    }
    push @infiles, $file;
    print "using $file...\n";
}





if ($merged==0) {
    $outfile = $infiles[0];
    $outfile =~ s/L\d{3}/MERGED/;

    die "cannot guess correct outputfile name...\n\n" if -e $outfile;



    print "merging ", scalar @infiles, " files into $outfile... \n";

    my $infiles = join(" ",@infiles);
    system "sambamba_v0.4.3 merge -t 12 -p $outfile $infiles";
    print "\nindexing...\n";
    system "sambamba_v0.4.3 index -t 12 -p $outfile > $outfile.bai";

    my $flg = $outfile;
    $flg =~s/bam/flagstat/;
    print "\nflagstat...\n";
    system "sambamba_v0.4.3 flagstat -t 12 -p $outfile > $flg";

}

if (!-e "$outfile\_QMreport.pdf") {
    print "Qualimap...\n";
    #system "/hpc/cog_bioinf/common_scripts/qualimap-build-11-11-13/qualimap bamqc --java-mem-size=30G -outformat PDF --outfile $outfile\_QMreport.pdf -c -nt 12 -gd HUMAN -gff /hpc/cog_bioinf/ENRICH/SS_exome_v5_S04380110_Covered_qualimap.bed -bam \$PWD/$outfile --outdir \$PWD";
}
print "\ndone!\n\n"; 


