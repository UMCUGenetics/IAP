#!/usr/bin/perl -w
use strict;
use POSIX qw(tmpnam);
#usage script OUTFILE infile1 infile2 ...
#my $outfile = shift(@ARGV);
#my $infiles = join(" ", @ARGV);

my $pwd = `pwd`; 
chomp($pwd);
chdir ($pwd);


foreach my $sample (<*>) {
    next unless -d $sample; 
    chdir $sample;
    


    my @infiles;
    my $outfile;
    my $merged = 0;
    my $command = '';

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

	if (-e $outfile) {
	    warn "cannot guess correct outputfile name or file exists: $!...\n\n";
	    chdir '../';
	}



	print "merging ", scalar @infiles, " files into $outfile... \n";

	my $infiles = join(" ",@infiles);
	#system "sambamba_v0.4.3 merge -t 12 -p $outfile $infiles";
	$command = "sambamba_v0.4.3 merge -t 12 $outfile $infiles\n";
	print "\nindexing...\n";
	#system "sambamba_v0.4.3 index -t 12 -p $outfile > $outfile.bai";
	$command .= "sambamba_v0.4.3 index -t 12 $outfile > $outfile.bai\n";

	my $flg = $outfile;
	$flg =~s/bam/flagstat/;
	print "\nflagstat...\n";
	#system "sambamba_v0.4.3 flagstat -t 12 -p $outfile > $flg";
	$command .= "sambamba_v0.4.3 flagstat -t 12 $outfile > $flg\n";

    }else{
	print "nothing to do...\n";
	chdir '../';
	next;
    }
    
    cluster($command,"$pwd/$sample");
    print "\ndone!\n\n"; 
    
    chdir '../';  
}



sub cluster {
    my $jid2 = tmpnam();
    $jid2=~s/\/tmp\/file//;
    my $comm = shift;
    my $pwd = shift;
    
    open OUT, ">$pwd/MERGE_$jid2.sh" or die "cannot open file $pwd/MERGE_$jid2.sh\n\n";
    print OUT "#!/bin/bash\n\n";
    print OUT "cd $pwd\n";
    print OUT "$comm\n";
    
    system "qsub -q veryshort -pe threaded 8 -o $pwd/reads -e $pwd/reads -N MERGE_$jid2 $pwd/MERGE_$jid2.sh";
    return $jid2;


}