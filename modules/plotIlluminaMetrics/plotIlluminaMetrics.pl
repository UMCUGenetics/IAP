#!/usr/bin/perl
### Plot Illumina Metrics
### Robert Ernst
### 15-01-2013

use warnings;
use strict;
use Cwd            qw( abs_path );
use File::Basename qw( dirname );

die usage() if @ARGV == 0;
my @samples = @ARGV;

#get dir where plotHSMetric is installed
my $rootDir = dirname(abs_path($0));
my $runDir = abs_path();
my $runName = (split("/", $runDir))[-1];

### Parse HSMetrics files if present
my @files = grep -f, <*/QCStats/*HSMetrics\.txt>;
if ( scalar(@files) != 0 ) {
    my $fileName = $runName.".HSMetric_summary.txt";
    open(SUMFILE, ">", $fileName) || die ("Can't create $fileName");
    my $printedHeader = 0;

    ### Generate metric summary
    print "Parsing HSMetric files \n";
    foreach my $file (@files) {
	print "\t Parsing: ". $file . "\n";
	open(FILE, $file) || die ("Can't open $file");
	my $baitIntervals;
	my $targetIntervals;

	#Processing headerlines -> beginning with # or empty lines.
        while(<FILE> =~ /(^\s*#)(.*)/ || <FILE> eq "") {
	    my $headerLine = $2;
	    if ($headerLine =~ m/BAIT_INTERVALS=(\S*).TARGET_INTERVALS=(\S*)/){ #grep bait and target interval file settings.
		$baitIntervals = $1;
		$targetIntervals = $2;
	    }
	}

        #Processing table
	my $tableHeader =  <FILE>;
	unless ($printedHeader) { #print table header once.
	    print SUMFILE "sample \t samplePath \t baitIntervals \t targetIntervals \t". $tableHeader; 
	    $printedHeader = 1;
	}
	
	my $line = <FILE>; #grep statistics, here we assume one row with statistics.
	my $sample = (split("/",$file))[0];
	my $samplePath = (split("_",$file))[0];
	print SUMFILE $sample ."\t". $samplePath ."\t". $baitIntervals ."\t". $targetIntervals ."\t". $line;
    }
}

### Parse WGSMetrics files if present
@files = grep -f, <*/QCStats/*WGSMetrics\.txt>;
if ( scalar(@files) != 0 ) {
    my $fileName = $runName.".WGSMetric_summary.txt";
    open(SUMFILE, ">", $fileName) || die ("Can't create $fileName");
    my $printedHeader = 0;
    
    print "Parsing WGSMetric files \n";
    foreach my $file (@files) {
	print "\t Parsing: ". $file . "\n";
	open(FILE, $file) || die ("Can't open $file");
    
    #Processing headerlines -> beginning with # or empty lines.
        while(<FILE> =~ /(^\s*#)(.*)/ || <FILE> eq "") {}
        
        my $tableHeader =  <FILE>;
	unless ($printedHeader) { #print table header once.
	    print SUMFILE "sample \t samplePath \t". $tableHeader; 
	    $printedHeader = 1;
	}
	
	my $line = <FILE>; #grep statistics, here we assume one row with statistics.
	my $sample = (split("/",$file))[0];
	my $samplePath = (split("_",$file))[0];
	print SUMFILE $sample ."\t". $samplePath ."\t". $line;
    }
}


### Run R plot script and markdown to generate pdf
system "Rscript $rootDir/plotIlluminaMetrics.R ".$rootDir." ". $runName." ".join(" ",@samples);
system "rm plotIlluminaMetrics.md";
system "rm -r pdfFigures";
system "rm plotIlluminaMetrics.tex"; #clean up tex to pdf conversion.
system "mv plotIlluminaMetrics.pdf ".$runName.".picardMetrics.pdf";


sub usage{
    warn <<END;
    Usage: perl plotIlluminaMetrics.pl [sample_1 sample_2 sample_n]
END
    exit;
}