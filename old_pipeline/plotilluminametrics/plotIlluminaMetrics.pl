#!/usr/bin/perl
### Plot Illumina Metrics
### Robert Ernst
### 15-01-2013

use warnings;
use strict;

#get dir where plotHSMetric is installed
use Cwd            qw( abs_path );
use File::Basename qw( dirname );
my $rootDir = dirname(abs_path($0));
my $runDir = abs_path();
my $runName = (split("/", $runDir))[-1];

### Setup variables
my $fileName = $runName.".HSMetric_summary.txt";
open(SUMFILE, ">", $fileName) || die ("Can't create $fileName");
my @files = grep -f, <*/*HSMetrics\.txt>;
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
        print SUMFILE "sample \t sampleShort \t baitIntervals \t targetIntervals \t". $tableHeader; 
        $printedHeader = 1;
    }
    
    my $line = <FILE>; #grep statistics, here we assume one row with statistics.
    my $sample = (split("/",$file))[0];
    my $sampleShort = (split("_",$file))[0];
    print SUMFILE $sample ."\t". $sampleShort ."\t". $baitIntervals ."\t". $targetIntervals ."\t". $line;
}

### Run R plot script and markdown to generate pdf
`Rscript $rootDir/plotIlluminaMetrics.R $fileName $rootDir $runName`;
`rm plotIlluminaMetrics.md`;
`rm -r pdfFigures`;
`rm plotIlluminaMetrics.tex`; #clean up tex to pdf conversion.
`mv plotIlluminaMetrics.pdf $runName".picardMetrics.pdf"`;