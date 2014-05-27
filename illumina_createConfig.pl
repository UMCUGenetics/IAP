#!/usr/bin/perl -w
##################################################################################################################################################
###This script is designed to generate a .config based on user input
###
###
###Author: R.F.Ernst
###Latest change: Skeleton
###
###TODO: 
##################################################################################################################################################
use strict;
use Cwd            qw( abs_path );
use File::Basename qw( dirname );
use Getopt::Long;
use File::Path qw(make_path);

### Variables ###
my $settingsDir = dirname(abs_path($0))."/settings";

### Check usage ###
interactive() if @ARGV == 0;

### Parse options ###
my $iniFile;
my $outputDir;
my @rawDataDirs;
my $mail;
my $help;
my $run;

GetOptions ("iniFile|i=s" => \$iniFile,
	    "outputDir|o=s" => \$outputDir,
	    "rawDataDir|r=s" => \@rawDataDirs,
	    "mail|m=s" => \$mail,
	    "help|h" => \$help,
	    "run" => \$run)
or die usage();

if ($help || ! $iniFile || ! $outputDir || ! @rawDataDirs || ! $mail ) { usage() }

### Non interactive mode ###
$iniFile = $settingsDir."/".$iniFile;
createConfig($iniFile,$outputDir,\@rawDataDirs,$mail,$run);

### Interactive mode ###
sub interactive{
    print "Using interactive mode \n";
    print "Avaible setting files:\n";
    my @iniFiles = getIniFiles($settingsDir);
    
    # Settings file
    print "Choose setting file [index]: ";
    chomp(my $iniIndex = <STDIN>);
    if ($iniIndex eq "" || ! $iniFiles[$iniIndex] ) { die "Please provide a correct ini index number." }
    my $iniFile = $settingsDir ."/". $iniFiles[$iniIndex];
    
    #Output dir
    print "Output dir: ";
    chomp($outputDir = <STDIN>); # no tab completion
    if ($outputDir eq "") { die "Please provide a correct output directory." }
    
    #Raw data dir -> add while loop to allow for multiple raw data dirs.
    print "Raw data project dir: ";
    chomp(my $rawDataDir = <STDIN>); # no tab completion
    if($rawDataDir eq "") { die "Please provide a correct raw data directory." } #check for existence
    push(@rawDataDirs, $rawDataDir);
    
    #Output dir
    print "Mail address: ";
    chomp($mail = <STDIN>);
    if ($mail eq "") { die "Please provide a correct mail address." }
    
    #Create config
    createConfig($iniFile,$outputDir,\@rawDataDirs,$mail,$run);
}

### Parse and print available ini files ###
sub getIniFiles{
    my $iniDir = shift;
    my @iniFiles;
    my $iniIndex = -1;

    opendir (INIDIR, $iniDir) or die "Can't open $iniDir";
    while (my $iniFile = readdir(INIDIR)) {
	next unless ($iniFile =~ /\.ini$/); #skip non .ini files
	push(@iniFiles, $iniFile);
	$iniIndex ++;
	print "\t$iniIndex: \t $iniFile\n";
    }
    closedir(INIDIR);

    return(@iniFiles);
}

### Create config file ###
sub createConfig {
    my $iniFile = $_[0];
    my $outputDir = $_[1];
    my @rawDataDirs = @{$_[2]};
    my $mail = $_[3];
    my $run = $_[4];

    my $configFile = $outputDir."/settings.config";

    if(! -e $outputDir){
	make_path($outputDir) or die "Couldn't create directory: $outputDir\n";
    }
    # Create settings.config file in outputDir
    open CONFIG, ">$configFile" or die "cannot open file $configFile \n";
    print CONFIG "### SETTINGS ###\n";
    print CONFIG "INIFILE\t$iniFile\n";
    print CONFIG "OUTPUT_DIR\t$outputDir\n";
    print CONFIG "MAIL\t$mail\n";

    print CONFIG "\n### FASTQ FILES ###";
    #Find fastq files for each rawDataDir
    foreach my $rawDataDir (@rawDataDirs){
	if(! -e $rawDataDir) { die "$rawDataDir does not exist." }
	print CONFIG "\n# $rawDataDir\n";
	my @files = glob($rawDataDir."/*{/,}*.fastq.gz");
	foreach my $file (@files){ print CONFIG "FASTQ\t$file\n" }
    }

    close CONFIG;
    
    ###Run pipeline if -run is specified 
    if($run) {
	### run pipeline
	my $pipeline = dirname(abs_path($0))."/illumina_pipeline.pl";
	system "perl $pipeline $configFile";
    }
}

### Help information ###
sub usage{
    print "Usage: perl illumina_createConfig.pl\n\n";
    print "Advanced usage: \n";
    print "perl illumina_createConfig.pl -i|-iniFile settings.ini -o|-outputDir /path/to/outputDir -r|-rawDataDir /hiseq/140305_D00267_0081_AH8DB2ADXX/Unaligned/Project_1/ -r|-rawDataDir /hiseq/140305_D00267_0081_AH8DB2ADXX/Unaligned/Project_2 -m|-mail example\@mail.nl [-run]\n\n";
    print "Available ini files:\n";
    getIniFiles($settingsDir);
    exit;
}

