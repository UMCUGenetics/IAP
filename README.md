## IAP
Illumina variant calling pipeline. 

## Download
Seperate releases can be downloaded here: https://github.com/CuppenResearch/IAP/releases or use git clone:
```bash
git clone git@github.com:CuppenResearch/IAP.git
```

## Usage
IAP is configured using ini files and on run/analysis level using a config file. The idea is to have one ini file per run/analysis type (e.g. exome sequencing). Every setting can be reconfigured in the run/analysis config file. All ini files are located in the [settings subfolder](https://github.com/CuppenResearch/IAP/tree/master/settings). The run/analysis config is created using the illumina_createConfig script and is stored in the ouput directory.
#### View available ini files
```bash
perl illumina_createConfig.pl
```
#### Create config
```bash
perl illumina_createConfig.pl -i <filename.ini> -o </path/to/output_dir> (-f /path/to/fastq_dir OR -b /path/to/bam_dir OR -v /path/to/vcfFile.vcf) -m your@mail.com
```
#### Run pipeline
```bash
perl illumina_pipeline.pl /path/to/output_dir/settings.config>
```

## Dependencies
#### Core tools
- Opengrid engine
- Perl/dev
- Python/dev
- R/dev
- Java 1.7/jre/dev
 
#### Bio tools
- [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [BWA](http://bio-bwa.sourceforge.net/)
- [Sambamba](http://lomereiter.github.io/sambamba/)
- [bamMetrics](https://github.com/CuppenResearch/bamMetrics)
- [GATK (genomeanalysis toolkit) and GATK QUEUE >= 3.2-2](https://www.broadinstitute.org/gatk/)
- [Picard >= 1.119](http://broadinstitute.github.io/picard/) 
- [SnpEff / SnpSift](http://snpeff.sourceforge.net/)
- [Samtools](http://www.htslib.org/)
- [Vcftools](http://vcftools.sourceforge.net/)
- [IGVtools](https://www.broadinstitute.org/igv/igvtools)
- [Contra](http://contra-cnv.sourceforge.net/)
- [FREEC](http://bioinfo-out.curie.fr/projects/freec/)
- [Varscan](http://varscan.sourceforge.net/)
- [Strelka](https://sites.google.com/site/strelkasomaticvariantcaller/)
- [Freebayes](https://github.com/ekg/freebayes)

#### Perl modules
- strict
- POSIX
- Getopt::Long
- FindBin
- File::Path
- Number::Format

#### R packages
- ggplot2
- knitr
- markdown
- reshape
- xtable
- tools
- brew
