import argparse

import pysam
import numpy

def print_sex(bam):
    """
    Print sex based on chr x ratio
    
    Args:
	bam (str): Path to bam file
    """

    idxstats = pysam.idxstats(bam)
    chr_ratio = []
    # Calculate read / chromosome length ratio per chromosome
    for chr in idxstats[0:24]:
	chr = chr.strip('\n').split('\t')
	chr_length = float(chr[1])
	chr_mapped = float(chr[2])
	ratio = chr_mapped / chr_length
	chr_ratio.append(ratio)
	
    chr_ratio_std = numpy.std(chr_ratio)
    chr_ratio_mean = numpy.mean(chr_ratio)

    chr_x = idxstats[22].strip('\n').split('\t')
    chr_x_ratio = float(chr_x[2]) / float(chr_x[1])

    if ( (chr_x_ratio > chr_ratio_mean - (2 * chr_ratio_std)) and (chr_x_ratio < chr_ratio_mean + (2 * chr_ratio_std)) ):
	print 'female'
    elif (chr_x_ratio < chr_ratio_mean - (2 * chr_ratio_std)):
	print 'male'
    else:
	print "unkown"

if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=100, width=200),
	description = 'Determine and print sex based on bam file')
    parser.add_argument('-b','--bam', type=str, action='append', help='Bam file', required=True)
    args = parser.parse_args()
    
    for bam in args.bam:
	print_sex(bam)