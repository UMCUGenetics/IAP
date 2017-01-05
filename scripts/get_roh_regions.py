"""Scripts used to analyse bcftools roh output."""
import sys
import argparse
from itertools import islice
from operator import itemgetter


def list_mean(list):
    """Calculate the mean of a list with numbers."""
    mean = float(sum(list)) / float(len(list))
    return mean


def get_roh_regions(roh_file):
    """Create roh regions."""
    roh_regions = []
    try:
        f = open(roh_file, 'r')
    except IOError:
        sys.exit("Can't open: {0}".format(roh_file))
    else:
        with f:
            # skip header rows == 4, bcftools 1.3
            for line in islice(f, 4):
                next

            # Get first roh_region
            chr, pos, state, qual = next(f).rstrip().split('\t')
            roh_region = {
                'chr': chr,
                'start': int(pos),
                'end': int(pos),
                'state': int(state),
                'qual': [float(qual)]
            }

            for line in f:
                chr, pos, state, qual = line.rstrip().split('\t')
                if (roh_region['chr'] == chr and roh_region['state'] == int(state)):
                    roh_region['end'] = int(pos)
                    roh_region['qual'].append(float(qual))
                else:
                    # Calculate mean qual and length
                    roh_region['qual'] = list_mean(roh_region['qual'])
                    roh_region['len'] = roh_region['end'] - roh_region['start']
                    roh_regions.append(roh_region)
                    roh_region = {
                        'chr': chr,
                        'start': int(pos),
                        'end': int(pos),
                        'state': int(state),
                        'qual': [float(qual)]
                    }

            # Calculate average qual and length
            roh_region['qual'] = list_mean(roh_region['qual'])
            roh_region['len'] = roh_region['end'] - roh_region['start']
            roh_regions.append(roh_region)

    return roh_regions


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=lambda prog: argparse.HelpFormatter(
            prog, max_help_position=100, width=200),
        description='Generate ROH regions from bcftools roh output.')

    parser.add_argument('roh_file', type=str, help='bcftools ROH file')
    parser.add_argument('-s', '--sort_by_size', action='store_true',
                        help='Sort output by region size')
    parser.add_argument('-r', '--roh_only', action='store_true',
                        help='Only output ROH regions')
    args = parser.parse_args()

    roh_regions = get_roh_regions(args.roh_file)

    if args.sort_by_size:
        roh_regions.sort(key=itemgetter('len'), reverse=True)

    # Print output
    print "#Chrom\tStart\tEnd\tLength\tState\tMeanQual"
    for roh_region in roh_regions:
        if not args.roh_only or roh_region['state'] == 1:
            print "{chr}\t{start}\t{end}\t{len}\t{state}\t{qual:.1f}".format(
                chr=roh_region['chr'],
                start=roh_region['start'],
                end=roh_region['end'],
                len=roh_region['len'],
                state=roh_region['state'],
                qual=roh_region['qual'],
            )
