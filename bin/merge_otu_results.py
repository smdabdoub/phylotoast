#!/usr/bin/env python

'''
Created on Feb 4, 2013

Author: Shareef Dabdoub

Distributing sequence data across the cluster for OTU picking results
in a set of result files that need to be merged into a single pick otus
result.
'''
import sys
import os
import argparse
from collections import defaultdict


def merge_results(results_FNs):
    otus = defaultdict(list)
    for fn in results_FNs:
        with open(fn, 'rU') as resultF:
            for line in resultF.readlines():
                line = line.split()
                otus[line[0]].extend(line[1:])
    return otus


def handle_program_options():
    parser = argparse.ArgumentParser(description="Distributing sequence data \
                                     across the cluster for OTU picking results\
                                     in a set of result files that need to be \
                                     merged into a single pick otus result.")
    parser.add_argument('pick_otus_results', nargs='+',
                        help="The result files from multiple runs of a pick \
                              otus script that need to be merged.")
    parser.add_argument('-o', '--output_fn', default='seqs_otus.txt',
                        help='The name of the file the merged results will be \
                              written to.')
    parser.add_argument('-v', '--verbose', action='store_true')

    return parser.parse_args()


def main():
    args = handle_program_options()

    for file in args.pick_otus_results:
        if not os.path.exists(file):
            sys.exit('Error! {} could not be found.'.format(file))

    otus = merge_results(args.pick_otus_results)

    with open(args.output_fn, 'w') as outF:
        for otuID in otus:
            line = '{otu_id}\t{sample_ids}\n'
            outF.write(line.format(otu_id=otuID,
                                   sample_ids='\t'.join(otus[otuID])))

    if args.verbose:
        print '{} files merged'.format(len(args.pick_otus_results))
        print '{} otus found'.format(len(otus))


if __name__ == '__main__':
    main()
