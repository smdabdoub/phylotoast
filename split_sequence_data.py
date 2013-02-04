#!/usr/bin/env python

'''
Created on Jan 31, 2013

@author: Shareef Dabdoub


Split an input FASTA-formatted file into a user-specified number of 
smaller files such that the data is evenly distributed among them.
'''
import util
import argparse
from itertools import cycle
import os
import os.path as osp


def split_data(fastaFN, partitions):
    buckets = [[] for _ in xrange(partitions)]
    ring = cycle(xrange(partitions))
    
    with open(fastaFN, 'rU') as inF:
        for seq in util.parseFASTA(inF):
            buckets[ring.next()].append(seq)
    
    return buckets


def handle_program_options():
    parser = argparse.ArgumentParser(description="Split an input \
                                     FASTA-formatted sequence file into a \
                                     user-specified number of smaller files \
                                     such that the sequence data is evenly \
                                     distributed among them.")
    parser.add_argument('-i','--input_fasta_fn', required=True, 
                        help="The sequence data file to be split up.")    
    parser.add_argument('-n', '--num_output_files', type=int, default=2,
                        help="The number of files the input data should be \
                              split into.")
    parser.add_argument('-o', '--output_dir', default='.',
                        help="The location to write the split data files.")
    parser.add_argument('-v', '--verbose', action='store_true')
    
    return parser.parse_args()
        

    
def main():
    args = handle_program_options()    
    buckets = split_data(args.input_fasta_fn, args.num_output_files)
    util.ensure_dir(args.output_dir)
    
    # write out split files
    for i, bucket in enumerate(buckets):
        with open(osp.join(args.output_dir,'%i.fna'%i), 'w') as outF:
            outF.write(''.join(['>{0.id} {0.descr}\n{0.data}\n'.format(r) 
                                                             for r in bucket]))
    if args.verbose:
        msg = '{0} files generated with ~{1} sequences per file.'
        print msg.format(len(buckets), len(buckets[0]))


if __name__ == '__main__':
    main()