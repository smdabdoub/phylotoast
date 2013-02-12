#!/usr/bin/env python

'''
Created on Feb 11, 2013

@author: Shareef Dabdoub

Parse the OTU-sequence data in two steps. First remove any OTUs that occur in 
less than a user-defined percent of samples (default 1%). Second, remove any 
OTUs that make up less than a user-defined percentage of the overall 
sequences (default 0.01%)
'''
import argparse




def handle_program_options():
    parser = argparse.ArgumentParser(description="Parse the OTU-sequence data \
                                     in two steps. First remove any OTUs that \
                                     occur in less than a user-defined percent\
                                     of samples (default 1%). Second, remove \
                                     any OTUs that make up less than a \
                                     user-defined percentage of the overall \
                                     sequences (default 0.01%)")
                                     
    parser.add_argument('-i', '--seqs_otus_fp',
                        help="The output from the pick otus step, e.g. \
                              seqs_otus.txt")
    parser.add_argument('-o', '--output_pruned_otus_fn',
                        default='seqs_otus_pruned.txt', 
                        help='The name of the file the pruned results will be \
                              written to.')
    parser.add_argument('-v', '--verbose', action='store_true')
    
    return parser.parse_args()


def main():
    args = handle_program_options()








if __name__ == '__main__':
    main()