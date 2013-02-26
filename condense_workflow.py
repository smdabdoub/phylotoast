#!/usr/bin/env python
'''
Created on Feb 25, 2013

@author: Shareef Dabdoub
'''
import argparse
import subprocess

def handle_program_options():
    parser = argparse.ArgumentParser(description="This workflow script will run\
                                     all three steps of the OTU condensing \
                                     pipeline automatically with the default \
                                     output file settings.")
    parser.add_argument('-i', '--assigned_taxonomy_fn', required=True,
                        help="The taxonomy file output by the assign_taxonomy \
                              script.")
    parser.add_argument('-r', '--rep_set_fn', required=True,
                        help="The set of representative sequences.")
    parser.add_argument('-s', '--seqs_otus_fn', required=True,
                        help="The list of OTU IDs and their associated \
                              sequence IDs.")
    parser.add_argument('-L', '--phylogenetic_level', default='s',
                        choices = ['k','p','c','o','f','g','s'],
                        help="Set the phylogenetic level at which to define \
                              OTUs for condensing and downstream processing.\
                              Defaults to species level.")
    
    

    parser.add_argument('-v', '--verbose', action='store_true')
    
    return parser.parse_args()


def main():
    args = handle_program_options()

    output = []
    
    if args.verbose:
        print "Condensing OTUs:\n"
        
        print "Step 1: Condensing assigned taxonomy file...\n"
    output.append(subprocess.check_output(['otu_condense.py', 
                                           '-i', args.assigned_taxonomy_fn,
                                           '-l', args.phylogenetic_level,
                                           '-v']))
    if args.verbose:
        print "output:\n"
        print output[0]
    
        print "Step 2: Condensing representative set...\n"
    output.append(subprocess.check_output(['filter_rep_set.py', 
                                           '-r',args.rep_set_fn,
                                           '-u',
                                           'condensed_assigned_taxonomy.txt',
                                           '-v']))
    if args.verbose:
        print "output:\n"
        print output[1]
        
        print "Step 3: Condensing pick otus output file...\n"
    output.append(subprocess.check_output(['pick_otus_condense.py', 
                                           '-s', args.seqs_otus_fn,
                                           '-n', 'nonunique_otu_matrix.txt',
                                           '-v']))
    
    if args.verbose:
        print "output:\n"
        print output[2]



if __name__ == '__main__':
    main()