#!/usr/bin/env python
'''
Created on Nov 27, 2012

@author: Shareef M Dabdoub
@author Matthew Mason

Step 2 of the condensing process.
'''
import argparse

from Bio import SeqIO

def parse_unique_otus(inF):
    """
    Create a list of the OTU IDs from the input file.
    
    :@type inF: file
    :@param inF: The unique OTU file
    
    :@rtype: list
    :@return: The sequence IDs associated with unique OTUs
    """
    return {line.split('\t')[0] for line in inF}
        

def filter_rep_set(inF, otuSet):
    """
    Parse the rep set file and remove all sequences not associated with unique
    OTUs.
    
    :@type inF: file
    :@param inF: The representative sequence set
    
    :@rtype: list
    :@return: The set of sequences associated with unique OTUs 
    """
    seqs = []
    for record in SeqIO.parse(inF, "fasta"):
        if record.id in otuSet:
            seqs.append(record)
    return seqs


def handle_program_options():
    parser = argparse.ArgumentParser(description="Step 2 of the condensing \
                                     process. Filter the representative \
                                     sequence set to include only those \
                                     sequences that map to unique OTUs.")
    parser.add_argument('-r', '--rep_set_fn', required=True,
                        help="The set of representative sequences.")
    parser.add_argument('-u', '--unique_otus_fn', required=True, 
                        help="The condensed assigned taxonomy file.")
    
    parser.add_argument('-o', '--output_filtered_rep_set_fn', 
                        default='condensed_rep_set.fna',
                        help="The filtered representative set. By default \
                              outputs to condensed_rep_set.fna")
    parser.add_argument('-v', '--verbose', action='store_true')
    
    return parser.parse_args()
        

    
def main():
    args = handle_program_options()
    
    with open(args.unique_otus_fn, 'rU') as uoF:
        otuSet = parse_unique_otus(uoF)

    with open(args.rep_set_fn, 'rU') as rsF:
        records = filter_rep_set(rsF, otuSet)

    SeqIO.write(records, args.output_filtered_rep_set_fn, "fasta")
    
    if args.verbose:
        print '%i sequences associated with unique OTUs' % len(records)
        print 
        print 'Output written to: %s' % args.output_filtered_rep_set_fn


if __name__ == '__main__':
    main()