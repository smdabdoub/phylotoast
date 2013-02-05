#!/usr/bin/env python

'''
Created on Dec 21, 2012

Picking OTUs by BLAST with QIIME against an annotated reference set (such as 
GreenGenes results in OTU IDs that are sequence IDs from the reference set 
which can then be looked up directly from the reference set taxonomy table 
(e.g. greengenes_tax.txt), thus bypassing the need for the assign_taxonomy 
methods (BLAST, RDP, RTAX)

@author: Shareef Dabdoub
'''
import argparse

from Bio import SeqIO

def parse_taxonomy_table(idtaxFN):
    idtax = {}
    with open(idtaxFN,'rU') as idtxF:
        for line in idtxF:
            ID, tax = line.strip().split('\t')
            idtax[ID] = tax
    
    return idtax


def handle_program_options():
    parser = argparse.ArgumentParser(description="Assign taxonomy to a rep \
                                     set of OTUs that were chosen by BLAST \
                                     from an annotated database.")
    parser.add_argument('-i','--rep_set_fp', required=True, 
                        help="The set of representative sequences.")
    parser.add_argument('-t','--id_to_taxonomy_fp', required=True, 
                        help="Path to tab-delimited file mapping sequences to \
                        assigned taxonomy.")
    
    parser.add_argument('-o', '--assigned_taxonomy_fp', 
                        default='assigned_taxonomy.txt',
                        help="The path to the result file. By default \
                              outputs to assigned_taxonomy.txt")
    parser.add_argument('-v', '--verbose', action='store_true')
    
    return parser.parse_args()


def main():
    args = handle_program_options()
    
    # input the ID to Taxonomy table and the rep set
    taxids = parse_taxonomy_table(args.id_to_taxonomy_fp)
    rep_set = SeqIO.to_dict(SeqIO.parse(args.rep_set_fp, 'fasta'))
    
    # write out the assigned taxonomy file
    with open(args.assigned_taxonomy_fp, 'w') as outF:
        for taxid in rep_set:
            line = '{0}\t{1}\t{2}\t{0}\n'.format(taxid, taxids[taxid], 0.0)
            outF.write(line)

if __name__ == '__main__':
    main()