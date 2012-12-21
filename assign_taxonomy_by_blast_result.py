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


def handle_program_options():
    parser = argparse.ArgumentParser(description="Assign taxonomy to a rep \
                                     set of OTUs that were chosen by BLAST \
                                     from an annotated database.")
    parser.add_argument('rep_set_fp', help="The set of representative sequences.")
    parser.add_argument('id_to_taxonomy_fp', help="Path to tab-delimited file \
                                      mapping sequences to assigned taxonomy.")
    
    parser.add_argument('-o', '--assigned_taxonomy_fp', 
                        default='assigned_taxonomy.tx',
                        help="The path to the result file. By default \
                              outputs to assigned_taxonomy.txt")
    parser.add_argument('-v', '--verbose', action='store_true')
    
    return parser.parse_args()


if __name__ == '__main__':
    pass