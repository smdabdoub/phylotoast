#!/usr/bin/env python

import argparse
import util

def sOTU_name(tax):
    """
    Determine a simple Genus-species identifier for an OTU, if possible. 
    If OTU is not identified to the species level, name it as Unclassified (familly/genus/etc...)
    """
    tax = tax.split('; ')
    for i, lvl in enumerate(tax):
        lvl = lvl.strip()
        if i < len(tax) - 1 and len(tax[i + 1].strip()) == 3:
            if tax[i].strip()[0] == 'g':
                return lvl.split('_')[-1] + '_spp.'
            else:
                return 'Unclassified_' + lvl.split('_')[-1]
        elif i == len(tax) - 1:
            name = lvl.split('_')[-1]
            if lvl[0] == 's':
                name = tax[i-1].split('_')[-1] + '_' + name
            return name
            

def handle_program_options():
    parser = argparse.ArgumentParser(description="Convert a list of OTU IDs to a list of\
                                                  OTU IDs paired with Genus_species \
                                                  identifiers.")
    parser.add_argument('-i', '--otu_id_fp', required=True,
                        help="A text file containing a list (one per line) of OTU IDs.")
    parser.add_argument('-t', '--taxonomy_fp', required=True,
                        help="A file associating OTU ID with a full taxonomic specifier.")
                        
    parser.add_argument('-o', '--output_fp', default='otu_taxonomy.txt',
                        help="A new file containing a list OTU IDs paired with \
                              of short taxonomic identifiers.")

#    parser.add_argument('-v', '--verbose', action='store_true')

    return parser.parse_args()


def main():
    args = handle_program_options()
    
    with open(args.otu_id_fp, 'rU') as otuF:
        otu_ids = [line.strip() for line in otuF.readlines()]
    taxa = util.parse_taxonomy_table(args.taxonomy_fp)
    
    with open(args.output_fp, 'w') as outF:
        for ID in otu_ids:
            if ID in taxa:
                outF.write('{}\t{}\n'.format(ID, sOTU_name(taxa[ID])))
            else:
                print 'OTU ID {} not found in supplied taxonomy file'.format(ID)
    
if __name__ == '__main__':
    main()