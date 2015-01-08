#!/usr/bin/env python
# coding: utf-8
'''
Take a subset BIOM table (e.g. from a core calculation) and a
representative set (repset) FASTA file and create a new repset
restricted to the OTUs in the BIOM table.

@author: Shareef M Dabdoub
'''
import argparse
import json
from phylotoast import util


def handle_program_options():
    parser = argparse.ArgumentParser(description="Take a subset BIOM table\
                                     (e.g. from a core calculation) and a\
                                     representative set (repset) FASTA file and\
                                     create a new repset restricted to the OTUs\
                                     in the BIOM table.")
    parser.add_argument('-i', '--biom_fp', required=True,
                        help="Path to a biom-format file with OTU-Sample\
                              abundance data.")
    parser.add_argument('-r', '--repset_fp', required=True, 
                        help='Path to a FASTA-format file containing the\
                              representative set of OTUs')
    parser.add_argument('-o', '--repset_out_fp', default='filtered_repset.fna',
                        help="Path to the new restricted repset file")
    
    return parser.parse_args()


def main():
    args = handle_program_options()

    with open(args.biom_fp, 'rU') as bf:
        biom_otus = {row['id'] for row in  json.load(bf)['rows']}

    repset = util.parseFASTA(args.repset_fp)

    with open(args.repset_out_fp, 'w') as out_f:
        fasta_str = ">{} {}\n{}\n"
        for seq in repset:
            if seq.id in biom_otus:
                out_f.write(fasta_str.format(seq.id, seq.descr, seq.data))


if __name__ == '__main__':
    main()
