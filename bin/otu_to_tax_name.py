#!/usr/bin/env python

import argparse
import sys
from phylotoast import otu_calc as otuc, util


def handle_program_options():
    parser = argparse.ArgumentParser(description="Convert a list of OTU IDs to a list of\
                                                  OTU IDs paired with Genus_species \
                                                  identifiers.")
    parser.add_argument('-i', '--otu_id_fp', required=True,
                        help="Either a text file containing a list (one per line) \
                        of OTU IDs, or a tab-separated (classic) BIOM-format file.")
    parser.add_argument('-t', '--taxonomy_fp', required=True,
                        help="A file associating OTU ID with a full taxonomic specifier.")

    parser.add_argument('-o', '--output_fp', default='otu_taxonomy.txt',
                        help="For a list input, a new file containing a list of OTU IDs \
                              and their corresponding short taxonomic identifiers \
                              separated by tabs. For a BIOM file input, a new \
                              mapping file with all the OTU IDs replaced by the short\
                              identifier.")

    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.otu_id_fp):
            pass
    except IOError as ioe:
        sys.exit('\nError with file containing OTUIDs/BIOM format:{}\n'.format(ioe))

    with open(args.otu_id_fp, 'rU') as otuF:
        otu_ids = [line.strip().split('\t') for line in otuF.readlines()]
    taxa = util.parse_taxonomy_table(args.taxonomy_fp)

    with open(args.output_fp, 'w') as outF:
        for entry in otu_ids:
            if isinstance(entry, list):
                # check for comments in BIOM files
                if not entry[0][0] == '#':
                    ID = entry[0]
                else:
                    outF.write('{}\n'.format('\t'.join(entry)))
                    continue
            # instead of a BIOM file, a line-by-line list of OTU IDs
            else:
                ID = entry

            if ID in taxa:
                named_ID = otuc.otu_name(taxa[ID].split('; '))
            else:
                print 'Error: OTU ID {} not found in supplied taxonomy file; stopping...'.format(ID)
                return

            # write out to file
            out_str = '{}\t{}\n'
            if isinstance(entry, list):
                outF.write(out_str.format(named_ID, '\t'.join(entry[1:])))
            else:
                outF.write(out_str.format(ID, named_ID))

if __name__ == '__main__':
    main()
