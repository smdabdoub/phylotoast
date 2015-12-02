#!/usr/bin/env python
# coding: utf-8
import argparse
import sys


def handle_program_options():
    parser = argparse.ArgumentParser(description="This filter allows for the \
                                     removal of sequences not contained within \
                                     a user-specified list of Sample IDs. This \
                                     script examines each OTU and removes any \
                                     sequences not originating from the \
                                     specified set of allowed Sample IDs. Any \
                                     empty OTUs that result are removed.")
    parser.add_argument('-i', '--otu_map', required=True,
                        help="path to the input OTU map (i.e., the output from\
                        pick_otus.py) [REQUIRED]")
    parser.add_argument('-k', '--samples_to_keep_fp', required=True,
                        help="path to the file containing Sample IDs to keep \
                              in the new OTU map. One Sample ID per line.")
    parser.add_argument('-o', '--output_otu_map_fp', required=True,
                        help="path to the output filtered OTU map")
    parser.add_argument('-v', '--verbose', action='store_true')

    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.otu_map):
            pass
    except IOError as ioe:
        sys.exit(
            '\nError with input OTU map filepath:{}\n'
            .format(ioe)
        )

    try:
        with open(args.samples_to_keep_fp):
            pass
    except IOError as ioe:
        sys.exit(
            '\nError with file containing SampleID to retain:{}\n'
            .format(ioe)
        )

    seqs_otus = {}
    with open(args.otu_map, 'rU') as otuF:
        for line in otuF:
            line = line.strip().split('\t')
            seqs_otus[line[0]] = line[1:]

    with open(args.samples_to_keep_fp, 'rU') as inF:
        keep = frozenset([line.strip() for line in inF])

    with open(args.output_otu_map_fp, 'w') as outF:
        for otu in seqs_otus:
            seqs = [entry for entry in seqs_otus[otu] if entry[:entry.rfind('_')] in keep]
            if seqs:
                outF.write('{otu}\t{seqs}\n'.format(otu=otu, seqs='\t'.join(seqs)))
            else:
                if args.verbose:
                    print otu, 'removed'

if __name__ == '__main__':
    sys.exit(main())
