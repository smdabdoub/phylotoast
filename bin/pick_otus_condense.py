#!/usr/bin/env python
'''
Created on Nov 27, 2012

Author: Shareef M Dabdoub
Author Matthew Mason

Step 3 of the condensing process.
'''
import argparse
from collections import defaultdict


def condense_otus(otuF, nuniqueF):
    """
    Traverse the input otu-sequence file, collect the non-unique OTU IDs and
    file the sequences associated with then under the unique OTU ID as defined
    by the input matrix.

    :@type otuF: file
    :@param otuF: The output file from QIIME's pick_otus.py
    :@type nuniqueF: file
    :@param nuniqueF: The matrix of unique OTU IDs associated to the list of
                      non-unique OTU IDs they replaced.

    :@rtype: dict
    :@return: The new condensed table of unique OTU IDs and the sequence IDs
              associated with them.
    """
    uniqueOTUs = set()
    nuOTUs = {}

    # parse non-unique otu matrix
    for line in nuniqueF:
        line = line.split()
        uOTU = line[0]
        for nuOTU in line[1:]:
            nuOTUs[nuOTU] = uOTU
        uniqueOTUs.add(uOTU)

    otuFilter = defaultdict(list)
    # parse otu sequence file
    for line in otuF:
        line = line.split()
        otuID, seqIDs = line[0], line[1:]
        if otuID in uniqueOTUs:
            otuFilter[otuID].extend(seqIDs)
        elif otuID in nuOTUs:
            otuFilter[nuOTUs[otuID]].extend(seqIDs)

    return otuFilter


def handle_program_options():
    parser = argparse.ArgumentParser(description="Step 3 of the condensing \
                                     process. Condense the QIIME pick_otus.py \
                                     script output by moving the sequences \
                                     associated with non-unique OTUs to OTU \
                                     IDs that were identified as unique.")
    parser.add_argument('-s', '--seqs_otus', required=True,
                        help="The list of OTU IDs and their associated \
                              sequence IDs.")

    parser.add_argument('-n', '--non_unique_otu_matrix', required=True,
                        help="The list of unique OTU IDs associated with the \
                              OTU IDs they replaced.")

    parser.add_argument('-o', '--condensed_seqs_otus_file',
                        default='condensed_seqs_otus.txt',
                        help="The condensed set of OTU IDs and the matching \
                              sequences. By default outputs to \
                              condensed_seqs_otus.txt")
    parser.add_argument('-v', '--verbose', action='store_true')

    return parser.parse_args()


def main():
    args = handle_program_options()

    with open(args.seqs_otus, 'rU') as sotuF, open(args.non_unique_otu_matrix, 'rU') as nuotuF:
        filteredOTUs = condense_otus(sotuF, nuotuF)

    with open(args.condensed_seqs_otus_file, 'w') as outF:
        for otuID, seqIDs in filteredOTUs.iteritems():
            outF.write("{0}\t{1}\n".format(otuID, '\t'.join(seqIDs)))

    if args.verbose:
        print 'Output written to {}'.format(args.condensed_seqs_otus_file)


if __name__ == '__main__':
    main()
