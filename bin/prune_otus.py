#!/usr/bin/env python

'''
Created on Feb 11, 2013

Author: Shareef Dabdoub

Parse the OTU-sequence data in two steps. First remove any OTUs that occur in
less than a user-defined percent of samples (default 1%). Second, remove any
OTUs that make up less than a user-defined percentage of the overall
sequences (default 0.01%)
'''
import sys
import argparse
from collections import defaultdict
from phylotoast import util

def filter_by_sample_pct(otus, nsamples, pct, phyl_level):
    """
    Split the list of OTUs (and associated sequence ids) into two lists:
    those occurring in more than some percentage of samples and those less than
    the cutoff.

    :type otus: dict
    :param otus: {otuid: [taxonomy, [sequence IDs]]}
    :type nsamples: int
    :param nsamples: The total number of samples in the data set
    :type pct: float
    :param pct: The cutoff percentage for inclusion in the filtered
                 set of OTUs
    :type phyl_level: str
    :param phyl_level: The phylogenetic level (e.g. family, group, etc...) at
                       which to combine OTU counts for thresholding. One of
                       the following: ['k','p','c','o','f','g','s']

    :rtype: tuple
    :return: Two dicts: the OTU IDs and sequence IDs above and below the
              percentage threshold.
    """
    if phyl_level not in ['k', 'p', 'c', 'o', 'f', 'g', 's']:
        phyl_level = 's'
    nsamples = float(nsamples)
    sample_counts = defaultdict(set)
    # count the number of sequences per OTU
    for otuid in otus:
        phyl = util.split_phylogeny(otus[otuid][0], phyl_level)
        samples = {seqid.split('_')[0] for seqid in otus[otuid][1]}
        sample_counts[phyl].update(samples)
    sample_counts = {phyl: len(sample_counts[phyl])/nsamples
                     for phyl in sample_counts}

    # separate OTUs
    above = {}
    below = {}
    for otuid in otus:
        phyl = util.split_phylogeny(otus[otuid][0], phyl_level)
        if sample_counts[phyl] >= pct:
            above[otuid] = otus[otuid]
        else:
            below[otuid] = [sample_counts[phyl], '',
                            otus[otuid][0], otus[otuid][1]]

    return above, below


def filter_by_sequence_pct(otus, nseqs, pct, phyl_level):
    """
    Split the list of OTUs (and associated sequence ids) into two lists:
    those occurring associated with more than some percentage of total sequences
    and those less than the cutoff.

    :type otus: dict
    :param otus: {otuid: [taxonomy, [sequence IDs]]}
    :type nseqs: int
    :param nseqs: The total number of sequences in the data set
    :type pct: float
    :param pct: The cutoff percentage for inclusion in the filtered
                 set of OTUs
    :type phyl_level: str
    :param phyl_level: The phylogenetic level (e.g. family, group, etc...) at
                       which to combine OTU counts for thresholding. One of
                       the following: ['k','p','c','o','f','g','s']

    :rtype: tuple
    :return: Two dicts: the OTU IDs and sequence IDs above and below the
              percentage threshold.
    """
    if phyl_level not in ['k', 'p', 'c', 'o', 'f', 'g', 's']:
        phyl_level = 's'
    seq_counts = defaultdict(int)
    nseqs = float(nseqs)
    # gather counts
    for oid in otus:
        phyl = util.split_phylogeny(otus[oid][0], phyl_level)
        seq_counts[phyl] += len(otus[oid][1])
    seq_counts = {phyl: seq_counts[phyl]/nseqs for phyl in seq_counts}

    # separate OTUs
    above = {}
    below = {}
    for otuid in otus:
        phyl = util.split_phylogeny(otus[otuid][0], phyl_level)
        if seq_counts[phyl] >= pct:
            above[otuid] = otus[otuid]
        else:
            below[otuid] = ['', seq_counts[phyl],
                            otus[otuid][0], otus[otuid][1]]

    return above, below


def gather_otus_samples(inFN):
    otus = {}
    with open(inFN, 'rU') as seqsF:
        for line in seqsF:
            line = line.strip().split('\t')
            otus[line[0]] = line[1:]
    samples = set()
    nseqs = 0
    for seqs in otus.values():
        samples.update({seq.split('_')[0] for seq in seqs})
        nseqs += len(seqs)

    return otus, len(samples), nseqs


def assign_taxonomy(otus, taxFN):
    fsOTUs = frozenset(otus)
    with open(taxFN) as taxF:
        return {otu: tax for otu, tax in (line.split('\t') for line in taxF)
                if otu in fsOTUs}


def handle_program_options():
    parser = argparse.ArgumentParser(description="Parse the OTU-sequence data \
                                     in two steps. First remove any OTUs that \
                                     occur in less than a user-defined percent\
                                     of samples (default 5%). Second, remove \
                                     any OTUs that make up less than a \
                                     user-defined percentage of the overall \
                                     sequences (default 0.01%)")

    parser.add_argument('-i', '--seqs_otus_fn', required=True,
                        help="The output from the pick OTUs step, e.g. \
                              seqs_otus.txt")
    parser.add_argument('-t', '--id_to_taxonomy_fn', required=True,
                        help="Path to tab-delimited file mapping sequences to \
                        assigned taxonomy.")
    parser.add_argument('-p', '--percent_of_samples', type=float,
                        default=0.05, help="OTUs that occur in less than this \
                                            percent of samples will be removed.\
                                            Default is 5 percent.")
    parser.add_argument('-s', '--percent_of_sequences', type=float,
                        default=0.0001, help="OTUs that occur in less than \
                                              this percent of total sequences \
                                              will be removed. \
                                              Default is 0.01 percent.")
    parser.add_argument('-l', '--phylogenetic_level', default='g',
                        choices=['k', 'p', 'c', 'o', 'f', 'g', 's'],
                        help="Set the phylogenetic level at which to join \
                              OTUs for consideration in pruning. \
                              Default is 'g'(group).")
    parser.add_argument('-o', '--output_pruned_otus_fn',
                        default='seqs_otus_pruned.txt',
                        help="The main output file that will contain the \
                              remaining OTUs and sequence IDs.")
    parser.add_argument('--output_removed_otus_fn',
                        default='removed_otus.txt',
                        help="The file to write out the set of OTUs that were \
                              removed by the filter.")
    parser.add_argument('-v', '--verbose', action='store_true')

    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.seqs_otus_fn):
            pass
    except IOError as ioe:
        sys.exit('\nError with output file from pick OTUs step:{}\n'.format(ioe))

    try:
        with open(args.id_to_taxonomy_fn):
            pass
    except IOError as ioe:
        sys.exit('\nError with file mapping seqences to asssigned taxonomy:{}\n'.format(ioe))

    seqs_otus, nsamples, nseqs = gather_otus_samples(args.seqs_otus_fn)
    otu_taxa = assign_taxonomy(seqs_otus.keys(), args.id_to_taxonomy_fn)

    otus = {}
    for otuid, seqids in seqs_otus.iteritems():
        otus[otuid] = (otu_taxa[otuid], seqids)

    above, below = filter_by_sample_pct(otus, nsamples,
                                        args.percent_of_samples,
                                        args.phylogenetic_level)

    above, below2 = filter_by_sequence_pct(above, nseqs,
                                           args.percent_of_sequences,
                                           args.phylogenetic_level)

    above2, below3 = filter_by_sequence_pct({boid: below[boid][2:]
                                            for boid in below},
                                            nseqs,
                                            args.percent_of_sequences,
                                            args.phylogenetic_level)
    below = {boid: below[boid] for boid in below3}
    below.update(below2)
    above.update(above2)

    with open(args.output_pruned_otus_fn, 'w') as outF:
        for otuid, item in above.iteritems():
            outF.write('{0}\t{1}\n'.format(otuid, '\t'.join(item[1])))

    with open(args.output_removed_otus_fn, 'w') as outF:
        outF.write('OTU ID\tSample%\tSeq %\tSequence IDs\n')
        for oid, item in below.iteritems():
            seqpct = '{seqpct:.4f}' if item[0] != '' else '     '
            samplepct = '{samplepct:.2G}' if item[1] != '' else '     '
            line = '{otuid}\t'+seqpct+'\t'+samplepct+'\t{seqs}\n'
            outF.write(line.format(otuid=oid, seqpct=item[0],
                                   samplepct=item[1], seqs='\t'.join(item[3])))

    if args.verbose:
        print 'Input: \t{} total samples'.format(nsamples)
        print '\t{} total sequences\n'.format(nseqs)
        print 'From a total of {} input otus'.format(len(otus))
        print '{} otus remain '.format(len(above))
        print '{} otus removed'.format(len(below))

        phyl_map = {'k': 'kingdoms', 'p': 'phyla', 'c': 'classes', 'o': 'orders',
                    'f': 'families', 'g': 'genera', 's': 'species'}
        phyls = {util.split_phylogeny(otus[oid][0], args.phylogenetic_level)
                 for oid in otus}
        print '\nFrom the {} total {}'.format(len(phyls),
                                              phyl_map[args.phylogenetic_level])
        phyl_above = {util.split_phylogeny(otus[aoid][0], args.phylogenetic_level)
                      for aoid in above}
        phyl_below = {util.split_phylogeny(otus[boid][0], args.phylogenetic_level)
                      for boid in below}
        above_abundance = sum([len(item[1]) for item in above.values()])
        below_abundance = sum([len(below[boid][3]) for boid in below])
        report = ('{0} {1} ({2:.4G}%) were {3}, and account for {4:.4G}% of' +
                  ' all sequence data ({5} sequences)')
        print report.format(len(phyl_above), phyl_map[args.phylogenetic_level],
                            len(phyl_above)/float(len(phyls))*100, 'kept',
                            above_abundance/float(nseqs)*100, above_abundance)
        print report.format(len(phyl_below), phyl_map[args.phylogenetic_level],
                            len(phyl_below)/float(len(phyls))*100, 'removed',
                            below_abundance/float(nseqs)*100, below_abundance)


if __name__ == '__main__':
    main()
