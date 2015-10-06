#!/usr/bin/env python
'''
This is a quick and dirty script that fits a need I have,
but should be taken more as an example to build on for
your specific needs.

Author: Shareef M Dabdoub
'''
import argparse
import sys
import json
from collections import defaultdict


def levelData(counts, seqCount, level='g'):
    lvlData = sorted([(k, counts[k]) for k in counts if level+'_' in k], key=lambda x: x[1], reverse=True)
    return ['{}: {:4.2f}%'.format(p[0], p[1]/seqCount*100) for p in lvlData]


def lookup(s, d):
    for idx, item in enumerate(d):
        if s in item:
            return d[idx]
    return s+': None'


def summarize_taxa(biom):
    """
    Given an abundance table, group the counts by every
    taxonomic level.
    """
    tamtcounts = defaultdict(int)
    tot_seqs = 0.0

    for row, col, amt in biom['data']:
        tot_seqs += amt
        rtax = biom['rows'][row]['metadata']['taxonomy']
        for i, t in enumerate(rtax):
            t = t.strip()
            if i == len(rtax)-1 and len(t) > 3 and len(rtax[-1]) > 3:
                t = 's__'+rtax[i-1].strip().split('_')[-1]+'_'+t.split('_')[-1]
            tamtcounts[t] += amt

    lvlData = {lvl: levelData(tamtcounts, tot_seqs, lvl) for lvl in ['k', 'p', 'c', 'o', 'f', 'g', 's']}

    return tot_seqs, lvlData


def general_taxon_summary(lvlData, tot_seqs):
    print 'Total sequence count: %i\n' % tot_seqs

    print 'Top families:\n'
    print lvlData['f'][:20]

    print '\nTop genera:\n'
    print lvlData['g'][:20]

    print '\nTop species:\n'
    print lvlData['s'][:20]


def oral_taxon_summary(lvlData):
    print '\nRed Complex:'
    print lookup('s__Porphyromonas_gingivalis', lvlData['s'])
    print lookup('s__Treponema_denticola', lvlData['s'])
    print lookup('s__Tannerella_forsythia', lvlData['s'])

    print '\nOrange Complex: '
    print lookup('s__Fusobacterium_nucleatum', lvlData['s'])
    print lookup('s__Prevotella_intermedia', lvlData['s'])
    print lookup('s__Prevotella_nigrescens', lvlData['s'])
    print lookup('s__Peptostreptococcus_micros', lvlData['s'])

    print '\nYellow Complex: '
    print lookup('s__Streptococcus_sanguis', lvlData['s'])
    print lookup('s__Streptococcus_oralis', lvlData['s'])
    print lookup('s__Streptococcus_mitis', lvlData['s'])
    print lookup('s__Streptococcus_gordonii', lvlData['s'])
    print lookup('s__Streptococcus_intermedius', lvlData['s'])

    print '\nGreen Complex: '
    print lookup('s__Aggregatibacter_actinomycetemcomitans', lvlData['s'])
    print lookup('s__Campylobacter_concisus', lvlData['s'])
    print lookup('s__Eikenella_corrodens', lvlData['s'])
    print lookup('g__Capnocytophaga', lvlData['g'])

    print '\nPurple Complex: '
    print lookup('s__Veillonella_parvula', lvlData['s'])
    print lookup('s__Actinomyces_odontolyticus', lvlData['s'])
    print lookup('s__Selenomonas_noxia', lvlData['s'])
    print lookup('s__Actinomyces_naeslundii', lvlData['s'])


def handle_program_options():
    parser = argparse.ArgumentParser(description="Print a taxonomic summary of \
                                     a given BIOM abundance table.")
    parser.add_argument('-i', '--otu_table', required=True,
                        help="The biom-format file with OTU-Sample abundance \
                              data.")
    parser.add_argument('--oral_taxa', action='store_true',
                        help="If specified, print an additional summary of taxa \
                              important to the oral microbiome (color complexes).")

    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.otu_table):
            pass
    except IOError as ioe:
        sys.exit(
            '\nError with OTU-Sample abundance data file:{}\n'
            .format(ioe)
        )

    with open(args.otu_table, 'rU') as bF:
        biom = json.load(bF)
        tot_seqs, lvlData = summarize_taxa(biom)

        # print summaries
        general_taxon_summary(lvlData, tot_seqs)

        if args.oral_taxa:
            oral_taxon_summary(lvlData)


if __name__ == '__main__':
    main()
