#!/usr/bin/env python

from __future__ import division
import argparse
from collections import defaultdict

def multi_symmdiff(*sets):
    d = set.symmetric_difference(*sets[:2])
    for s in sets[2:]:
        d = set.symmetric_difference(d, s)
    return d
    
def combine_primer_otu(otuid, primers):
    seqs = []
    for p in primers:
        if otuid in primers[p]:
            seqs.extend(primers[p][otuid])
    return seqs

def handle_program_options():
    parser = argparse.ArgumentParser(description="Combine multi-primer pick \
                                     OTUs results files into a single results \
                                     file while at the same time averaging \
                                     sequence counts per sample for OTUs \
                                     shared between the primer-set results.\
                                     See reference: Kumar PS et al. (2011) \
                                     doi:10.1371/journal.pone.0020956")
    parser.add_argument('--p1', required=True,
                        help="Primer-set 1 seqs_otus results files.")
    parser.add_argument('--p2', required=True,
                        help="Primer-set 2 seqs_otus results files.")
    
    parser.add_argument('-o', '--output_fp', default='combined_seqs_otus.txt',
                        help="The combined seqs_otus file that has been \
                              averaged by shared OTU entries.\
                              Default: combined_seqs_otus.txt")
    
    parser.add_argument('-v', '--verbose', action='store_true',
                        help="Print detailed information about script \
                              operation.")
    
    return parser.parse_args()


def main():
    args = handle_program_options()
    
    # gather input seq_otus data
    primers = {}
    for i,soFP in enumerate([args.p1, args.p2]):
        with open(soFP, 'rU') as inF:
            otus = {}
            for line in inF:
                line = line.strip().split('\t')
                otuid, seqs = line[0], line[1:]
                otus[otuid] = seqs
            primers[i] = otus.copy()
    
    primer_otus = [set(primers[p].keys()) for p in primers]
    shared_otus = set.intersection(*primer_otus)
    nonshared_otus = multi_symmdiff(*primer_otus)
    
    outF = open(args.output_fp, 'w')
    
    for otuid in nonshared_otus:
        seqs = combine_primer_otu(otuid, primers)
        outF.write('{}\t{}\n'.format(otuid, '\t'.join(seqs)))
    
    for otuid in shared_otus:
        samples = {}
        # collect sequences by sample id into bins per primer
        for p in primers:
            for s in primers[p][otuid]:
                sid = s.split('_')[0]
                if sid not in samples:
                    samples[sid] = {}
                if p not in samples[sid]:
                    samples[sid][p] = []
                samples[sid][p].append(s)

        if args.verbose:
            print 'Shared OTU:', otuid
        
        # write sequence information for OTU
        outF.write(otuid+'\t')
        for i,sid in enumerate(samples):
            sout = []
            # sample is in both primers, only keep 1/2 of sequences to adjust for bias
            if len(samples[sid]) == 2:
                for p in samples[sid]:
                    sout.extend(samples[sid][p])
                lsb = len(sout)
                sout = sout[:int(round(len(sout)/2))]
                lsa = len(sout)
                print '\tSample {} reduced from {} to {} sequences'.format(sid, lsb, lsa)
            else:
                sout = samples[sid][samples[sid].keys()[0]]
                print '\tSample {}: {} sequences'.format(sid, len(sout))
            
            outF.write('\t'.join(sout))
            if i < len(samples) - 1:
                outF.write('\t')
        outF.write('\n')

    outF.close()
    

if __name__ == '__main__':
    main()