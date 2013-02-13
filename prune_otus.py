#!/usr/bin/env python

'''
Created on Feb 11, 2013

@author: Shareef Dabdoub

Parse the OTU-sequence data in two steps. First remove any OTUs that occur in 
less than a user-defined percent of samples (default 1%). Second, remove any 
OTUs that make up less than a user-defined percentage of the overall 
sequences (default 0.01%)
'''
import argparse

def filter_by_sample_pct(otus, sample_count, pct, phyl_level):
    """
    Split the list of OTUs (and associated sequence ids) into two lists:
    those occurring in more than some percentage of samples and those less than
    the cutoff.
    
    :@type otus: dict
    :@param otus: {otuid: [taxonomy, [sequence IDs]]}
    :@type sample_count: int
    :@param sample_count: The total number of samples in the data set
    :@type pct: float
    :@param pct: The cutoff percentage for inclusion in the filtered 
                 set of OTUs
    :@type phyl_level: str
    :@param phyl_level: The phylogenetic level (e.g. family, group, etc...) at 
                       which to combine OTU counts for thresholding. One of 
                       the following: ['k','p','c','o','f','g','s']
                       
    :@rtype: tuple
    :@return: Two dicts: the OTU IDs and sequence IDs above and below the 
              percentage threshold.
    """
    if phyl_level not in ['k','p','c','o','f','g','s']:
        phyl_level = 's'
    sample_count = float(sample_count)
    sample_counts = {split_phylogeny(item[0], phyl_level):0
                       for item in otus.values()}
    # count the number of sequences per OTU
    for otuid in otus:
        phyl = split_phylogeny(otus[otuid][0], phyl_level)
        samples = {seqid.split('_')[0] for seqid in otus[otuid][1]}
        sample_counts[phyl] += len(samples)
    sample_counts = {phyl:sample_counts[phyl]/sample_count 
                       for phyl in sample_counts}

    above = {oid:otus[oid][1] for oid in otus 
             if sample_counts[split_phylogeny(otus[oid][0],phyl_level)] >= pct}
    below = {}
    for otuid in otus:
        phyl = split_phylogeny(otus[otuid][0], phyl_level)
        sample_pct = sample_counts[phyl]/sample_count
        if sample_pct < pct:
            below[otuid] = [sample_pct, '', otus[otuid]]
    
    return above, below


def filter_by_sequence_pct(otus, pct):
    seq_counts = {oid:sum(otus[oid]) for oid in otus}
    nseqs = float(sum(seq_counts.values()))
    above = {oid:otus[oid] for oid in otus if seq_counts[oid]/nseqs >= pct}
    below = {}
    for oid in otus:
        seq_pct = seq_counts[oid]/nseqs
        if seq_counts[oid]/nseqs < pct:
            below[oid] = ['', seq_pct, otus[oid]] 
    
    return above, below


def gather_otus_samples(inFN):
    otus = {}
    with open(inFN, 'rU') as seqsF:
        for line in seqsF:
            line = line.strip().split('\t')
            otus[line[0]] = line[1:]
    samples = set()
    for seqs in otus.values():
        samples.update({seq.split('_')[0] for seq in seqs})


def assign_taxonomy(otus, taxFN):
    fsOTUs = frozenset(otus)
    with open(taxFN) as taxF:
        return {otu:tax for otu,tax in (line.split('\t') for line in taxF) 
                  if otu in fsOTUs}

        
def split_phylogeny(p, level='s'):
    level = level+'__'
    result = p.split(level)
    return result[0]+level+result[1].split(';')[0]


def handle_program_options():
    parser = argparse.ArgumentParser(description="Parse the OTU-sequence data \
                                     in two steps. First remove any OTUs that \
                                     occur in less than a user-defined percent\
                                     of samples (default 1%). Second, remove \
                                     any OTUs that make up less than a \
                                     user-defined percentage of the overall \
                                     sequences (default 0.01%)")
                                     
    parser.add_argument('-i', '--seqs_otus_fp', required=True,
                        help="The output from the pick OTUs step, e.g. \
                              seqs_otus.txt")
    parser.add_argument('-t','--id_to_taxonomy_fp', required=True, 
                        help="Path to tab-delimited file mapping sequences to \
                        assigned taxonomy.")
    parser.add_argument('-p', '--percent_of_samples', type=float,
                        default=0.01, help="OTUs that occur in less than this \
                                         percent of samples will be removed. \
                                         Default is 1%.")
    parser.add_argument('-s', '--percent_of_sequences', type=float,
                        default=0.0001, help="OTUs that occur in less than \
                                              this percent of total sequences \
                                              will be removed. \
                                              Default is 0.01%")
    parser.add_argument('-l', '--phylogenetic_level', default='g',
                        choices = ['k','p','c','o','f','g','s'],
                        help="Set the phylogenetic level at which to join \
                              OTUs for consideration in pruning. \
                              Default is Group.")
    parser.add_argument('-o', '--output_pruned_otus_fn',
                        default='seqs_otus_pruned.txt', 
                        help="The main output file that will contain the \
                              remaining OTUs and sequence IDs.")
    parser.add_argument('--output_removed_otus_fn', 
                        default='removed_otus.txt',
                        help="The file to write out the set of OTUs that were \
                              removed by the filter.")
    #parser.add_argument('-v', '--verbose', action='store_true')
    
    return parser.parse_args()


def main():
    args = handle_program_options()
    
    seqs_otus, sampleIDs = gather_otus_samples(args.seqs_otus_fp)
    otu_taxa = assign_taxonomy(seqs_otus.keys(), args.id_to_taxonomy_fp)
    
    otus = {}
    for otuid, seqids in seqs_otus.iteritems():
        otus[otuid] = (otu_taxa[otuid], seqids)
        

    above, below = filter_by_sample_pct(otus, len(sampleIDs), 
                                        args.percent_of_samples,
                                        args.phylogenetic_level)
    above, below2 = filter_by_sequence_pct(above, args.percent_of_sequences)
    below.update(below2)
    
    with open(args.output_pruned_otus_fn, 'w') as outF:
        for otuid, seqids in above.iteritems():
            outF.write('{0}\t{1}\n'.format(otuid, '\t'.join(seqids)))
    
    with open(args.output_removed_otus_fn, 'w') as outF:
        for oid in otus:
            line = '{0}\t{1}\t{2}\t{3}'
            outF.write(line.format(oid, otus[oid][0], otus[oid][1], 
                                   '\t'.join(otus[oid][2])))


if __name__ == '__main__':
    main()