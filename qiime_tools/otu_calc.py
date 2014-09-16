from __future__ import division
# python libs
import ast
from collections import (OrderedDict,defaultdict)
import json
# 3rd party
from fuzz import FuzzySet
import numpy as np
# local
import util


def fuzzy_lookup(orig, keys):
    """
    Return the intersection of a fuzzy set and a 
    collection of keys (presumably a subset)
    """
    fset = FuzzySet()
    for k in keys:
        if k in orig:
            fset.add(orig[k])
    return fset


def sdi(fset):
    """
    Calculate the Shannon Diversity Index.
    H = -sum(p*ln(p)) 
    where p is the relative abundance of a single OTU in the set.
    Note that the Equitability Index ( $E_H = H / H_{max}$ ) could be easily 
    calculated from the returned array by:
    >>>diversities = sdi(fset)
    >>>equitabilities = diversities/max(diversities)
    
    :@type fset: FuzzySet
    :@param fset: The set of OTUs and their relative abundance values
    :@rtype: float
    :@return: The Shannon Diversity Index.
    """
    p = np.array([e.mu for e in fset])
    return -1*np.sum(p*np.log(p))



def otu_name_biom(biom_row):
    """
    Given an OTU row from a BIOM table, determine a Genus-species identifier
    from the taxonomic specifier (see otu_name() method)
    """
    return otu_name(biom_row['metadata']['taxonomy'])
    
    
def otu_name(tax):
    """
    Determine a simple Genus-species identifier for an OTU, if possible. 
    If OTU is not identified to the species level, name it as 
    Unclassified (familly/genus/etc...).
    
    :@param tax: A list of qiime-style taxonomy identifiers, e.g. 
                 ['k__Bacteria', u'p__Firmicutes', u'c__Bacilli', ...
    """
    extract_name = lambda lvl: '_'.join(lvl.split('_')[2:])
    for i, lvl in enumerate(tax):
        lvl = lvl.strip()
        if i < len(tax) - 1 and len(tax[i + 1].strip()) == 3:
            if tax[i].strip()[0] == 'g':
                return extract_name(lvl) + '_spp.'
            else:
                return 'Unclassified_' + extract_name(lvl)
        elif i == len(tax) - 1:
            name = extract_name(lvl)
            if lvl[0] == 's':
                name = extract_name(tax[i-1]) + '_' + name
            return name

def load_core_file(core_fp):
    with open(core_fp, 'rU') as in_f:
        return {otu_name(ast.literal_eval(line.split('\t')[1])) for line in in_f.readlines()[1:]}


def assign_otu_membership(biom):
    """
    Determines the number of OTUs associated with samples using
    fuzzy sets with membership amount determined by relative abundance.
    """
    samples = defaultdict(FuzzySet)
    smax = calculate_total_abundance(biom)
    
    for row,col,amt in biom['data']:
        otu_name = otu_name_biom(biom['rows'][row])
        sampleID = biom['columns'][col]['id']
        
        samples[sampleID].add(otu_name, amt/smax[sampleID])

    return samples


def print_membership(entry):
    """
    Given an entry from a Sample dictionary (see assign_otu_membership) of fuzzy 
    otu membership, pretty-print the members as percentages.
    """
    data = [str(e).split(' \\ ') for e in entry]
    data = [(name, float(mu)) for name, mu in data]
    data = ['{0}: {1:4.2%}'.format(name, mu) for name, mu in sorted(data, key=lambda x:x[1])]
    for e in data:
        print e


def calculate_total_abundance(biom):
    """
    Calculates the total abundance for each sample ID
    """
    smax = defaultdict(int)
    for row,col,amt in biom['data']:
        sampleID = biom['columns'][col]['id']
        smax[sampleID] += amt
    return smax


def relative_abundance(biom):
    """
    Given a BIOM table, calculate per-sample relative abundance for 
    each OTU.
    :@type biom: dict (translated json string)
    :@param biom: BIOM-formatted OTU/Sample abundance data
    """
    rel_abundance = defaultdict(lambda: defaultdict(dict))
    smax = calculate_total_abundance(biom)
    
    for row,col,amt in biom['data']:
        otuName = otu_name_biom(biom['rows'][row])
        sampleID = biom['columns'][col]['id']

        rel_abundance[otuName][sampleID] = amt/smax[sampleID]
    
    return rel_abundance


def biom_summary(biom, id_match):
    tamtcounts = defaultdict(int)
    tot_seqs = 0.0
    
    for row, col, amt in biom['data']:
        sampleID = biom['columns'][col]['id']
        if not id_match in sampleID:
            continue
        tot_seqs += amt
        rtax = biom['rows'][row]['metadata']['taxonomy']
        for i,t in enumerate(rtax):
            t = t.strip()
            if i == len(rtax)-1 and len(t) > 3 and len(rtax[-1]) > 3:
                t = 's__'+rtax[i-1].strip().split('_')[-1]+'_'+t.split('_')[-1]
            tamtcounts[t] += amt        
    
    print 'Total sequence count: %i\n' % tot_seqs
    
    def levelData(counts, seqCount, level='g'):
        lvlData = sorted([(k, counts[k]) for k in counts if level+'_' in k], key=lambda x: x[1], reverse=True)
        return ['{}: {:4.2f}%'.format(p[0], p[1]/seqCount*100) for p in lvlData]
    
    def lookup(s, d):
        for idx, item in enumerate(d):
            if s in item:
                return d[idx]
        return s+': None'
    
    print 'Top families:\n'
    lvlDataF = levelData(tamtcounts, tot_seqs, 'f')
    print lvlDataF[:20]
            
    print '\nTop genera:\n'
    lvlDataG = levelData(tamtcounts, tot_seqs)
    print lvlDataG[:20]
            
    print '\nTop species:\n'
    lvlData = levelData(tamtcounts, tot_seqs, 's')
    print lvlData[:20]
    
    
    print '\nRed Complex:'
    print lookup('s__Porphyromonas_gingivalis', lvlData)
    print lookup('s__Treponema_denticola', lvlData)
    print lookup('s__Tannerella_forsythia', lvlData)
    
    print '\nOrange Complex: '
    print lookup('s__Fusobacterium_nucleatum', lvlData)
    print lookup('s__Prevotella_intermedia', lvlData)
    print lookup('s__Prevotella_nigrescens', lvlData)
    print lookup('s__Peptostreptococcus_micros', lvlData)
    
    print '\nYellow Complex: '
    print lookup('s__Streptococcus_sanguis', lvlData)
    print lookup('s__Streptococcus_oralis', lvlData)
    print lookup('s__Streptococcus_mitis', lvlData)
    print lookup('s__Streptococcus_gordonii', lvlData)
    print lookup('s__Streptococcus_intermedius', lvlData)
    
    print '\nGreen Complex: '
    print lookup('s__Aggregatibacter_actinomycetemcomitans', lvlData)
    print lookup('s__Campylobacter_concisus', lvlData)
    print lookup('s__Eikenella_corrodens', lvlData)
    print lookup('g__Capnocytophaga', lvlDataG)
    
    print '\nPurple Complex: '
    print lookup('s__Veillonella_parvula', lvlData)
    print lookup('s__Actinomyces_odontolyticus', lvlData)
    print lookup('s__Selenomonas_noxia', lvlData)
    print lookup('s__Actinomyces_naeslundii', lvlData)
