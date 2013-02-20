'''
Created on Feb 19, 2013

@author: Shareef Dabdoub

This module provides methods for calculating various metrics with regards to 
each OTU in an input OTU abundance table. This is currently used by iTol.py 
to offload the different methods.
'''
from collections import defaultdict

def relative_abundance(biom, sampleIDs):
    ra = {item['id']: defaultdict(int) for item in biom['columns'] 
           if item['id'] in sampleIDs}
    totals = defaultdict(float)
    
    for row, col, amt in biom['data']:
        otuID = biom['rows'][row]['id']
        sampleID = biom['columns'][col]['id']
        
        if sampleID in sampleIDs:
            ra[sampleID][otuID] = amt
            totals[sampleID] += amt
    
    return {sid: {oid: ra[sid][oid] / totals[sid] for oid in ra[sid]} 
              for sid in ra}


def mean_otu_pct_abundance(ra, otuIDs):
    sids = ra.keys()
    otumeans = defaultdict(int)
    
    for oid in otuIDs:
        otumeans[oid] = sum([ra[sid][oid] for sid in sids 
                               if oid in ra[sid]]) / len(sids) * 100
    
    return otumeans

def MRA(biom, sampleIDs):
    """
    Calculate the mean relative abundance.
    """
    ra = relative_abundance(biom, sampleIDs)
    otuIDs = {row['id'] for row in biom['rows']}
    return mean_otu_pct_abundance(ra, otuIDs)


def raw_abundance(biom, sampleIDs):
    """
    Calculate the total number of sequences in each OTU
    """
    results = defaultdict(int)
    for row, col, amt in biom['data']:
        otuID = biom['rows'][row]['id']
        sampleID = biom['columns'][col]['id']
        
        if sampleID in sampleIDs:
            results[otuID] += amt
            
    return results





























