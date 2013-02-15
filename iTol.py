# coding: utf-8
import util
from collections import defaultdict
import json

def otu_name(biom_row):
	tax = biom_row['metadata']['taxonomy']
	for i,lvl in enumerate(tax):
		if i < len(tax)-1 and len(tax[i+1].strip()) == 3:
			return 'Unclassified_'+lvl.split('_')[-1]
		elif i == len(tax)-1:
			return lvl.split('_')[-1]


def relative_abundance(biom, filterIDs):
	ra = {item['id']: defaultdict(int) for item in biom['columns'] 
	       if item['id'] in filterIDs}
	totals = defaultdict(float)
	
	for row,col,amt in biom['data']:
		oname = otu_name(biom['rows'][row])
		sampleID = biom['columns'][col]['id']
		
		if sampleID in filterIDs:
			ra[sampleID][oname] = amt
			totals[sampleID] += amt
	
	return {sid: {oid: ra[sid][oid]/totals[sid] for oid in ra[sid]} for sid in ra}


def mean_otu_pct_abundance(ra, otuids):
	sids = ra.keys()
	otumeans = defaultdict(int)
	
	for oid in otuids:
		otumeans[oid] = sum([ra[sid][oid] for sid in sids if oid in ra[sid]])/len(sids)
	
	return otumeans

def newick_replace_otuids(tree, biom):
	"""
	Replace the OTU ids in the Newick phylogenetic tree format with truncated 
	OTU names 
	"""
	for row in biom['rows']:
		tree = tree.replace(row['id'], otu_name(row))
	return tree

with open('data/genus_condensed_pruned_otu_table.biom') as bF:
	biom = json.loads(bF.readline())
with open('data/genus_condensed_pruned_rep_set.tre') as treF, open('data/genus_condensed_pruned_rep_set_itol.tre','w') as ntreF:
	ntreF.write(newick_replace_otuids(treF.readline(), biom))
hmap = util.parse_map_file('data/HealthySmokersMapping.txt')

# remove tooth data samples
sids = frozenset({sid for sid in hmap if hmap[sid][3] == 'Smoker'})
nsids = frozenset({sid for sid in hmap if hmap[sid][3] != 'Smoker'})

sra = relative_abundance(biom, sids)
nsra = relative_abundance(biom, nsids)

all_otus = {otu_name(item) for item in biom['rows']}

s_otu_means = mean_otu_pct_abundance(sra, all_otus)
ns_otu_means = mean_otu_pct_abundance(nsra, all_otus)


with open('iTol_table.txt', 'w') as itolF:
	itolF.write('LABELS\tSmokers\tNon-Smokers\n')
	itolF.write('COLORS\t#ff0000\t#00ff00\n')
	
	for oname in all_otus:
		s_avg = s_otu_means[oname]*100 if oname in s_otu_means else 0.0
		ns_avg = ns_otu_means[oname]*100 if oname in ns_otu_means else 0.0
		row = '{name}\t{s:.2f}\t{ns:.2f}\n'
		itolF.write(row.format(name=oname, s=s_avg, ns=ns_avg))

