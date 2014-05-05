#!/usr/bin/env python

from __future__ import division
import argparse
from collections import defaultdict
import json

def otu_name(biom_row):
    """
    Determine a simple Genus-species identifier for an OTU, if possible. 
    If OTU is not identified to the species level, name it as Unclassified (familly/genus/etc...)
    """
    tax = biom_row['metadata']['taxonomy']
    for i, lvl in enumerate(tax):
        lvl = lvl.strip()
        if i < len(tax) - 1 and len(tax[i + 1].strip()) == 3:
            if tax[i].strip()[0] == 'g':
                return lvl.split('_')[-1] + '_spp.'
            else:
                return 'Unclassified_' + lvl.split('_')[-1]
        elif i == len(tax) - 1:
            name = lvl.split('_')[-1]
            if lvl[0] == 's':
                name = tax[i-1].split('_')[-1] + '_' + name
            return name


def calculate_total_abundance(biom):
    """
    Calculates the total abundance for each sample ID
    """
    smax = defaultdict(int)
    for row,col,amt in biom['data']:
		otuID = biom['rows'][row]['id']
		sampleID = biom['columns'][col]['id']
		
		smax[sampleID] += amt
    return smax


def output_relative_abundance(biom, outFN):
	abundance = defaultdict(lambda: defaultdict(dict))
	sids = []
	smax = calculate_total_abundance(biom)
	
	for row,col,amt in biom['data']:
		otuName = otu_name(biom['rows'][row])
		sampleID = biom['columns'][col]['id']

		sids.append(sampleID)
		abundance[otuName][sampleID] = str(amt/smax[sampleID])

	with open(outFN, 'w') as outF:
		sids = sorted(set(sids), key=lambda x:x[-3:])
		sids.insert(0, 'OTU')
		outF.write('\t'.join(sids)+'\n')
		for row in biom['rows']:
			otuName = otu_name(row)
			outF.write('{}\t'.format(otuName))
			sabd = []
			for col in biom['columns']:
				sampleID = col['id']
				sabd.append(abundance[otuName][sampleID] if sampleID in abundance[otuName] else '0')
			outF.write('\t'.join(sabd)+'\n')


def handle_program_options():
    parser = argparse.ArgumentParser(description="Convert a BIOM file of OTU abundance \
                                                  data into a CSV of relative abundance \
                                                  data.")
    parser.add_argument('-i', '--input_biom_fp',
                        help="The BIOM file path.")
    parser.add_argument('-o', '--output_csv_fp', default='relative_abundance.csv',
                        help="A CSV table of relative OTU abundance data.")
    parser.add_argument('-v', '--verbose', action='store_true',
                        help="Print detailed information about script \
                              operation.")
    
    return parser.parse_args()


def main():
    args = handle_program_options()
    
    with open(args.input_biom_fp, 'rU') as inF:
        biom = json.loads(inF.readline())
    output_relative_abundance(biom, args.output_csv_fp)
    
if __name__ == '__main__':
    main()
