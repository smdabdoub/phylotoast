#!/usr/bin/env python
'''
Created on Feb 8, 2012

@author: Shareef M Dabdoub
'''
import biom_calc as bc
import util
import argparse
from collections import namedtuple, OrderedDict
import json


def otu_name(biom_row):
    tax = biom_row['metadata']['taxonomy']
    for i, lvl in enumerate(tax):
        lvl = lvl.strip()
        if i < len(tax) - 1 and len(tax[i + 1].strip()) == 3:
            return 'Unclassified_' + lvl.split('_')[-1]
        elif i == len(tax) - 1:
            name = lvl.split('_')[-1]
            if lvl[0] == 's':
                name = tax[i-1].split('_')[-1] + '_' + name 
            return name


def newick_replace_otuids(tree, biom):
    """
    Replace the OTU ids in the Newick phylogenetic tree format with truncated 
    OTU names 
    """
    for row in biom['rows']:
        tree = tree.replace(row['id'], otu_name(row))
    return tree


# Meant to contain all the data necessary for calculating a single column of 
# an iTol data table
DataCategory = namedtuple('DataCategory', 'sids results')

def gather_categories(imap, header, categories=None):
    """
    Find the user specified categories in the map and create a dictionary 
    to contain the relevant data for each type within the categories. Multiple
    categories will have their types combined such that each possible 
    combination will have its own entry in the dictionary.
    
    :@type imap: dict
    :@param imap: The input mapping file data keyed by SampleID
    :@type header: list
    :@param header: The header line from the input mapping file. This will be 
                    searched for the user-specified categories
    :@type categories: list
    :@param categories: The list of user-specified categories from the mapping 
                        file
    :@rtype: OrderedDict
    :@return: A sorted dictionary keyed on the combinations of all the types 
    		  found within the user-specified categories. Each entry will 
    		  contain an empty DataCategory namedtuple. If no categories are 
    		  specified, a single entry with the key 'default' will be returned  
    """
    if categories is None:
        return {'default': DataCategory(set(imap.keys()), {})}
    
    cat_ids = [header.index(cat) for cat in categories if cat in header
                                                         and '=' not in cat]
    conditions = {}
    for cat in categories:
        if '=' in cat and cat.split('=')[0] in header:
            conditions[header.index(cat.split('=')[0])] =  cat.split('=')[1]

    if not cat_ids and not conditions:
        return {'default': DataCategory(set(imap.keys()), {})}
    
    if not cat_ids and conditions:
        sids = {sid for sid,row in imap.iteritems() if 
                            all([row[c] == conditions[c] for c in conditions])}
        return {'default': DataCategory(sids, {})}
    
    table = {}
    for sid,row in imap.iteritems():
        if all([row[c] == conditions[c] for c in conditions]):
            key = '_'.join([row[cid] for cid in cat_ids])
            if not key in table:
                table[key] = DataCategory(set(), {})
            table[key].sids.add(sid)
            
    return OrderedDict(sorted(table.items(), key=lambda t: t[0]))


def handle_program_options():
    parser = argparse.ArgumentParser(description="Create files appropriate for\
                                     use in the iTol visualization program by \
                                     using the abundance data from a \
                                     biom-format file and groups specified in \
                                     a QIIME mapping file. The program also \
                                     modifies a Newick-format phylogenetic \
                                     tree file to use proper taxonomic names \
                                     instead of OTU IDs for useful display in \
                                     iTol.")
    parser.add_argument('-i', '--otu_table', required=True,
                        help="The biom-format file with OTU-Sample abundance \
                              data.")
    parser.add_argument('-m', '--mapping', required=True,
                        help="The mapping file specifying group information \
                              for each sample.")
    parser.add_argument('-t', '--input_tree', default='',
                        help="A phylogenetic tree in Newick format to be \
                              modified by exchanging the OTU ID node names for\
                              taxonomic names.")
    parser.add_argument('-e', '--output_tre', default='iTol.tre',
                        help="The output .tre file")
    parser.add_argument('-o', '--output_itol_table', default='iTol_table.txt',
                        help="Other than a phylogenetic tree, the main input \
                              to iTol is a dataset file containing some \
                              representation of the abundance of every OTU \
                              across the specified data groups. This program \
                              provides multiple calculation methods. See the \
                              --analysis_metric option for details.")
    parser.add_argument('-c', '--map_categories', default=None,
                        help="Any mapping categories, such as treatment type, \
                              that will be used to group the data in the \
                              output iTol table. For example, one category \
                              with three types will result in three data \
                              columns in the final output. Two categories with\
                              three types each will result in six data \
                              columns. Default is no categories and all the \
                              data will be treated as a single group.")
    parser.add_argument('-a', '--analysis_metric', default='MRE', 
                        choices=['MRA','NMRA','raw'],
                        help="Specifies which metric is calculated on the \
                              abundance data in the OTU table. Available \
                              options: MRE - mean relative abundance \
                              (Abundance data is normalized by total sample \
                              abundance, then averaged across OTU), NMRE - \
                              normalized mean relative abundance (MRE \
                              normalized by the total MRE across the groups \
                              as specified in --map_categories), raw (outputs \
                              the actual sequence abundance data for \
                              each OTU).")
#    parser.add_argument('-v', '--verbose', action='store_true')
    
    return parser.parse_args()


def main():
    args = handle_program_options()
    
    # input data
    with open(args.otu_table) as bF:
        biom = json.loads(bF.readline())
    imap = util.parse_map_file(args.mapping)
    with open(args.mapping, 'rU') as mapF:
        map_header = mapF.readline()[1:].split()
        
    # rewrite tree file with otu names
    if args.input_tree:
        with open(args.input_tree) as treF, open(args.output_tre, 'w') as outF:
            outF.write(newick_replace_otuids(treF.readline(), biom))
    
    oid_rows = {row['id']:row for row in biom['rows']}
    
    # calculate analysis results
    categories = None
    if args.map_categories is not None:
        categories = args.map_categories.split(',')
    
    groups = gather_categories(imap, map_header, categories)
    for group in groups.values():
        if args.analysis_metric in ['MRA','NMRA']:
            results = bc.MRA(biom, group.sids)
        elif args.analysis_metric == 'raw':
            results = bc.raw_abundance(biom, group.sids)

        group.results.update({otu_name(oid_rows[oid]):results[oid] 
                                    for oid in results})
        
    # write iTol data set file
    with open('iTol_table.txt', 'w') as itolF:
        itolF.write('LABELS\t' + '\t'.join(groups.keys())+'\n')
        itolF.write('COLORS\t{}\n'.format('\t'.join(['#ff0000' 
											    for _ in range(len(groups))])))
        all_otus = frozenset({otu_name(row) for row in biom['rows']})
        
        for oname in all_otus:
            row = ['{name}']#\t{s:.2f}\t{ns:.2f}\n'
            row_data = {'name': oname}
            msum = 0
            for name, group in groups.iteritems():
                row.append('{{{}:.5f}}'.format(name))
                if oname in group.results:
                    row_data[name] = group.results[oname] 
                else:
                    row_data[name] = 0.0
                msum += row_data[name]
            # normalize avg relative abundance data
            if args.analysis_metric == 'NMRA' and msum > 0:
                row_data.update({key:data/msum for key,data in row_data.items() 
                                   if key != 'name'})
            
            itolF.write('\t'.join(row).format(**row_data) + '\n')
            
if __name__ == '__main__':
    main()
