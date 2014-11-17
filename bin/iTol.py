#!/usr/bin/env python
'''
Created on Feb 8, 2012

@author: Shareef M Dabdoub
'''
from qiime_tools import biom_calc as bc, otu_calc as oc, util
import argparse
from collections import namedtuple, OrderedDict
import json
import re


def find_otu(otuid, tree):
    """
    Find an OTU ID in a Newick-format tree.
    Return the starting position of the ID or None if not found.
    """
    for m in re.finditer(otuid, tree):
        before, after = tree[m.start()-1], tree[m.start()+len(otuid)]
        if before in ['(', ',', ')'] and after in [':', ';']:
            return m.start()
    return None


def newick_replace_otuids(tree, biom):
    """
    Replace the OTU ids in the Newick phylogenetic tree format with truncated
    OTU names
    """
    for row in biom['rows']:
        otu_loc = find_otu(row['id'], tree)
        if otu_loc is not None:
            tree = tree[:otu_loc] + oc.otu_name_biom(row) + tree[otu_loc+len(row['id']):]
        else:
            print 'ID not found:', row['id']
    return tree


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
    parser.add_argument('-a', '--analysis_metric', default='MRA',
                        choices=['MRA', 'NMRA', 'raw'],
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
    map_header, imap = util.parse_map_file(args.mapping)

    # rewrite tree file with otu names
    if args.input_tree:
        with open(args.input_tree) as treF, open(args.output_tre, 'w') as outF:
            tree = treF.readline()
            if "'" in tree:
                tree = tree.replace("'", '')
            outF.write(newick_replace_otuids(tree, biom))

    oid_rows = {row['id']: row for row in biom['rows']}

    # calculate analysis results
    categories = None
    if args.map_categories is not None:
        categories = args.map_categories.split(',')

    groups = util.gather_categories(imap, map_header, categories)
    for group in groups.values():
        if args.analysis_metric in ['MRA', 'NMRA']:
            results = bc.MRA(biom, group.sids)
        elif args.analysis_metric == 'raw':
            results = bc.transform_raw_abundance(biom, sampleIDs=group.sids,
                                                 sample_abd=False)

        group.results.update({oc.otu_name_biom(oid_rows[oid]): results[oid]
                             for oid in results})

    # write iTol data set file
    with open(args.output_itol_table, 'w') as itolF:
        itolF.write('LABELS\t' + '\t'.join(groups.keys())+'\n')
        itolF.write('COLORS\t{}\n'.format('\t'.join(['#ff0000'
                    for _ in range(len(groups))])))
        all_otus = frozenset({oc.otu_name_biom(row) for row in biom['rows']})

        for oname in all_otus:
            row = ['{name}']        # \t{s:.2f}\t{ns:.2f}\n'
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
                row_data.update({key: data/msum
                                for key, data in row_data.items()
                                if key != 'name'})

            itolF.write('\t'.join(row).format(**row_data) + '\n')

if __name__ == '__main__':
    main()
