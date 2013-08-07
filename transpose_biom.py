#!/usr/bin/env python

'''
Created on Jul 1, 2013

@author: Shareef Dabdoub
'''
import argparse
import copy
from collections import defaultdict
import json
import os.path as osp
import sys
import util

def split_by_category(biom_cols, mapping, category_id):
    """
    Split up the column data in a biom table by mapping category value.
    """
    columns = defaultdict(list)
    for i,col in enumerate(biom_cols):
        columns[mapping[col['id']][category_id]].append((i,col))
    
    return columns

def handle_program_options():
    parser = argparse.ArgumentParser(description="Transpose a BIOM-format file\
                                     so that the matrix is sample by species.")
    parser.add_argument('-i','--input_biom_fp', required=True, 
                        help="The BIOM-format file.")
    parser.add_argument('-m', '--mapping', required=True,
                        help="The mapping file specifying group information \
                              for each sample.")
    parser.add_argument('-c', '--map_category', default=None,
                        help="A mapping category, such as TreatmentType, \
                              that will be used to split the data into \
                              separate BIOM files; one for each value found\
                              in the category.")
    parser.add_argument('-o','--output_biom_fp', default='transposed.biom', 
                        required=True, help="The BIOM-format file to write.")

    parser.add_argument('-v', '--verbose', action='store_true')
    
    return parser.parse_args()

def main():
    args = handle_program_options()
    out_fp, ext = osp.splitext(args.output_biom_fp)

    with open(args.input_biom_fp) as bF:
        biom = json.loads(bF.readline())
    
    mapping = util.parse_map_file(args.mapping)
    with open(args.mapping, 'rU') as mF:
        header = mF.readline().strip().split('\t')
    
    try:
        category_id = header.index(args.map_category)
    except ValueError:
        sys.exit('Category {} not found in supplied mapping file.'.format(args.map_category))
        
    values = {mapping[sid][category_id] for sid in mapping}
    biom_copies = {value: copy.deepcopy(biom) for value in values}
    split_samples = split_by_category(biom['columns'], mapping, category_id)
    for cat_val in biom_copies:
        biom_copies[cat_val]['data'] = []
        biom_copies[cat_val]['rows'], biom_copies[cat_val]['columns'] = [item[1] for item in split_samples[cat_val]], biom_copies[cat_val]['rows']
        sample_ids = [item[0] for item in split_samples[cat_val]]
        
        for i in xrange(len(biom['data'])):
            if biom['data'][i][1] in sample_ids:
                row, col, amt = biom['data'][i]
                biom_copies[cat_val]['data'].append([sample_ids.index(col), row, amt])
        
        biom_copies[cat_val]['shape'] = [len(biom_copies[cat_val]['rows']),
                                         len(biom_copies[cat_val]['columns'])]

        with open(out_fp + '_' + cat_val + ext, 'w') as outF:
            outF.write(json.dumps(biom_copies[cat_val]))


if __name__ == '__main__':
    main()



















