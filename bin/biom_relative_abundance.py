#!/usr/bin/env python
"""
Given a BIOM table, calculate per-sample relative abundance for
each OTU and write out to a tab-separated file listing OTUs as
rows and Samples as columns.
"""
from __future__ import division
import sys
import biom
import argparse
from phylotoast import biom_calc as bc
from phylotoast import otu_calc as oc


def write_relative_abundance(rel_abd, biomf, out_fn, sort_by=None):
    """
    Given a BIOM table, calculate per-sample relative abundance for
    each OTU and write out to a tab-separated file listing OTUs as
    rows and Samples as columns.
    :type biom: biom object
    :param biom: BIOM-formatted OTU/Sample abundance data
    :type out_fn: str
    :param out_fn: The full path to the desired output file.
    :type sort_by: function
    :param sort_by: A function acting as a sorting key that will determine
                     the order in which the Sample IDs appear as columns in
                     the output file.
    """
    with open(out_fn, 'w') as out_f:
        sids = sorted(set(biomf.ids()), key=sort_by)
        out_f.write('#OTU ID\t{}\n'.format('\t'.join(sids)))

        for otuid in biomf.ids(axis="observation"):
            otuName = oc.otu_name(biomf.metadata(otuid, "observation")\
                      ["taxonomy"])
            sabd = [str(rel_abd[sid][otuid]) 
                    if sid in rel_abd and otuid in rel_abd[sid] else '0' 
                    for sid in sids]
            out_f.write('{}\t{}\n'.format(otuName, '\t'.join(sabd)))


def handle_program_options():
    parser = argparse.ArgumentParser(description="Convert a BIOM file of OTU abundance \
                                                  data into a TSV file of relative \
                                                  abundance data.")
    parser.add_argument('-i', '--input_biom_fp',
                        help="BIOM file path.")
    parser.add_argument('-o', '--output_tsv_fp', default='relative_abundance.tsv',
                        help="A TSV table of relative OTU abundance data.")
    parser.add_argument('--stabilize_variance', action='store_true',
                        help="Apply the variance-stabilizing arcsine square\
                              root transformation to the OTU proportion data.")

    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.input_biom_fp):
            pass
    except IOError as ioe:
        sys.exit('\nError in BIOM file path: {}\n'.format(ioe))

#     with open(args.input_biom_fp, 'rU') as in_f:
    biomf = biom.load_table(args.input_biom_fp)
    norm_biomf = biomf.norm(axis="observation", inplace=False)
    rel_abd = bc.relative_abundance(norm_biomf)
    if args.stabilize_variance:
        rel_abd = bc.arcsine_sqrt_transform(rel_abd)

    write_relative_abundance(rel_abd, norm_biomf, args.output_tsv_fp)

if __name__ == '__main__':
    main()
