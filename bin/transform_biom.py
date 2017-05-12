#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Apply various transforms to data in a BIOM-format
table.
"""
from __future__ import absolute_import, division, print_function

# standard library imports
import argparse
from gzip import open as gzip_open
# 3rd party imports
import biom
import numpy as np

try:
    import h5py
    HAVE_H5PY = True
except ImportError:
    HAVE_H5PY = False

# local imports
import phylotoast


def write_biom(biom_tbl, output_fp, fmt="hdf5", gzip=False):
    """
    Write the BIOM table to a file.
    :type biom_tbl: biom.table.Table
    :param biom_tbl: A BIOM table containing the per-sample OTU counts and metadata
                  to be written out to file.
    :type output_fp str
    :param output_fp: Path to the BIOM-format file that will be written.
    :type fmt: str
    :param fmt: One of: hdf5, json, tsv. The BIOM version the table will be
                output (2.x, 1.0, 'classic').
    """
    opener = open
    mode = 'w'
    if gzip and fmt != "hdf5":
        if not output_fp.endswith(".gz"):
            output_fp += ".gz"
        opener = gzip_open
        mode = 'wt'

    # HDF5 BIOM files are gzipped by default
    if fmt == "hdf5":
        opener = h5py.File

    gen_str = "PhyloToAST v{} (phylotoast.org)".format(phylotoast.__version__)
    biom_tbl.generated_by = gen_str

    with opener(output_fp, mode) as biom_f:
        if fmt == "json":
            biom_tbl.to_json(biom_tbl.generated_by, direct_io=biom_f)
        elif fmt == "tsv":
            biom_f.write(biom_tbl.to_tsv())
        else:
            biom_tbl.to_hdf5(biom_f, biom_tbl.generated_by)

    return output_fp


def relative_abd(biom_tbl):
    return biom_tbl.norm(inplace=False)


def log10(biom_tbl):
    log10_tfm = lambda data, id_, md: np.nan_to_num(np.log10(data))
    return biom_tbl.transform(log10_tfm, inplace=False)


def relative_abd_log10(biom_tbl):
    tbl_norm = relative_abd(biom_tbl)
    return log10(tbl_norm)


def arcsin_sqrt(biom_tbl):
    """
    Applies the arcsine square root transform to the
    given BIOM-format table
    """
    arcsint = lambda data, id_, md: np.arcsin(np.sqrt(data))

    tbl_relabd = relative_abd(biom_tbl)
    tbl_asin = tbl_relabd.transform(arcsint, inplace=False)

    return tbl_asin


transforms = {"arcsin_sqrt": arcsin_sqrt, 
              "ra": relative_abd, 
              "log10": log10, 
              "ra_log10": relative_abd_log10}

def handle_program_options():
    """Parses the given options passed in at the command line."""
    parser = argparse.ArgumentParser(description="This script applies various"
                                     "transforms to the data in a given "
                                     "BIOM-format table and outputs a new"
                                     "BIOM table with the transformed data.")
    parser.add_argument("-i", "--biom_table_fp", required=True,
                        help="Path to the input BIOM-format table. [REQUIRED]")
    parser.add_argument("-t", "--transform", default="arcsin_sqrt", 
                        choices=transforms.keys(),
                        help="The transform to apply to the data. Default: "
                             "arcsine square root.")
    parser.add_argument('--fmt', default="hdf5", 
                        choices=["hdf5", "json", "tsv"],
                        help="Set the output format of the BIOM table.\
                              Default is HDF5.")
    parser.add_argument('--gzip', action='store_true',
                        help="Compress the output BIOM table with gzip.\
                              HDF5 BIOM (v2.x) files are internally\
                              compressed by default, so this option\
                              is not needed when specifying --fmt hdf5.")
    parser.add_argument("-o", "--output_fp", required=True,
                        help="Output path for the transformed BIOM table."
                             "[REQUIRED]")
    parser.add_argument('-v', '--verbose', action='store_true')


    return parser.parse_args()


def main():
    args = handle_program_options()

    biom_tbl = biom.load_table(args.biom_table_fp)
    tbl_tform = transforms[args.transform](biom_tbl)

    write_biom(tbl_tform, args.output_fp, fmt=args.fmt, gzip=args.gzip)

    if args.verbose:
        print("Transformed table written to: {}".format(args.output_fp))



if __name__ == '__main__':
    main()