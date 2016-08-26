#!/usr/bin/env python
"""
Abstract: Filter biom file on both 'sample' and 'observation' axes, given a
          list of sampleIDs to retain.
Author: Akshay Paropkari
Date: 02/15/2016
"""
import sys
import argparse
importerrors = []
try:
    import biom
    from biom.util import biom_open as bo
except ImportError as ie:
    importerrors.append(ie)
try:
    import pandas as pd
except ImportError as ie:
    importerrors.append(ie)
if len(importerrors) > 0:
    for err in importerrors:
        print("Import Error: {}".format(err))
    sys.exit()


def handle_program_options():
    parser = argparse.ArgumentParser(description="Filter biom file on both \
                                                 'sample' and 'observation' \
                                                 axes, given a list of \
                                                 sampieIDs to retain.")
    parser.add_argument("input_biom_fnh", help="BIOM file path.")
    parser.add_argument("output_biom_fnh", default="filtered_otu_table.biom",
                        help="Filtered biom output file.")
    parser.add_argument("mapping_fnh",
                        help="Mapping file with sampleIDs to retain in it. The"
                             " '#SampleID' column will be used to get the list"
                             " of all ids to retain.")
    parser.add_argument("-fo", "--filter_otuids_fnh", default="./otuids_to_discard.txt",
                        help="Path to file to write out the list of OTUIDs not present "
                             "in any SampleIDs in mapping file. This output is usually "
                             "used to filter out unwanted otuids from .tre file. If not "
                             "given, the discarded OTUIDs list will be saved in the "
                             "current working directory.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Error check input file
    try:
        biomf = biom.load_table(args.input_biom_fnh)
    except IOError as ioe:
        sys.exit("\nError in BIOM file path: {}\n".format(ioe))
    try:
        mapf = pd.read_csv(args.mapping_fnh, sep="\t")
    except IOError as ioe:
        sys.exit("\nError in mapping file path: {}\n".format(ioe))

    # Get filtered biom data
    keep_sampleIDs = list(mapf["#SampleID"])  # because biom.filter() expects a list of ids
    sid_filtered_biomf = biomf.filter(keep_sampleIDs, inplace=False)
    print("\n{} sampleIDs retained from original biom file.".format(len(keep_sampleIDs)))
    obs_abd_sum = sid_filtered_biomf.sum(axis="observation")
    otuids = sid_filtered_biomf.ids("observation")
    abd_sum = {a: b for a, b in zip(otuids, obs_abd_sum)}

    # Run checks for redundant otus
    redundant_otuids = [otu for otu, abd in abd_sum.items() if abd == 0]
    otuid_filtered_biomf = sid_filtered_biomf.filter(redundant_otuids,
                                                     "observation",
                                                     invert=True,
                                                     inplace=False)
    print("{} otuIDs filtered out of the original biom file.\n"
          .format(len(redundant_otuids)))

    # Write out files
    with bo(args.output_biom_fnh, "w") as rth:
        otuid_filtered_biomf.to_hdf5(rth, "Filtered OTU Table.")
    with open(args.filter_otuids_fnh, "w") as yui:
        for otuid in redundant_otuids:
            yui.write("{}\n".format(otuid))


if __name__ == "__main__":
    sys.exit(main())
