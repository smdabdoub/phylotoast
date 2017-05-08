#!/usr/bin/env python

import sys
import csv
import argparse
from phylotoast import otu_calc as otuc
importerrors = []
try:
    import biom
except ImportError:
    importerrors.append("biom")
if len(importerrors) > 0:
    for err in importerrors:
        print("Please install missing module: {}".format(err))
    sys.exit()


def handle_program_options():
    parser = argparse.ArgumentParser(description="Convert a list of OTU IDs to a list of "
                                                 "OTU IDs paired with Genus_Species "
                                                 "identifiers and perform reverse lookup,"
                                                 " if needed.")
    parser.add_argument("-i", "--otu_table", required=True,
                        help="Input biom file format OTU table. [REQUIRED]")
    parser.add_argument("-oid", "--otu_id_fp", required=True,
                        help="A single or multi-column (tab-separated) file containing "
                        "the OTU IDs to be converted in the first column. [REQUIRED]")
    parser.add_argument("-o", "--output_fp",
                        help="For a list input, a new file containing a list of OTU IDs "
                              "and their corresponding short taxonomic identifiers "
                              "separated by tabs. For a BIOM file input, a new "
                              "mapping file with all the OTU IDs replaced by the short "
                              "identifier.[REQUIRED]")
    parser.add_argument("--reverse_lookup", action="store_true",
                        help="Get OTUIDs from genus-species OTU name.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Read biom format file
    try:
        biomf = biom.load_table(args.otu_table)
    except IOError as ioe:
        sys.exit("\nError with biom format file (-i): {}\n".format(ioe))

    # Read in otus file data
    try:
        with open(args.otu_id_fp, "rU") as inf:
            csvr = csv.reader(inf, delimiter="\t")
            otu_ids = [line[0] for line in csvr]
    except IOError as ioe:
        sys.exit("\nError with file containing OTUs:{}\n".format(ioe))

    output = {}
    for val, idx, md in biomf.iter(axis="observation"):
        name = otuc.otu_name(md["taxonomy"])
        if args.reverse_lookup:
            if name in otu_ids:
                output[name] = idx   # Get otuids from otu names
        else:
            if idx in otu_ids:
                output[idx] = name  # Get otu name from otu IDs
    with open(args.output_fp, "w") as outf:
        outf.write("Input\tOutput\n")
        for k, v in output.items():
            outf.write("{0}\t{1}\n".format(k, v))


if __name__ == "__main__":
    sys.exit(main())
