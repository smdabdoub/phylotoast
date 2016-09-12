#!/usr/bin/env python

import sys
import argparse
from phylotoast import otu_calc as otuc
importerrors = []
try:
    import biom
except ImportError:
    importerrors.append("biom")
try:
    import pandas as pd
except ImportError:
    importerrors.append("pandas")
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
                        help="Either a text file containing a list (one per line) "
                        "of OTU IDs, or a tab-separated (classic) BIOM-format file. "
                        "[REQUIRED]")
    parser.add_argument("-o", "--output_fp", default="converted_otus.txt",
                        help="For a list input, a new file containing a list of OTU IDs "
                              "and their corresponding short taxonomic identifiers "
                              "separated by tabs. For a BIOM file input, a new "
                              "mapping file with all the OTU IDs replaced by the short "
                              "identifier.")
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
        otu_fnh = pd.read_csv(args.otu_id_fp, sep="\n")
    except IOError as ioe:
        sys.exit("\nError with file containing OTUs:{}\n".format(ioe))
    else:
        otu_ids = otu_fnh.iloc[:, 0].values  # get only first column data (otus)

    output = []
    for val, idx, md in biomf.iter(axis="observation"):
        name = otuc.otu_name(md["taxonomy"])
        if args.reverse_lookup:
            if name in otu_ids:
                output.append(idx)   # Get otuids from otu names
        else:
            if idx in otu_ids:
                output.append(name)  # Get otu name from otu IDs
    name_S = pd.Series(output)
    name_S.to_csv(args.output_fp, sep="\n", index=False)


if __name__ == "__main__":
    sys.exit(main())
