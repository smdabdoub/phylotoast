#!/usr/bin/env python
"""
Created on Feb 8, 2012

Author: Shareef M Dabdoub
"""
import sys
import re
import argparse
from phylotoast import biom_calc as bc, otu_calc as oc, util
try:
    import biom
except ImportError:
    sys.exit("Please install missing module: {}.".format("biom-format"))


def find_otu(otuid, tree):
    """
    Find an OTU ID in a Newick-format tree.
    Return the starting position of the ID or None if not found.
    """
    for m in re.finditer(otuid, tree):
        before, after = tree[m.start()-1], tree[m.start()+len(otuid)]
        if before in ["(", ",", ")"] and after in [":", ";"]:
            return m.start()
    return None


def newick_replace_otuids(tree, biomf):
    """
    Replace the OTU ids in the Newick phylogenetic tree format with truncated
    OTU names
    """
    for val, id_, md in biomf.iter(axis="observation"):
        otu_loc = find_otu(id_, tree)
        if otu_loc is not None:
            tree = tree[:otu_loc] + \
                   oc.otu_name(md["taxonomy"]) + \
                   tree[otu_loc + len(id_):]
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
    parser.add_argument("-i", "--otu_table", required=True,
                        help="The biom-format file with OTU-Sample abundance \
                              data.")
    parser.add_argument("-m", "--mapping", required=True,
                        help="The mapping file specifying group information \
                              for each sample.")
    parser.add_argument("-t", "--input_tree", default="",
                        help="A phylogenetic tree in Newick format to be \
                              modified by exchanging the OTU ID node names for\
                              taxonomic names.")
    parser.add_argument("-e", "--output_tre", default="iTol.tre",
                        help="The output .tre file")
    parser.add_argument("-o", "--output_itol_table", default="iTol_table.txt",
                        help="Other than a phylogenetic tree, the main input \
                              to iTol is a dataset file containing some \
                              representation of the abundance of every OTU \
                              across the specified data groups. This program \
                              provides multiple calculation methods. See the \
                              --analysis_metric option for details.")
    parser.add_argument("-c", "--map_categories", default=None,
                        help="Any mapping categories, such as treatment type, \
                              that will be used to group the data in the \
                              output iTol table. For example, one category \
                              with three types will result in three data \
                              columns in the final output. Two categories with\
                              three types each will result in six data \
                              columns. Default is no categories and all the \
                              data will be treated as a single group.")
    parser.add_argument("-a", "--analysis_metric", default="MRA",
                        choices=["MRA", "NMRA", "raw"],
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
    parser.add_argument("--stabilize_variance", action="store_true",
                        default=False,
                        help="Apply the variance-stabilizing arcsine square\
                              root transformation to the OTU proportion data.")

    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.otu_table):
            pass
    except IOError as ioe:
        sys.exit(
            "\nError with OTU_Sample abundance data file:{}\n"
            .format(ioe)
        )

    try:
        with open(args.mapping):
            pass
    except IOError as ioe:
        sys.exit(
            "\nError with mapping file:{}\n"
            .format(ioe)
        )

    # input data
    biomf = biom.load_table(args.otu_table)
    map_header, imap = util.parse_map_file(args.mapping)

    # rewrite tree file with otu names
    if args.input_tree:
        with open(args.input_tree) as treF, open(args.output_tre, "w") as outF:
            tree = treF.readline()
            if "'" in tree:
                tree = tree.replace("'", '')
            outF.write(newick_replace_otuids(tree, biomf))

    oid_rows = {id_: md["taxonomy"]
                for val, id_, md in biomf.iter(axis="observation")}

    # calculate analysis results
    categories = None
    if args.map_categories is not None and args.analysis_metric != "raw":
        categories = args.map_categories.split(",")

    # set transform if --stabilize_variance is specfied
    tform = bc.arcsine_sqrt_transform if args.stabilize_variance else None

    groups = util.gather_categories(imap, map_header, categories)
    for group in groups.values():
        if args.analysis_metric in ["MRA", "NMRA"]:
            results = bc.MRA(biomf, group.sids, transform=tform)
        elif args.analysis_metric == "raw":
            results = bc.transform_raw_abundance(biomf, sampleIDs=group.sids,
                                                 sample_abd=False)
        group.results.update({oc.otu_name(oid_rows[oid]): results[oid]
                             for oid in results})

    # write iTol data set file
    with open(args.output_itol_table, "w") as itolF:
        if args.analysis_metric == "raw":
            itolF.write("DATASET_GRADIENT\nSEPARATOR TAB\n")
            itolF.write("DATASET_LABEL\tLog Total Abundance\n")
            itolF.write("COLOR\t#000000\n")
            itolF.write("LEGEND_TITLE\tLog Total Abundance\n")
            itolF.write("LEGEND_SHAPES\t1\n")
            itolF.write("LEGEND_COLORS\t#000000\n")
            itolF.write("LEGEND_LABELS\tLog Total Abundance\n")
            itolF.write("COLOR_MIN\t#FFFFFF\n")
            itolF.write("COLOR_MAX\t#000000\n")
        else:
            itolF.write("DATASET_MULTIBAR\nSEPARATOR TAB\n")
            itolF.write("DATASET_LABEL\tNMRA\n")
            itolF.write("FIELD_COLORS\t{}\n".format("\t".join(["#ff0000"
                        for _ in range(len(groups))])))
            itolF.write("FIELD_LABELS\t" + "\t".join(groups.keys())+"\n")
            itolF.write("LEGEND_TITLE\tNMRA\n")
            itolF.write("LEGEND_SHAPES\t{}\n".format("\t".join(["1"
                        for _ in range(len(groups))])))
            itolF.write("LEGEND_COLORS\t{}\n".format("\t".join(["#ff0000"
                        for _ in range(len(groups))])))
            itolF.write("LEGEND_LABELS\t" + "\t".join(groups.keys())+"\n")
            itolF.write("WIDTH\t300\n")
        itolF.write("DATA\n")
        all_otus = frozenset({oc.otu_name(md["taxonomy"])
                              for val, id_, md in
                              biomf.iter(axis="observation")})

        for oname in all_otus:
            row = ["{name}"]        # \t{s:.2f}\t{ns:.2f}\n"
            row_data = {"name": oname}
            msum = 0
            for name, group in groups.iteritems():
                row.append("{{{}:.5f}}".format(name))
                if oname in group.results:
                    row_data[name] = group.results[oname]
                else:
                    row_data[name] = 0.0
                msum += row_data[name]
            # normalize avg relative abundance data
            if args.analysis_metric == "NMRA" and msum > 0:
                row_data.update({key: data/msum
                                for key, data in row_data.items()
                                if key != "name"})
            itolF.write("\t".join(row).format(**row_data) + "\n")

if __name__ == "__main__":
    main()
