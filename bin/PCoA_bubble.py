#!/usr/bin/env python
"""
Create a series of PCoA plots where the marker size varies by relative
abundance of a particular OTU

Author: Shareef M Dabdoub
"""
from __future__ import division
import os
import sys
import biom
import argparse
try:
    import matplotlib.pyplot as plt
except ImportError as ie:
    sys.exit("Import Error. Please install missing module: {}".format(ie))
from phylotoast import util, graph_util as gu, biom_calc as bc, otu_calc as oc


def calculate_xy_range(data):
    xr = [float("inf"), float("-inf")]
    yr = [float("inf"), float("-inf")]

    for cat in data:
        pc1, pc2 = data[cat]["pc1"], data[cat]["pc2"]
        if pc1:
            xr[0] = min(min(pc1), xr[0])
            xr[1] = max(max(pc1), xr[1])
        if pc2:
            yr[0] = min(min(pc2), yr[0])
            yr[1] = max(max(pc2), yr[1])

    return xr, yr


def plot_PCoA(cat_data, otu_name, unifrac, names, colors, xr, yr, outDir,
              save_as, plot_style):
    """
    Plot PCoA principal coordinates scaled by the relative abundances of
    otu_name.
    """
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(111)

    for i, cat in enumerate(cat_data):
        plt.scatter(cat_data[cat]["pc1"], cat_data[cat]["pc2"], cat_data[cat]["size"],
                    color=colors[cat], alpha=0.85, marker="o", edgecolor="black",
                    label=cat)
    lgnd = plt.legend(loc="best", scatterpoints=3, fontsize=13)
    for i in range(len(colors.keys())):
        lgnd.legendHandles[i]._sizes = [80]  # Change the legend marker size manually
    plt.title(" ".join(otu_name.split("_")), style="italic")
    plt.ylabel("PC2 (Percent Explained Variance {:.3f}%)".format(float(unifrac["varexp"][1])))
    plt.xlabel("PC1 (Percent Explained Variance {:.3f}%)".format(float(unifrac["varexp"][0])))
    plt.xlim(round(xr[0]*1.5, 1), round(xr[1]*1.5, 1))
    plt.ylim(round(yr[0]*1.5, 1), round(yr[1]*1.5, 1))
    if plot_style:
        gu.ggplot2_style(ax)
        fc = "0.8"
    else:
        fc = "none"
    fig.savefig(os.path.join(outDir, "_".join(otu_name.split())) + "." + save_as,
                facecolor=fc, edgecolor="none", format=save_as,
                bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)


def handle_program_options():
    parser = argparse.ArgumentParser(description="Create a series of Principal\
                                     Coordinate plots for each OTU in an \
                                     input list where the plot points are \
                                     varied in size by the relative abundance \
                                     of the OTU relative to either Sample or\
                                     the total contribution of the OTU to the \
                                     data set.")
    parser.add_argument("-i", "--otu_table", required=True,
                        help="The biom-format file with OTU-Sample abundance \
                              data.")
    parser.add_argument("-m", "--mapping", required=True,
                        help="The mapping file specifying group information \
                              for each sample.")
    parser.add_argument("-pc", "--pcoa_fp", required=True,
                        help="Principal Coordinates Analysis file. \
                              Eg. unweighted_unifrac_pc.txt, or any other\
                              output from principal_coordinates.py.")
    parser.add_argument("-b", "--group_by", required=True,
                        help="Column name in mapping file specifying group\
                              information.")
    parser.add_argument("-c", "--colors", default=None,
                        help="A column name in the mapping file containing\
                              hexadecimal (#FF0000) color values that will\
                              be used to color the groups. Each sample ID must\
                              have a color entry.")
    parser.add_argument("-ids", "--otu_ids_fp", required=True,
                        help="Path to a file containing one OTU ID per line.\
                              One plot will be created for each OTU.")
    parser.add_argument("-o", "--output_dir", default=".",
                        help="The directory to output the PCoA plots to.")
    parser.add_argument("-s", "--save_as", default="svg",
                        help="The type of image file for PCoA plots. By\
                              default, files will be saved in SVG format.")
    parser.add_argument("--scale_by", default=1000, type=float,
                        help="Species relative abundance is multiplied by this \
                              factor in order to make appropriate visible \
                              bubbles in the output plots. Default is 10000.")
    parser.add_argument("--ggplot2_style", action="store_true",
                        help="Apply ggplot2 styling to the figure.")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Displays species name as each is being plotted \
                              and stored to disk.")

    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.otu_table):
            pass
    except IOError as ioe:
        sys.exit("\nError with BIOM format file:{}\n".format(ioe))

    try:
        with open(args.pcoa_fp):
            pass
    except IOError as ioe:
        sys.exit("\nError with principal coordinates file:{}\n".format(ioe))

    try:
        with open(args.mapping):
            pass
    except IOError as ioe:
        sys.exit("\nError with mapping file:{}\n".format(ioe))

    if not os.path.exists(args.output_dir):
        try:
            os.mkdir(args.output_dir)
        except OSError as oe:
            if os.errno == 2:
                msg = ("One or more directories in the path provided for " +
                       "--output-dir ({}) do not exist. If you are specifying " +
                       "a new directory for output, please ensure all other " +
                       "directories in the path currently exist.")
                sys.exit(msg.format(args.output_dir))
            else:
                msg = ("An error occurred trying to create the output " +
                       "directory ({}) with message: {}")
                sys.exit(msg.format(args.output_dir, oe.strerror))

    # load the BIOM table
    biomtbl = biom.load_table(args.otu_table)

    # Read unifrac principal coordinates file
    unifrac = util.parse_unifrac(args.pcoa_fp)

    # Read otu data file
    otus = set()
    with open(args.otu_ids_fp, "rU") as nciF:
        for line in nciF.readlines():
            line = line.strip()
            otus.add(line)

    # Gather categories from mapping file
    header, imap = util.parse_map_file(args.mapping)
    try:
        category_idx = header.index(args.group_by)
    except ValueError:
        msg = "Error: Specified mapping category '{}' not found."
        sys.exit(msg.format(args.group_by))
    category_ids = util.gather_categories(imap, header, [args.group_by])
    color_map = util.color_mapping(imap, header, args.group_by, args.colors)
    rel_abd = bc.relative_abundance(biomtbl)
    rel_abd = bc.arcsine_sqrt_transform(rel_abd)

    # plot samples based on relative abundance of some OTU ID
    for otuid in otus:
        otuname = oc.otu_name(biomtbl.metadata(otuid, axis="observation")["taxonomy"])
        cat_data = {cat: {"pc1": [], "pc2": [], "size": []}
                    for cat in category_ids}

        for sid in unifrac["pcd"]:
            category = cat_data[imap[sid][category_idx]]
            try:
                size = rel_abd[sid][otuid] * args.scale_by
            except KeyError as ke:
                print("{} not found in {} sample.".format(ke, sid))
                continue
            category["pc1"].append(float(unifrac["pcd"][sid][0]))
            category["pc2"].append(float(unifrac["pcd"][sid][1]))
            category["size"].append(size)

        if args.verbose:
            print("Saving chart for {}".format(" ".join(otuname.split("_"))))
        xr, yr = calculate_xy_range(cat_data)
        plot_PCoA(cat_data, otuname, unifrac, color_map.keys(),
                  color_map, xr, yr, args.output_dir,
                  args.save_as, args.ggplot2_style)

if __name__ == "__main__":
    sys.exit(main())
