#!/usr/bin/env python
"""
Abstract: This script returns LDA plots, with samples/dots sized by relative abundances
          of input OTU(s).
Date: 07/18/2016
"""

import sys
import argparse
from os.path import join as pj
from phylotoast import biom_calc as bc, otu_calc as oc, graph_util as gu, util
errors = []
try:
    import biom
except ImportError as ie:
    errors.append(ie)
try:
    import pandas as pd
except ImportError as ie:
    errors.append(ie)
try:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rc("font", family="Arial")  # define font for figure text
    mpl.rc("xtick", labelsize=12)  # increase X axis ticksize
    mpl.rc("ytick", labelsize=12)  # increase Y axis ticksize
except ImportError as ie:
    errors.append(ie)
try:
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
except ImportError as ie:
    errors.append(ie)
if len(errors) != 0:
    for item in errors:
        print("Import Error:", item)
    sys.exit()


def run_LDA(df):
    """
    Run LinearDiscriminantAnalysis on input dataframe (df) and return
    transformed data, scalings and explained variance by discriminants.
    """
    # Prep variables for sklearn LDA
    X = df[range(1, df.shape[1])].values     # input data matrix
    y = df["Condition"].values               # data categories list

    # Calculate LDA
    sklearn_lda = LDA()
    X_lda_sklearn = sklearn_lda.fit_transform(X, y)
    try:
        exp_var = sklearn_lda.explained_variance_ratio_
    except AttributeError as ae:
        print("\n{}: explained variance cannot be computed.\nPlease check this GitHub PR:"
              " https://github.com/scikit-learn/scikit-learn/pull/6027".format(ae))
        return X_lda_sklearn, y, "NA"

    return X_lda_sklearn, y, exp_var


def handle_program_options():
    """Command line arguments."""
    parser = argparse.ArgumentParser(description="Create an LDA bubble plot from  either \
                                                 sample-grouped OTU abundance data or \
                                                 sample-wise distance matrix file.")
    parser.add_argument("-i", "--otu_table", required=True,
                        help="Input biom file format OTU table. [REQUIRED]")
    parser.add_argument("-m", "--map_fp", required=True,
                        help="Metadata mapping file. [REQUIRED]")
    parser.add_argument("-g", "--group_by", required=True,
                        help="A column name in the mapping file containing categorical \
                              values that will be used to identify groups. Each sample \
                              ID must have a group entry. Default is no categories and \
                              all the data will be treated as a single group. [REQUIRED]")
    parser.add_argument("-c", "--color_by", required=True,
                        help="A column name in the mapping file containing hexadecimal \
                              (#FF0000) color values that will be used to color the \
                              groups. Each sample ID must have a color entry. [REQUIRED]")
    parser.add_argument("-ids", "--otu_ids_fp", required=True,
                        help="Path to a file containing one OTU ID per line. One plot \
                              will be created for each OTU. [REQUIRED]")
    parser.add_argument("-dm", "--dist_matrix_file",
                        help="Input distance matrix file.")
    parser.add_argument("--save_lda_input",
                        help="Save a CSV-format file of the transposed LDA-input table \
                              to the file specifed by this option.")
    parser.add_argument("-od", "--output_dir", default=".",
                        help="The directory to save the LDA bubble plots to. By default, \
                              plots will be saved in current working directory.")
    parser.add_argument("--scale_by", default=1000, type=float,
                        help="Species relative abundance is multiplied by this factor in \
                              order to make appropriate visible bubbles in the output \
                              plots. Default scaling is 1000.")
    parser.add_argument("-s", "--save_as", default="svg",
                        help="The type of image file for LDA plots. By default, plots \
                              will be saved in 'svg' format.")
    parser.add_argument("--ggplot2_style", action="store_true",
                        help="Apply ggplot2 styling to the figure.")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Displays species name as each is being plotted and stored \
                              to disk.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    # Parse and read mapping file
    try:
        header, imap = util.parse_map_file(args.map_fp)
        category_idx = header.index(args.group_by)
    except IOError as ioe:
        err_msg = "\nError in metadata mapping filepath (-m): {}\n"
        sys.exit(err_msg.format(ioe))
    # Obtain group colors
    class_colors = util.color_mapping(imap, header, args.group_by, args.color_by)

    # Get otus for LDA bubble plots
    try:
        with open(args.otu_ids_fp) as oif:
            for line in oif.readlines():
                bubble_otus = line.strip().split("\r")
    except IOError as ioe:
        err_msg = "\nError in OTU IDs file (--bubble): {}\n"
        sys.exit(err_msg.format(ioe))

    # Load biom file and calculate relative abundance
    try:
        biomf = biom.load_table(args.otu_table)
    except IOError as ioe:
        err_msg = "\nError with biom format file (-d): {}\n"
        sys.exit(err_msg.format(ioe))

    # Get normalized relative abundances
    rel_abd = bc.relative_abundance(biomf)
    rel_abd = bc.arcsine_sqrt_transform(rel_abd)

    # Set up input for LDA calc and get LDA transformed data
    if args.dist_matrix_file:
        try:
            uf_data = pd.read_csv(args.dist_matrix_file, sep="\t", index_col=0)
        except IOError as ioe:
            err_msg = "\nError with unifrac distance matrix file (-d): {}\n"
            sys.exit(err_msg.format(ioe))
        uf_data.insert(0, "Condition", [imap[sid][category_idx] for sid in uf_data.index])
        sampleids = uf_data.index
        if args.save_lda_input:
            uf_data.to_csv(args.save_lda_input, sep="\t")
        # Run LDA
        X_lda, y_lda, exp_var = run_LDA(uf_data)
    else:
        df_rel_abd = pd.DataFrame(rel_abd).T
        df_rel_abd.insert(0, "Condition", [imap[sid][category_idx] for sid in df_rel_abd.index])
        sampleids = df_rel_abd.index
        if args.save_lda_input:
            df_rel_abd.to_csv(args.save_lda_input, sep="\t")
        # Run LDA
        X_lda, y_lda, exp_var = run_LDA(df_rel_abd)

    # Calculate position and size of SampleIDs to plot for each OTU
    for otuid in bubble_otus:
        otuname = oc.otu_name(biomf.metadata(otuid, axis="observation")["taxonomy"])
        plot_data = {cat: {"x": [], "y": [], "size": [], "label": []}
                     for cat in class_colors.keys()}
        for sid, data in zip(sampleids, X_lda):
            category = plot_data[imap[sid][category_idx]]
            try:
                size = rel_abd[sid][otuid] * args.scale_by
            except KeyError as ke:
                print("{} not found in {} sample.".format(ke, sid))
                continue
            category["x"].append(float(data[0]))
            category["y"].append(float(data[1]))
            category["size"].append(size)

        # Plot LDA bubble for each OTU
        fig = plt.figure(figsize=(14, 8))
        ax = fig.add_subplot(111)
        for i, cat in enumerate(plot_data):
            plt.scatter(plot_data[cat]["x"], plot_data[cat]["y"],
                        plot_data[cat]["size"], label=cat, color=class_colors[cat],
                        alpha=0.85, marker="o", edgecolor="k")
        if X_lda.shape[1] == 1:
            plt.ylim((0.5, 2.5))
        plt.title(" ".join(otuname.split("_")), style="italic", fontsize=13)
        try:
            plt.xlabel("LD1 (Percent Explained Variance: {:.3f}%)".format(exp_var[0]*100),
                       fontsize=13, labelpad=15)
        except:
            plt.xlabel("LD1", fontsize=13, labelpad=15)
        try:
            plt.ylabel("LD2 (Percent Explained Variance: {:.3f}%)".format(exp_var[1]*100),
                       fontsize=13, labelpad=15)
        except:
            plt.ylabel("LD2", fontsize=13, labelpad=15)

        lgnd = plt.legend(loc="best", scatterpoints=3, fontsize=13)
        for i in range(len(class_colors.keys())):
            lgnd.legendHandles[i]._sizes = [80]  # Change the legend marker size manually

        # Set style for LDA bubble plots
        if args.ggplot2_style:
            gu.ggplot2_style(ax)
            fc = "0.8"
        else:
            fc = "none"

        # Save LDA bubble plots to output directory
        if args.verbose:
            print("Saving chart for {}".format(" ".join(otuname.split("_"))))
        fig.savefig(pj(args.output_dir, "_".join(otuname.split())) + "." + args.save_as,
                    facecolor=fc, edgecolor="none", dpi=300,
                    bbox_inches="tight", pad_inches=0.2)
        plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
