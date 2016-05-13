#!/usr/bin/env python
import sys
import argparse
from warnings import catch_warnings as cw
from phylotoast import util, biom_calc as bc, otu_calc as oc
errors = []
try:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
except ImportError as ie:
    errors.append(ie)
try:
    import numpy as np
except ImportError as ie:
    errors.append(ie)
try:
    import pandas as pd
except ImportError as ie:
    errors.append(ie)
try:
    import biom
except ImportError as ie:
    errors.append(ie)
try:
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
except ImportError as ie:
    errors.append(ie)
if len(errors) != 0:
    for item in errors:
        print "Import Error:", item
    sys.exit()


def get_relative_abundance(biomfile):
    """
    Return relative abundance from a OTU table. OTUIDs are converted to their
    genus-species identifier.
    """
    biomf = biom.load_table(biomfile)
    norm_biomf = biomf.norm(inplace=False)
    rel_abd = {}
    for sid in norm_biomf.ids():
        rel_abd[sid] = {}
        for otuid in norm_biomf.ids("observation"):
            otuname = oc.otu_name(norm_biomf.metadata(otuid, axis="observation")["taxonomy"])
            abd = norm_biomf.get_value_by_ids(otuid, sid)
            rel_abd[sid][otuname] = abd
    ast_rel_abd = bc.arcsine_sqrt_transform(rel_abd)
    return ast_rel_abd


def plot_LDA(X_lda, y_lda, class_colors, exp_var, out_fp=""):
    """
    Plot transformed LDA data.
    """
    cats = class_colors.keys()
    group_lda = {c: [] for c in cats}
    fig = plt.figure(figsize=(15, 10))
    for i, target_name in zip(range(len(cats)), cats):
        cat_x = X_lda[:, 0][y_lda == target_name]
        if X_lda.shape[1] == 1:
            cat_y = np.ones((cat_x.shape[0], 1)) + i
        else:
            cat_y = X_lda[:, 1][y_lda == target_name]
        group_lda[target_name].append(cat_x)
        group_lda[target_name].append(cat_y)
        plt.scatter(x=cat_x, y=cat_y, label=target_name,
                    color=class_colors[target_name],
                    alpha=0.85, s=250, edgecolors="k")
    mpl.rc("font", family="Arial")  # define font for figure text
    mpl.rc('xtick', labelsize=12)  # increase X axis ticksize
    mpl.rc('ytick', labelsize=12)  # increase Y axis ticksize
    if X_lda.shape[1] == 1:
        plt.ylim((0.5, 2.5))
    plt.xlabel("LD1 (Percent Explained Variance: {:.3f}%)".format(exp_var[0]*100), fontsize=16)
    plt.ylabel("LD2 (Percent Explained Variance: {:.3f}%)".format(exp_var[1]*100), fontsize=16)
    leg = plt.legend(loc="best", frameon=True, framealpha=1, fontsize=16)
    leg.get_frame().set_edgecolor('k')

    # save or display result
    if out_fp:
        plt.savefig(out_fp, facecolor="white", edgecolor="none", dpi=300,
                    bbox_inches="tight", pad_inches=0.1)
    else:
        plt.show()


def run_LDA(df):
    """
    Run LinearDiscriminantAnalysis on input dataframe (df) and return
    transformed data, scalings and
    """
    # Prep variables for sklearn LDA
    X = df[range(1, df.shape[1])].values     # input data matrix
    y = df["Condition"].values               # data categories list

    # Calculate LDA
    sklearn_lda = LDA()
    X_lda_sklearn = sklearn_lda.fit_transform(X, y)
    exp_var = sklearn_lda.explained_variance_ratio_

    return X_lda_sklearn, y, exp_var


def handle_program_options():
    parser = argparse.ArgumentParser(description="Create an LDA plot from\
                                     sample-grouped OTU data. It is necessary\
                                     to remove the header cell '#OTU ID'\
                                     before running this program.")
    parser.add_argument("-i", "--input_data_type", required=True,
                        choices=["biom", "unifrac_dm"],
                        default="unifrac_dm",
                        help="Specify if the input file is biom file format OTU \
                              table or unifrac distance matrix. If biom file is \
                              provided, the arc-sine transformed relative abundances \
                              eill be used as input whereas, if unifrac distance matrix.\
                              is given, unifrac distances will be used as input to LDA.\
                              [REQUIRED]")
    parser.add_argument("-d", "--input_data", required=True,
                        help="Input biom file format or unifrac distance matrix \
                              (output of beta_diversity.py from QIIME). [REQUIRED]")
    parser.add_argument("-m", "--map_fp", required=True,
                        help="Metadata mapping file. [REQUIRED]")
    parser.add_argument("-g", "--group_by", required=True,
                        help="A column name in the mapping file containing\
                              categorical values that will be used to identify \
                              groups. Each sample ID must have a group entry. \
                              Default is no categories and all the data will be \
                              treated as a single group.")
    parser.add_argument("-c", "--color_by", required=True,
                        help="A column name in the mapping file containing\
                              hexadecimal (#FF0000) color values that will\
                              be used to color the groups. Each sample ID must\
                              have a color entry.")
    parser.add_argument("--bubble", action="store_true",
                        help="If set, provide a file with 1 OTU name per line \
                              for bubble plotting. OTU name must be condensed to \
                              genus-species identifier. Default parameter value \
                              will not plot bubble plots.")
    parser.add_argument("--save_lda_input",
                        help="Save a CSV-format file of the transposed LDA-input\
                              table to the file specifed by this option.")
    parser.add_argument("--plot_title", default=None,
                        help="Plot title. Default is no title.")
    parser.add_argument("-o", "--out_fp", default="",
                        help="The path and file name to save the plot under.\
                              If specified, the figure will be saved directly\
                              instead of opening a window in which the plot \
                              can be viewed before saving")
    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.map_fp):
            pass
    except IOError as ioe:
        err_msg = "\nError in metadata mapping filepath (-m): {}\n"
        sys.exit(err_msg.format(ioe))

    if args.bubble:
        try:
            with open(args.bubble) as hojiehr:
                bubble_otus = [line.strip() for line in hojiehr.readlines()]
        except IOError as ioe:
            err_msg = "\nError in OTU name list file (--bubble): {}\n"
            sys.exit(err_msg.format(ioe))

    # Parse and read mapping file
    header, imap = util.parse_map_file(args.map_fp)
    class_colors = util.color_mapping(imap, header, args.group_by, args.color_by)

    if args.input_data_type == "unifrac_dm":
        try:
            with open(args.input_data):
                pass
        except IOError as ioe:
            err_msg = "\nError with unifrac distance matrix file (-d): {}\n"
            sys.exit(err_msg.format(ioe))
        uf_data = pd.read_csv(args.input_data, sep="\t", index_col=0)
        uf_data.insert(0, "Condition", [imap[sid][header.index(args.group_by)]
                                        for sid in uf_data.index])
        if args.save_lda_input:
            uf_data.to_csv(args.save_lda_input, sep="\t")
        # Run LDA
        X_lda, y_lda, exp_var = run_LDA(uf_data)
        # Plot LDA Bubble
        plot_LDA(X_lda, y_lda, class_colors, exp_var, out_fp=args.out_fp)
    else:
        try:
            rel_abd = get_relative_abundance(args.input_data)  # load biom file
        except ValueError as ve:
            err_msg = "\nError with biom format file (-d): {}\n"
            sys.exit(err_msg.format(ve))
        df_rel_abd = pd.DataFrame(rel_abd).T
        df_rel_abd.insert(0, "Condition", [imap[sid][header.index(args.group_by)]
                                           for sid in df_rel_abd.index])
        if args.save_lda_input:
            df_rel_abd.to_csv(args.save_lda_input, sep="\t")
        # Run LDA
        X_lda, y_lda, exp_var = run_LDA(df_rel_abd)
        # Plot LDA
        plot_LDA(X_lda, y_lda, class_colors, exp_var, out_fp=args.out_fp)


if __name__ == "__main__":
    sys.exit(main())
