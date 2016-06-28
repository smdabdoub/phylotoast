#!/usr/bin/env python
"""
Calculate and plot alpha diversity for two or more sample categories.
"""
import sys
import csv
import argparse
import os.path as osp
from collections import defaultdict
from phylotoast import graph_util as gu, util as putil
importerrors = []
try:
    import biom
except ImportError as ie:
    importerrors.append(ie)
try:
    import scipy.stats as stats
except ImportError as ie:
    importerrors.append(ie)
try:
    from skbio.diversity import alpha
except ImportError as ie:
    importerrors.append(ie)
try:
    from matplotlib import pyplot as plt, gridspec
except ImportError as ie:
    importerrors.append(ie)
if len(importerrors) != 0:
    for item in importerrors:
        print "Import Error. Please install missing module:", item
    sys.exit()


def gather_samples(biomT):
    return {sid: biomT.data(sid).astype(int) for sid in biomT.ids()}


def calc_diversity(method, parsed_mapf, biom, cats, cats_index):
    counts = {cat: [] for cat in cats}
    sample_ids = []

    for sid, sample_counts in gather_samples(biom).items():
        sample_ids.append(sid)
        if sid in parsed_mapf:
            counts[parsed_mapf[sid][cats_index]].append((sid, sample_counts))

    div_calc = {cat: {count[0]: method(count[1]) for count in counts}
                for cat, counts in counts.items()}

    return div_calc, sample_ids


def print_MannWhitneyU(div_calc):
    """
    Compute the Mann-Whitney U test for unequal group sample sizes.
    """
    try:
        x = div_calc.values()[0].values()
        y = div_calc.values()[1].values()
    except:
        return "Error setting up input arrays for Mann-Whitney U Test. Skipping "\
               "significance testing."
    T, p = stats.mannwhitneyu(x, y)
    print "\nMann-Whitney U test statistic:", T
    print "Two-tailed p-value: {}".format(2 * p)


def print_KruskalWallisH(div_calc):
    """
    Compute the Kruskal-Wallis H-test for independent samples. A typical rule is that
    each group must have at least 5 measurements.
    """
    calc = defaultdict(list)
    try:
        for k1, v1 in div_calc.iteritems():
            for k2, v2 in v1.iteritems():
                calc[k1].append(v2)
        print calc
    except:
        return "Error setting up input arrays for Kruskal-Wallis H-Test. Skipping "\
               "significance testing."
    h, p = stats.kruskal(*calc.values())
    print "\nKruskal-Wallis H-test statistic for {} groups: {}".format(str(len(div_calc)), h)
    print "p-value: {}".format(p)


def plot_group_diversity(diversities, grp_colors, title, diversity_type, out_dir, plot_ext):
    fig_div = plt.figure(figsize=(21, 7))
    grid = gridspec.GridSpec(1, 2)

    # Disease States Shannon Diversity plots
    ax_div = fig_div.add_subplot(grid[0, 0])

    for i, grp in enumerate(diversities):
        gu.plot_kde(diversities[grp].values(), ax_div, title, grp_colors[grp])

    ax_div.set_xlabel(diversity_type.title())
    ax_div.set_ylabel("Density")
    ax_div.legend([plt.Rectangle((0, 0), 1, 1, fc=color) for color in grp_colors.values()],
                  grp_colors.keys(), loc="best")

    fig_div.savefig(osp.join(out_dir, diversity_type+"."+plot_ext), facecolor="white",
                    edgecolor="none", bbox_inches="tight", pad_inches=0.2)


def write_diversity_metrics(data, sample_ids, fp=None):
    """
    Given a dictionary of diversity calculations (keyed by method)
    write out the data to a file.
    """
    if fp is None:
        fp = "./diversity_data.txt"

    with open(fp, "w") as outf:
        out = csv.writer(outf, delimiter="\t")
        out.writerow(["SampleID", "Group", "Calculation"])
        for group, d in data.iteritems():
            for sid, value in d.iteritems():
                out.writerow([sid, group, value])


def handle_program_options():
    """Parses the given options passed in at the command line."""
    parser = argparse.ArgumentParser(description="Calculate the alpha diversity\
                                     of a set of samples using one or more \
                                     metrics and output a kernal density \
                                     estimator-smoothed histogram of the \
                                     results.")
    parser.add_argument("-m", "--map_file",
                        help="QIIME mapping file.")
    parser.add_argument("-i", "--biom_fp",
                        help="Path to the BIOM table")
    parser.add_argument("-c", "--category",
                        help="Specific category from the mapping file.")
    parser.add_argument("-d", "--diversity", default=["shannon"], nargs="+",
                        help="The alpha diversity metric. Default \
                             value is 'shannon', which will calculate the Shannon\
                             entropy. Multiple metrics can be specified (space separated).\
                             The full list of metrics is available at:\
                             http://scikit-bio.org/docs/latest/generated/skbio.diversity.alpha.html.\
                             Beta diversity metrics will be supported in the future.")
    parser.add_argument("--color_by",
                        help="A column name in the mapping file containing\
                              hexadecimal (#FF0000) color values that will\
                              be used to color the groups. Each sample ID must\
                              have a color entry.")
    parser.add_argument("--plot_title", default="",
                        help="A descriptive title that will appear at the top \
                        of the output plot. Surround with quotes if there are\
                        spaces in the title.")
    parser.add_argument("-o", "--out_dir", default=".",
                        help="The directory plots will be saved to.")
    parser.add_argument("--image_type", default="png",
                        help="The type of image to save: png, svg, pdf, eps, etc...")
    parser.add_argument("--save_calculations",
                        help="Path and name of text file to store the calculated "
                        "diversity metrics.")
    parser.add_argument("--show_significance", action="store_false", help="Display "
                        "significance testing results. The results will be shown by "
                        "default.")
    parser.add_argument("--show_available_metrics", action="store_true",
                        help="Supply this parameter to see which alpha diversity metrics "
                             " are available for usage. No calculations will be performed"
                             " if this parameter is provided.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    if args.show_available_metrics:
        print "\nAvailable alpha diversity metrics:"
        return "\n".join(alpha.__all__)
    else:
        return "ERROR! Cannot print out alpha diversity metrics list."

    try:
        with open(args.map_file):
            pass
    except IOError as ioe:
            err_msg = "\nError opening QIIME mapping file: {}\n"
            sys.exit(err_msg.format(ioe))

    try:
        with open(args.biom_fp):
            pass
    except IOError as ioe:
            err_msg = "\nError opening BIOM table file: {}\n"
            sys.exit(err_msg.format(ioe))

    header, sample_map = putil.parse_map_file(args.map_file)
    biom_tbl = biom.load_table(args.biom_fp)
    if args.category not in header:
        sys.exit("Category '{}' not found".format(args.category))

    cat_idx = header.index(args.category)
    cat_vals = {entry[cat_idx] for entry in sample_map.values()}

    plot_title = args.plot_title

    colors = putil.color_mapping(sample_map, header, args.category, args.color_by)

    # Perform diversity calculations and density plotting
    for method in args.diversity:
        if method not in alpha.__all__:
            sys.exit("ERROR: Diversity metric not found: " + method)
        metric = eval("alpha."+method)
        div_calc, sample_ids = calc_diversity(metric, sample_map, biom_tbl,
                                              cat_vals, cat_idx)

        plot_group_diversity(div_calc, colors, plot_title, method,
                             args.out_dir, args.image_type)

        # calculate and print significance testing results
        if args.show_significance:
            print "Diversity significance testing: {}".format(method.title())
            if len(cat_vals) == 2:
                print_MannWhitneyU(div_calc)
            elif len(cat_vals) > 2:
                print_KruskalWallisH(div_calc)
            print
        else:
            continue

        if args.save_calculations:
            write_diversity_metrics(div_calc, sample_ids, args.save_calculations)


if __name__ == "__main__":
    sys.exit(main())
