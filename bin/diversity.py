#!/usr/bin/env python
"""
Calculate and plot for two sample categories: Shannon diversity,
Chao1 diversity, and a Jaccard similiarity distance matrix heatmap
"""
# this is only to suppress the pre 1.5.0 matplotlib numpy-related warning - DELETE
import warnings
warnings.filterwarnings("ignore")

import sys
import argparse
import csv
import itertools
import os.path as osp
importerrors = []
try:
    import biom
except ImportError as ie1:
    importerrors.append(ie1)
try:
    import scipy.stats as stats
except ImportError as ie2:
    importerrors.append(ie2)
try:
    from skbio.diversity import alpha, beta
except ImportError as ie3:
    importerrors.append(ie3)
try:
    import matplotlib
except ImportError as ie4:
    importerrors.append(ie4)
try:
    from brewer2mpl import qualitative
except ImportError as ie5:
    importerrors.append(ie5)
if len(importerrors) != 0:
    for item in importerrors:
        print 'Import Error:', item
    sys.exit()


#matplotlib.use("Agg")  # for use on headless server
from matplotlib import pyplot as plt, gridspec
from phylotoast import graph_util as gu, util as putil


def color_mapping(sample_map, header, group_column, color_column=None):
    """
    Determine color-category mapping. If color_column was specified, then
    map the category names to color values. Otherwise, use the brewer colors
    to automatically generate a set of colors for the group values.
    """
    group_colors = {}
    group_gather = putil.gather_categories(sample_map, header, [group_column])

    if color_column is not None:
        color_gather = putil.gather_categories(sample_map, header, [color_column])
        # match sample IDs between color_gather and group_gather
        for group in group_gather:
            for color in color_gather:
                # allow incomplete assignment of colors, if group sids overlap at
                # all with the color sids, consider it a match
                if group_gather[group].sids.intersection(color_gather[color].sids):
                    group_colors[group] = color
    else:
        bmap = qualitative.Paired[12]
        bcolors = itertools.cycle(bmap.hex_colors)
        for group in group_gather:
            group_colors[group] = bcolors.next()

    return group_colors

def gather_samples(biomT):
    return {sid: biomT.data(sid).astype(int) for sid in biomT.ids()}


def calc_diversity(method, parsed_mapf, biom, cats, cats_index):
    counts = {cat: [] for cat in cats}
    sample_ids = []

    for sid, sample_counts in gather_samples(biom).items():
        sample_ids.append(sid)
        if sid in parsed_mapf:
            counts[parsed_mapf[sid][cats_index]].append(sample_counts)

    div_calc = {cat: [method(count) for count in counts] for cat, counts in counts.items()}

    return div_calc, sample_ids


def print_WilcoxonSRT(x, y=None):
    """
    Compute the Wilcoxon Signed Rank Test
    """
    T,p = stats.wilcoxon(x, y)
    print "Wilcoxon Signed Rank Test:"
    print "p-value: {}".format(p)


def print_KruskalWallisH(div_calc):
    """
    Compute the Kruskal-Wallis H-test for independent samples
    """
    h, p = stats.kruskal(*div_calc)
    print "Kruskal-Wallis H-test for {} groups:".format(str(len(div_calc)))
    print "p-value: {}".format(p)


def plot_group_diversity(diversities, grp_colors, title, diversity_type, out_dir, plot_ext):
    fig_div = plt.figure(figsize=(21, 7))
    grid = gridspec.GridSpec(1, 2)

    # Disease States Shannon Diversity plots
    ax_div = fig_div.add_subplot(grid[0, 0])

    for i, grp in enumerate(diversities):
        gu.plot_kde(diversities[grp], ax_div, title, grp_colors[grp])

    ax_div.set_xlabel(diversity_type)
    ax_div.set_ylabel('Density')
    ax_div.legend([plt.Rectangle((0, 0), 1, 1, fc=color) for color in grp_colors.values()],
                  grp_colors.keys(), loc='best')

    fig_div.savefig(osp.join(out_dir, diversity_type+'.'+plot_ext), facecolor='white',
                    edgecolor='none', bbox_inches='tight', pad_inches=0.5)

def write_diversity_metrics(data, sample_ids, fp=None):
    """
    Given a dictionary of diversity calculations (keyed by method)
    write out the data to a file.
    """
    if fp is None:
        fp = "./diversity_data.txt"

    with open(fp, 'w') as outf:
        out = csv.writer(outf, delimiter='\t')
        out.writerow(['calculation', 'group'])
        for group in data:
            for entry in data[group]:
                out.writerow([entry, group])


def handle_program_options():
    """Parses the given options passed in at the command line."""
    parser = argparse.ArgumentParser(description="Calculate the alpha diversity\
                                     of a set of samples using one or more \
                                     metrics and output a kernal density \
                                     estimator-smoothed histogram of the \
                                     results.")
    parser.add_argument('-m', '--map_file',
                        help="QIIME mapping file.")
    parser.add_argument('-i', '--biom_fp',
                        help="Path to the BIOM table")
    parser.add_argument('-c', '--category',
                        help="Specific category from the mapping file.")
    parser.add_argument('-d', '--diversity', default=['shannon'], nargs='+',
                        help='The alpha diversity metric. Default \
                             value is "shannon", which will calculate the Shannon\
                             entropy. Multiple metrics can be specified (space separated).\
                             The full list of metrics is available at:\
                             http://scikit-bio.org/docs/latest/generated/skbio.diversity.alpha.html.\
                             Beta diversity metrics will be supported in the future.')
    parser.add_argument('--color_by',
                        help="A column name in the mapping file containing\
                              hexadecimal (#FF0000) color values that will\
                              be used to color the groups. Each sample ID must\
                              have a color entry.")
    parser.add_argument('--plot_title', default='',
                        help="A descriptive title that will appear at the top \
                        of the output plot. Surround with quotes if there are\
                        spaces in the title.")
    parser.add_argument('--x_label', default='diversity', nargs='+',
                        help="The name of the diversity metric to be displayed on the\
                        plot as the X-axis label. If multiple metrics are specified,\
                        then multiple entries for the X-axis label should be given.")
    parser.add_argument('-o', '--out_dir', default=".",
                        help="The directory plots will be saved to.")
    parser.add_argument('--image_type', default='png',
                        help="The type of image to save: PNG, SVG, PDF, EPS, etc...")
    parser.add_argument('--save_calculations',
                        help="A text file containing the calculated metrics will be stored\
                        to the file system. Default: 'diversity_data.txt'")
    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.map_file):
            pass
    except IOError as ioe:
            err_msg = '\nError opening QIIME mapping file: {}\n'
            sys.exit(err_msg.format(ioe))

    try:
        with open(args.biom_fp):
            pass
    except IOError as ioe:
            err_msg = '\nError opening BIOM table file: {}\n'
            sys.exit(err_msg.format(ioe))

    header, sample_map = putil.parse_map_file(args.map_file)
    biom_tbl = biom.load_table(args.biom_fp)
    if args.category not in header:
        sys.exit('Category \'{}\' not found'.format(args.category))

    cat_idx = header.index(args.category)
    cat_vals = {entry[cat_idx] for entry in sample_map.values()}

    plot_title = args.plot_title

    colors = color_mapping(sample_map, header, args.category, args.color_by)

    # Perform diversity calculations and density plotting
    for method, x_label in zip(args.diversity, args.x_label):
        if method not in alpha.__all__:
            sys.exit("ERROR: Diversity metric not found: " + method)
        metric = eval('alpha.'+method)
        div_calc, sample_ids = calc_diversity(metric, sample_map, biom_tbl, cat_vals, cat_idx)

        plot_group_diversity(div_calc, colors, plot_title, x_label,
                             args.out_dir, args.image_type)

        print "Diversity significance testing: {}".format(x_label)
        # calculate and print significance testing results
        if len(cat_vals) == 2:
            print_WilcoxonSRT(*div_calc.values())
        elif len(cat_vals) > 2:
            print_KruskalWallisH(div_calc.values())
        print

        if args.save_calculations:
            prefix = '_'.join(x_label.split())
            write_diversity_metrics(div_calc, sample_ids, osp.join(args.out_dir, args.save_calculations))


    # # Chao1 Diversity
    # if args.diversity is None or 'chao1' in args.diversity:
    #     plot_group_diversity(chao1, colors, 'Chao1 Diversity - '+plot_title,
    #                          'Chao1 Diversity', args.out_dir, args.image_type)

    # # Jaccard Similarity
    # if args.diversity is None or 'jaccard' in args.diversity:
    #     jaccard = beta.pw_distances_from_table(biom.load_table(args.biom_fp),
    #                                            metric='jaccard')
    #     fig_jaccard = jaccard.plot(title='Jaccard similarity heatmap - '+plot_title)
    #     fig_jaccard.savefig(osp.join(args.out_dir, 'jaccard.'+args.image_type),
    #                         facecolor='white', edgecolor='none')


if __name__ == '__main__':
    main()
