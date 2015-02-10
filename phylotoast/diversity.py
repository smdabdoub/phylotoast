#!/usr/bin/env python
"""
Calculate and plot for two sample categories: Shannon diversity,
Chao1 diversity, and a Jaccard similiarity distance matrix heatmap
"""
import sys
import argparse
import itertools
import os.path as osp
importerrors = []
try:
    import biom
except ImportError as ie:
    importerrors.append('biom-format')
try:
    from skbio.diversity import alpha, beta
except ImportError as ie:
    importerrors.append('scikit-bio')
try:
    from brewer2mpl import qualitative
except ImportError as ie:
    importerrors.append('brewer2mpl')
if len(importerrors) != 0:
    for item in importerrors:
        print 'Import Error. Please install missing module:', item
    sys.exit()
from matplotlib import pyplot as plt, gridspec
from phylotoast import graph_util as gu, util as putil


def gather_samples(biomT):
    return {sid: biomT.data(sid).astype(int) for sid in biomT.ids()}


def calc_diversities(parsed_mapf, biom, cats, cats_index):
    counts = {cat: [] for cat in cats}

    for sid, sample_counts in gather_samples(biom).items():
        counts[parsed_mapf[sid][cats_index]].append(sample_counts)

    sdiv = {cat: [alpha.shannon(count)
            for count in counts] for cat, counts in counts.items()}
    chao1 = {cat: [alpha.chao1(count) for count in counts]
             for cat, counts in counts.items()}

    return sdiv, chao1


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

    fig_div.savefig(osp.join(out_dir, title+'.'+plot_ext), facecolor='white',
                    edgecolor='none', bbox_inches='tight', pad_inches=0.5)


def handle_program_options():
    """Parses the given options passed in at the command line."""
    parser = argparse.ArgumentParser(description="Creates a KEGG map\
             visualization of user-specified pathways with path lines\
             optionally scaled by supplied abundance data.")
    parser.add_argument('map_file',
                        help="QIIME mapping file.")
    parser.add_argument('biom_file',
                        help="BIOM table file name")
    parser.add_argument('category',
                        help="Specific category from the mapping file.")
    parser.add_argument('-c', default=None, nargs='+',
                        choices=['shannon', 'chao1', 'jaccard'],
                        help='Choose the type of calculation needed. Default'
'                             value is "None", which will output all 3'
'                             types of calculations.')
    parser.add_argument('plot_title',
                        help="The name of a PDF file the pathway map will be\
                              written to.")
    parser.add_argument('out_dir',
                        help="The directory all plots will be saved to.")
    parser.add_argument('-p', '--image_type', default='png',
                        help="The type of image to save: PNG, SVG, etc...")
    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.map_file):
            pass
    except IOError as ioe:
            err_msg = '\nError in QIIME mapping file: {}\n'
            sys.exit(err_msg.format(ioe))

    try:
        with open(args.biom_file):
            pass
    except IOError as ioe:
            err_msg = '\nError in BIOM table file: {}\n'
            sys.exit(err_msg.format(ioe))

    header, sample_map = putil.parse_map_file(args.map_file)
    biom_tbl = biom.load_table(args.biom_file)
    if args.category not in header:
        raise(ValueError('Category \'{}\' not found'.format(args.category)))
    cat_idx = header.index(args.category)
    cat_vals = {entry[cat_idx] for entry in sample_map.values()}

    sdiv, chao1 = calc_diversities(sample_map, biom_tbl, cat_vals, cat_idx)

    plot_title = args.plot_title

    bmap = qualitative.Paired[12]
    bcolors = itertools.cycle(bmap.hex_colors)
    colors = {}
    for cat in cat_vals:
        colors[cat] = bcolors.next()

    # Shannon Diversity
    if args.c is None or 'shannon' in args.c:
        plot_group_diversity(sdiv, colors, 'Shannon Diversity - '+plot_title,
                             'Shannon Diversity', args.out_dir,
                             args.image_type)

    # Chao1 Diversity
    if args.c is None or 'chao1' in args.c:
        plot_group_diversity(chao1, colors, 'Chao1 Diversity - '+plot_title,
                             'Chao1 Diversity', args.out_dir, args.image_type)

    # Jaccard Similarity
    if args.c is None or 'jaccard' in args.c:
        jaccard = beta.pw_distances_from_table(biom.load_table(args.biom_file),
                                               metric='jaccard')
        fig_jaccard = jaccard.plot(title='Jaccard similarity heatmap - '+plot_title)
        fig_jaccard.savefig(osp.join(args.out_dir, 'jaccard.'+args.image_type),
                            facecolor='white', edgecolor='none')


if __name__ == '__main__':
    main()
