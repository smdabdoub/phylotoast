#!/usr/bin/env python
# coding: utf-8
"""
Given a set of core microbiome files, create a matching set of ovelapping
barplots that visualize which species belong to each core microbiome.
"""
from __future__ import absolute_import, division, print_function

import ast
import argparse
from collections import Counter, OrderedDict, defaultdict
import csv
import itertools as it
import os.path as osp
import sys

from phylotoast import otu_calc as oc, util

importerrors = []
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.collections import PolyCollection
    from matplotlib.ticker import MaxNLocator

except ImportError as ie:
    importerrors.append(ie)
if len(importerrors) != 0:
    for item in importerrors:
        print("Import Error. Please install missing module:", item)
    sys.exit()

__author__ = "Shareef M. Dabdoub"
__copyright__ = "Copyright 2016, Shareef M. Dabdoub"
__credits__ = ["Shareef M. Dabdoub", "Akshay Paropkari", 
               "Sukirth Ganesan", "Purnima Kumar"]
__license__ = "MIT"
__maintainer__ = "Shareef M. Dabdoub"
__email__ = "dabdoub.2@osu.edu"

fontsize = 20
font = {"weight": "bold", "size": fontsize}
mpl.rc("font", **font)


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def translate(rect, x, y, width=1):
    """
    Given four points of a rectangle, translate the
    rectangle to the specified x and y coordinates and,
    optionally, change the width.

    :type rect: list of tuples
    :param rect: Four points describing a rectangle.
    :type x: float
    :param x: The amount to shift the rectangle along the x-axis.
    :type y: float
    :param y: The amount to shift the rectangle along the y-axis.
    :type width: float
    :param width: The amount by which to change the width of the
                  rectangle.
    """
    return ((rect[0][0]+x, rect[0][1]+y), (rect[1][0]+x, rect[1][1]+y),
            (rect[2][0]+x+width, rect[2][1]+y), (rect[3][0]+x+width, rect[3][1]+y))


def load_core_file(core_fp):
    """
    For core OTU data file, returns Genus-species identifier for each data
    entry.
    :type core_fp: str
    :param core_fp: A file containing core OTU data.
    :rtype: str
    :return: Returns genus-species identifier based on identified taxonomical
             level.
    """
    core = {}
    with open(core_fp) as in_f:
        for line in in_f.read().splitlines():
            if not line.startswith("#"):
                otu_id, tax = line.split("\t")
                core[otu_id] = oc.otu_name(ast.literal_eval(tax))
    return core


def load_tsv_core(core_fp, skip_header=False):
    with open(core_fp) as inf:
        data = inf.read().splitlines()
        if skip_header:
            data = data[1:]
        return [line.split("\t")[0].strip() for line in data]


def plot_overlaps(otus, group_otus, group_colors, 
                  out_fp, fig_size=None, title="",
                  filter_common=False, table_fp=None):
    """
    Given a list of OTUs and a number of groups containing subsets of
    the OTU set, plot a presence/absence bar chart showing which species
    belong to which groups.
    
    :type otus: list
    :param otus: A list of OTU identifiers (names or otherwise) ordered
                 by greatest presence across the groups, i.e. those that
                 come first appear in all the groups, next come the OTUs
                 that appear in n-1 groups, etc...
    :type group_otus: OrderedDict
    :param group_otus: A dictionary of OTU identifiers (subset of otus)
                       keyed on group name (for display purposes) in
                       display order (bottom to top).
    :type group_colors: dict
    :param group_colors: Color assignment for each group.
    """
    def sort_order_group(sp):
        """
        Assign a score (for use with sorting) to each OTU based
        on the number of groups they occur in and the order
        within those groups (order priority as set by group_otus).
        """
        count = 0
        rank = 0
        in_prev = True
        max_penalty = len(group_otus)

        for i, grp in enumerate(group_otus):
            if sp in group_otus[grp]:
                count += 1
                if in_prev:
                    rank += 1
            else:
                rank -= max_penalty - i
                in_prev = False

        return count, rank, sp


    if filter_common:
        otus = [otu for otu in otus if sort_order_group(otu)[0] < len(group_otus)]
    otus = sorted(otus, key=sort_order_group, reverse=True)

    shared_table = defaultdict(list)
    
    fig, ax = plt.subplots(figsize=fig_size)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=len(otus), integer=True))

    # rectangle prototype modified for each plot marker
    base = [(0,0),(0,0.5),(0,0.5),(0,0)]
    y_step = 1
    x_step = 2

    bars = []
    bar_colors = []

    for i, grp in enumerate(group_otus):
        for j, otu in enumerate(otus):
            if otu in group_otus[grp]:
                bars.append(translate(base, j*x_step+0.5, i*y_step))
                bar_colors.append(group_colors[grp])
                shared_table[grp].append(otu)
            else:
                shared_table[grp].append("")

    black = (0,0,0,1)

    collection = PolyCollection(
        verts=bars,
        facecolors = bar_colors,
        edgecolors = (black,),
        linewidths = (1,),
        transOffset = ax.transData,
        zorder=3
        )

    ax.add_collection(collection)

    # ax.legend([plt.Rectangle((0, 0), 1, 1, fc=color) for color in group_colors.values()],
    #               group_colors.keys(), loc="best")
    
    # Title
    axttl = ax.title
    axttl.set_position([.5, 1.05])
    ax.set_title(title, {"fontsize": fontsize*1.5, "fontweight": "bold"})
    
    plt.xticks(range(1, len(otus)*x_step, x_step), otus, rotation="vertical")
    plt.yticks([i-0.75 for i in range(1, len(group_otus)*y_step+1, y_step)], 
               group_otus.keys(), rotation="horizontal")

    ax.margins(0.05)
    ax.yaxis.set_visible(True)
    ax.set_xlim((0, len(otus)*x_step))

    if table_fp:
        with open(table_fp, "w") as outf:
            csvw = csv.writer(outf, delimiter="\t")
            # write header
            csvw.writerow(shared_table.keys())
            # write data rows
            for i, _ in enumerate(otus):
                csvw.writerow([shared_table[grp][i] for grp in shared_table])

    # save or display result
    if out_fp:
        plt.savefig(out_fp, facecolors="0.9", edgecolor="none",
                    bbox_inches="tight", pad_inches=0.1)
    else:
        plt.show()



def handle_program_options():
    """Parses the given options passed in at the command line."""
    parser = argparse.ArgumentParser(description="Given a set of core "
                                     "microbiome files, create a matching set "
                                     "of overlapping barplots that visualize "
                                     "the species belonging to each core "
                                     "microbiome.")
    input_grp = parser.add_mutually_exclusive_group(required=True)
    input_grp.add_argument("-i", "--core_files", nargs="+",
                           help="Path to each core microbiome file (i.e. from "
                                "compute_core_microbiome.py) to visualize as "
                                "an overlap plot. NOTE: The files should be "
                                "given in the same order that the groups "
                                "appear in the mapping file.")
    input_grp.add_argument("-tsv", "--tsv_core_files", nargs="+",
                           help="Path to each core microbiome file in TSV "
                                "format. The first column is expected to have "
                                "the OTUs or other names that will be matched "
                                "between the group cores. All other columns "
                                "will be ignored. NOTE: The files should be "
                                "given in the same order that the groups "
                                "appear in the mapping file.")
    parser.add_argument("--skipheader", action="store_true",
                        help="If using TSV files (-tsv) for input, the header "
                             "line will be skipped when reading each file.")
    parser.add_argument("-m", "--map_fp",
                        help="Metadata mapping file.")
    parser.add_argument("-g", "--group_by", required=True,
                        help="A column name in the mapping file containing "
                              "categorical values that will be used to identify"
                              "groups. Each sample ID must have a group entry."
                              "Default is no categories and all the data will"
                              "be treated as a single group. [REQUIRED]")
    parser.add_argument("-c", "--color_by", required=True,
                        help="A column name in the mapping file containing\
                              hexadecimal (#FF0000) color values that will\
                              be used to color the groups. Each sample ID must\
                              have a color entry. [REQUIRED]")
    parser.add_argument("--filtercommon", action="store_true",
                        help="Specifying this option will hide OTUs that \
                              are shared among all groups.")
    parser.add_argument("--title", default="",
                        help="A descriptive title that will appear at the top \
                        of the output plot. Surround with quotes if there are\
                        spaces in the title.")
    parser.add_argument("--figsize", default=[10, 5], type=int, nargs=2,
                        help="Specify the 'width height' in inches for the "
                             "core overlap plot. By default, figure size is "
                             "10x5 inches.")
    parser.add_argument("-o", "--out_fp", default=None,
                        help="The path of the file to save the plot under. "
                              "If specified, the figure will be saved directly "
                              "instead of opening a window in which the plot "
                              "can be viewed before saving.")
    parser.add_argument("--table_fp", default=None,
                        help="Output the list of overlapping OTUs as they "
                             "appear in the figure into a TSV file.")
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

    # map groups to colors
    class_colors = util.color_mapping(imap, header, args.group_by, args.color_by)

    core_files = args.core_files
    tsv = False
    if args.core_files is None:
        core_files = args.tsv_core_files
        tsv = True

    # map each core file to its matching category in the mapping file
    group_cores = OrderedDict()
    for group, fp in zip(class_colors, core_files):
        if not tsv:
            core = load_core_file(fp)
            group_cores[group] = [name.replace("_", " ") for name in core.values()
                                    if not name.startswith("Unclassified")]
        else:
            group_cores[group] = load_tsv_core(fp, args.skipheader)

    # create the overlap set of OTUs and plot
    overlap = set()
    overlap.update(*group_cores.values())

    plot_overlaps(overlap, group_cores, class_colors, 
                  out_fp=args.out_fp, fig_size=args.figsize, title=args.title,
                  filter_common=args.filtercommon,
                  table_fp=args.table_fp)


if __name__ == '__main__':
    main()













