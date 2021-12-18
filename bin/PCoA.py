#!/usr/bin/env python
import argparse
from collections import OrderedDict
import itertools
import sys
from phylotoast import util, graph_util as gu
errors = []
try:
    from palettable.colorbrewer.qualitative import Set3_12
except ImportError as ie:
    errors.append("No module named palettable")
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
except ImportError as ie:
    errors.append(ie)
if len(errors) != 0:
    for item in errors:
        print("Import Error:", item)
    sys.exit()


def handle_program_options():
    parser = argparse.ArgumentParser(description="Create a 2D or 3D PCoA plot. By default"
                                     ", this script opens a window with the plot "
                                     "displayed if you want to change certain aspects of "
                                     "the plot (such as rotate the view in 3D mode). If "
                                     "the -o option is specified, the plot will be saved "
                                     "directly to an image without the initial display "
                                     "window.")
    parser.add_argument("-i", "--coord_fp", required=True,
                        help="Input principal coordinates filepath (i.e. resulting file "
                             "from principal_coordinates.py) [REQUIRED].")
    parser.add_argument("-m", "--map_fp", required=True,
                        help="Input metadata mapping filepath [REQUIRED].")
    parser.add_argument("-g", "--group_by", required=True,
                        help="Any mapping categories, such as treatment type, that will "
                        "be used to group the data in the output iTol table. For example,"
                        " one category with three types will result in three data columns"
                        " in the final output. Two categories with three types each will "
                        "result in six data columns. Default is no categories and all the"
                        " data will be treated as a single group.")
    parser.add_argument("-d", "--dimensions", default=2, type=int, choices=[2, 3],
                        help="Choose whether to plot 2D or 3D.")
    parser.add_argument("-c", "--colors", default=None,
                        help="A column name in the mapping file containing hexadecimal "
                        "(#FF0000) color values that will be used to color the groups. "
                        "Each sample ID must have a color entry.")
    parser.add_argument("-s", "--point_size", default=100, type=int,
                        help="Specify the size of the circles representing each of the "
                        "samples in the plot")
    parser.add_argument("--pc_order", default=[1, 2], type=int, nargs=2,
                        help="Choose which Principle Coordinates are displayed and in "
                        "which order, for example: 1 2. This option is only used when a "
                        "2D plot is specified.")
    parser.add_argument("--x_limits", type=float, nargs=2,
                        help="Specify limits for the x-axis instead of automatic setting "
                        "based on the data range. Should take the form: --x_limits -0.5 "
                        "0.5")
    parser.add_argument("--y_limits", type=float, nargs=2,
                        help="Specify limits for the y-axis instead of automatic setting "
                        "based on the data range. Should take the form: --y_limits -0.5 "
                        "0.5")
    parser.add_argument("--z_limits", type=float, nargs=2,
                        help="Specify limits for the z-axis instead of automatic setting "
                        "based on the data range. Should take the form: --z_limits -0.5 "
                        "0.5")
    parser.add_argument("--z_angles", type=float, nargs=2, default=[-134.5, 23.],
                        help="Specify the azimuth and elevation angles for a 3D plot.")
    parser.add_argument("-t", "--title", default="", help="Title of the plot.")
    parser.add_argument("--figsize", default=[14, 8], type=int, nargs=2,
                        help="Specify the 'width height' in inches for PCoA plots. "
                        "Default figure size is 14x8 inches")
    parser.add_argument("--font_size", default=12, type=int,
                        help="Sets the font size for text elements in the plot.")
    parser.add_argument("--label_padding", default=15, type=int,
                        help="Sets the spacing in points between the each axis and its "
                        "label.")
    parser.add_argument("--legend_loc", default="best", choices=['best','upper right','upper left',
                        'lower left', 'lower right', 'right', 'center left', 'center right',
                        'lower center', 'upper center', 'center', 'outside'], 
                        help="Sets the location of the Legend. Default: best.")
    parser.add_argument("--annotate_points", action="store_true",
                        help="If specified, each graphed point will be labeled with its "
                        "sample ID.")
    parser.add_argument("--ggplot2_style", action="store_true",
                        help="Apply ggplot2 styling to the figure.")
    parser.add_argument("-o", "--out_fp", default=None,
                        help="The path and file name to save the plot under. If specified"
                        ", the figure will be saved directly instead of opening a window "
                        "in which the plot can be viewed before saving.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.coord_fp):
            pass
    except IOError as ioe:
        err_msg = "\nError in input principal coordinates filepath (-i): {}\n"
        sys.exit(err_msg.format(ioe))

    try:
        with open(args.map_fp):
            pass
    except IOError as ioe:
        err_msg = "\nError in input metadata mapping filepath (-m): {}\n"
        sys.exit(err_msg.format(ioe))

    with open(args.coord_fp) as F:
        pcd = F.readlines()
    pcd = [line.split("\t") for line in pcd]

    map_header, imap = util.parse_map_file(args.map_fp)

    data_gather = util.gather_categories(imap, map_header,
                                         args.group_by.split(","))
    categories = OrderedDict([(condition, {"pc1": [], "pc2": [], "pc3": []})
                              for condition in data_gather.keys()])

    bcolors = itertools.cycle(Set3_12.hex_colors)
    if not args.colors:
        colors = [bcolors.next() for _ in categories]
    else:
        colors = util.color_mapping(imap, map_header,
                                    args.group_by, args.colors)
        colors = colors.values()

    parsed_unifrac = util.parse_unifrac(args.coord_fp)

    pco = args.pc_order
    if args.dimensions == 3:
        pco.append(3)

    pc1v = parsed_unifrac["varexp"][pco[0] - 1]
    pc2v = parsed_unifrac["varexp"][pco[1] - 1]
    if args.dimensions == 3:
        pc3v = parsed_unifrac["varexp"][pco[2] - 1]

    for sid, points in parsed_unifrac["pcd"].items():
        for condition, dc in data_gather.items():
            if sid in dc.sids:
                cat = condition
                break
        categories[cat]["pc1"].append((sid, points[pco[0] - 1]))
        categories[cat]["pc2"].append((sid, points[pco[1] - 1]))

        if args.dimensions == 3:
            categories[cat]["pc3"].append((sid, points[pco[2] - 1]))

    axis_str = "PC{} (Percent Explained Variance {:.3f}%)"
    # initialize plot
    fig = plt.figure(figsize=args.figsize)
    if args.dimensions == 3:
        ax = fig.add_subplot(111, projection="3d")
        ax.view_init(elev=args.z_angles[1], azim=args.z_angles[0])
        ax.set_zlabel(axis_str.format(3, pc3v), labelpad=args.label_padding)
        if args.z_limits:
            ax.set_zlim(args.z_limits)
    else:
        ax = fig.add_subplot(111)

    # plot data
    for i, cat in enumerate(categories):
        if args.dimensions == 3:
            ax.scatter(xs=[e[1] for e in categories[cat]["pc1"]],
                       ys=[e[1] for e in categories[cat]["pc2"]],
                       zs=[e[1] for e in categories[cat]["pc3"]],
                       zdir="z", c=colors[i], s=args.point_size, label=cat,
                       edgecolors="k")
        else:
            ax.scatter([e[1] for e in categories[cat]["pc1"]],
                       [e[1] for e in categories[cat]["pc2"]],
                       c=colors[i], s=args.point_size, label=cat, edgecolors="k")

        # Script to annotate PCoA sample points.
        if args.annotate_points:
            for x, y in zip(categories[cat]["pc1"], categories[cat]["pc2"]):
                ax.annotate(
                    x[0], xy=(x[1], y[1]), xytext=(-10, -15),
                    textcoords="offset points", ha="center", va="center",
                    )

    # customize plot options
    if args.x_limits:
        ax.set_xlim(args.x_limits)
    if args.y_limits:
        ax.set_ylim(args.y_limits)

    ax.set_xlabel(axis_str.format(pco[0], float(pc1v)), labelpad=args.label_padding)
    ax.set_ylabel(axis_str.format(pco[1], float(pc2v)), labelpad=args.label_padding)

    if args.legend_loc == "outside":
        leg = plt.legend(bbox_to_anchor=(1.05, 0.94), loc=2, borderaxespad=-2.0, 
                         scatterpoints=3, frameon=True, framealpha=1)
        # https://stackoverflow.com/a/45846024
        plt.gcf().canvas.draw()
        invFigure = plt.gcf().transFigure.inverted()
        lgd_pos = leg.get_window_extent()
        lgd_coord = invFigure.transform(lgd_pos)
        lgd_xmax = lgd_coord[1, 0]
        ax_pos = plt.gca().get_window_extent()
        ax_coord = invFigure.transform(ax_pos)
        ax_xmax = ax_coord[1, 0]
        shift = 1 - (lgd_xmax - ax_xmax)
        plt.gcf().tight_layout(rect=(0, 0, shift, 1))
    else:
        leg = plt.legend(loc=args.legend_loc, scatterpoints=3, frameon=True, framealpha=1)
    leg.get_frame().set_edgecolor('k')

    # Set the font characteristics
    font = {"family": "normal", "weight": "bold", "size": args.font_size}
    mpl.rc("font", **font)

    if args.title:
        ax.set_title(args.title)

    if args.ggplot2_style and not args.dimensions == 3:
        gu.ggplot2_style(ax)

    # save or display result
    if args.out_fp:
        fig.savefig(args.out_fp, facecolor="white", edgecolor="none", bbox_inches="tight",
                    pad_inches=0.2)
    else:
        plt.show()


if __name__ == "__main__":
    sys.exit(main())
