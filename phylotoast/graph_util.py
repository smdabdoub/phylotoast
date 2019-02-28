from __future__ import division
# python libs
import sys
# 3rd party
importerrors = []
try:
    import numpy as np
except ImportError as ie:
    importerrors.append(ie)
try:
    import statsmodels.nonparametric.kde as kde
except ImportError as ie:
    importerrors.append(ie)
try:
    import matplotlib as mpl
except ImportError as ie:
    importerrors.append(ie)

if len(importerrors) != 0:
    for item in importerrors:
        print('Import Error:', item)
    sys.exit()

from matplotlib.ticker import FuncFormatter, MaxNLocator, MultipleLocator
import matplotlib.pyplot as plt


def plot_kde(data, ax, title=None, color='r', fill_bt=True):
    """
    Plot a smoothed (by kernel density estimate) histogram.
    :type data: numpy array
    :param data: An array containing the data to be plotted

    :type ax: matplotlib.Axes
    :param ax: The Axes object to draw to

    :type title: str
    :param title: The plot title

    :type color: str
    :param color: The color of the histogram line and fill. Note that the fill
                  will be plotted with an alpha of 0.35.

    :type fill_bt: bool
    :param fill_bt: Specify whether to fill the area beneath the histogram line
    """
    if isinstance(data, list):
        data = np.asarray(data)
    e = kde.KDEUnivariate(data.astype(np.float))
    e.fit()
    ax.plot(e.support, e.density, color=color, alpha=0.9, linewidth=2.25)
    if fill_bt:
        ax.fill_between(e.support, e.density, alpha=.35, zorder=1,
                        antialiased=True, color=color)
    if title is not None:
        t = ax.set_title(title)
        t.set_y(1.05)

def ggplot2_style(ax):
    """
    Styles an axes to appear like ggplot2
    Must be called after all plot and axis manipulation operations have been
    carried out (needs to know final tick spacing)
    """
    #set the style of the major and minor grid lines, filled blocks
    ax.grid(True, 'major', color='w', linestyle='-', linewidth=1.4)
    ax.grid(True, 'minor', color='0.92', linestyle='-', linewidth=0.7)
    ax.patch.set_facecolor('0.85')
    ax.set_axisbelow(True)

    #set minor tick spacing to 1/2 of the major ticks
    ax.xaxis.set_minor_locator(MultipleLocator( (plt.xticks()[0][1]-plt.xticks()[0][0]) / 2.0 ))
    ax.yaxis.set_minor_locator(MultipleLocator( (plt.yticks()[0][1]-plt.yticks()[0][0]) / 2.0 ))

    #remove axis border
    for child in ax.get_children():
        if isinstance(child, mpl.spines.Spine):
            child.set_alpha(0)

    #restyle the tick lines
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(5)
        line.set_color("gray")
        line.set_markeredgewidth(1.4)

    #remove the minor tick lines
    for line in ax.xaxis.get_ticklines(minor=True) + ax.yaxis.get_ticklines(minor=True):
        line.set_markersize(0)

    #only show bottom left ticks, pointing out of axis
    mpl.rcParams['xtick.direction'] = 'out'
    mpl.rcParams['ytick.direction'] = 'out'
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')


    if ax.legend_ != None:
        lg = ax.legend_
        lg.get_frame().set_linewidth(0)
        lg.get_frame().set_alpha(0.5)
