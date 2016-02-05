from __future__ import division
# python libs
import sys
# 3rd party
importerrors = []
try:
    import statsmodels.nonparametric.kde as kde
except ImportError:
    importerrors.append('statsmodels')
if len(importerrors) != 0:
    for item in importerrors:
        print 'Import Error. Please install missing module:', item
    sys.exit()


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
    e = kde.KDEUnivariate(data)
    e.fit()
    ax.plot(e.support, e.density, color=color, alpha=0.9, linewidth=2.25)
    if fill_bt:
        ax.fill_between(e.support, e.density, alpha=.35, zorder=1,
                        antialiased=True, color=color)
    if title is not None:
        t = ax.set_title(title)
        t.set_y(1.05)
