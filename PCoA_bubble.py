#!/usr/bin/env python
'''
Create a series of PCoA plots where the marker size varies by relative 
abundance of a particular OTU

@author: Shareef M Dabdoub
'''
from __future__ import division
from util import parse_map_file
import argparse
from collections import OrderedDict, defaultdict
import json
import sys
import pylab as p

def otu_biom_entry_num(ID, biom, entry_type='rows'):
    for i, entry in enumerate(biom[entry_type]):
        if entry['id'] == ID:
            return i


def rel_abundance(otuID, sampleID, biom):
    otuRow = otu_biom_entry_num(otuID, biom)
    sCol = otu_biom_entry_num(sampleID, biom, 'columns')
    otuAbundance = 0
    sampleAbundance = 0
    
    for row, col, amt in biom['data']:
        if col == sCol:
            sampleAbundance += amt
        if col == sCol and row == otuRow:
            otuAbundance = amt
    
    if sampleAbundance == 0:
        return 0
    
    return otuAbundance / sampleAbundance * 15000



def parse_unifrac(unifracFN):
    """
    Parses the unifrac results file into a dictionary
    
    :@type unifracFN: str
    :@param unifracFN: The path to the unifrac results file
    :@rtype: dict
    :@return: A dictionary with keys: 'pcd' (principle coordinates data) which 
              is a dictionary of the data keyed by sample ID, 
              'eigvals' (eigenvalues), and 'varexp' (variation explained)
    """
    with open(unifracFN) as uF:
        unifrac = {'pcd':{}, 'eigvals':[], 'varexp':[]}

        lines = uF.readlines()[1:]
        for line in lines:
            if line == '\n': break
            line = line.split('\t')
            unifrac['pcd'][line[0]] = line[1:]
        
        unifrac['eigvals'] = lines[-2].split('\t')
        unifrac['varexp'] = lines[-1].split('\t')
        
def link_samples_to_categories(imap, category_idx):
    """
    Creates a dictionary of category types with all the associated Sample IDs
    
    :@type imap: dict
    :@param imap: The mapping data that lists category information
    :@type category_idx: int
    :@param category_idx: The index in the mapping data of the user-specified
                          category
    :@rtype: OrderedDict
    :@return: A dictionary of category types mapped to a list of corresponding 
              sample IDs
    """
    cat_types = defaultdict(list)
    for row in imap.values():
        cat_types[row[category_idx]].append(row[0])
    return OrderedDict(sorted(cat_types.items(), key=lambda t: t[0]))


def rstyle(ax): 
    """
    Styles axes to appear like ggplot2
    Must be called after all plot and axis manipulation operations have been 
    carried out (needs to know final tick spacing)
    messymind.net/2012/07/making-matplotlib-look-like-ggplot/
    """
    # set the style of the major and minor grid lines, filled blocks
    ax.grid(True, 'major', color='w', linestyle='-', linewidth=1.4)
    ax.grid(True, 'minor', color='0.92', linestyle='-', linewidth=0.7)
    ax.patch.set_facecolor('0.85')
    ax.set_axisbelow(True)
    
    # set minor tick spacing to 1/2 of the major ticks
    ax.xaxis.set_minor_locator(p.MultipleLocator((p.plt.xticks()[0][1] - p.plt.xticks()[0][0]) / 2.0))
    ax.yaxis.set_minor_locator(p.MultipleLocator((p.plt.yticks()[0][1] - p.plt.yticks()[0][0]) / 2.0))
    
    # remove axis border
    for child in ax.get_children():
        if isinstance(child, p.matplotlib.spines.Spine):
            child.set_alpha(0)

    # restyle the tick lines
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(5)
        line.set_color("gray")
        line.set_markeredgewidth(1.4)
    
    # remove the minor tick lines    
    for line in ax.xaxis.get_ticklines(minor=True) + ax.yaxis.get_ticklines(minor=True):
        line.set_markersize(0)
    
    # only show bottom left ticks, pointing out of axis
    p.rcParams['xtick.direction'] = 'out'
    p.rcParams['ytick.direction'] = 'out'
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    
    if ax.legend_ <> None:
        lg = ax.legend_
        lg.get_frame().set_linewidth(0)
        lg.get_frame().set_alpha(0.5)


def plot_PCoA(rh, rd, hh, hd, otu_name):
    p.matplotlib.rc('axes', edgecolor='black')
    p.matplotlib.rc('axes', facecolor='grey')
    fig = p.plt.figure()
    fig.set_figheight(8)
    fig.set_figwidth(10)
    
    p.scatter(rd['pc1'], rd['pc2'], rd['size'], color='#90EE90', marker='s')
    p.scatter(rd['zpc1'], rd['zpc2'], s=20, edgecolor='#90EE90', facecolor='none', marker='s')
    p.scatter(rh['pc1'], rh['pc2'], rh['size'], color='#FFA07A', marker='^')
    p.scatter(rh['zpc1'], rh['zpc2'], s=20, edgecolor='#FFA07A', facecolor='none', marker='^')
    p.scatter(hh['pc1'], hh['pc2'], hh['size'], color='green', marker='o')
    p.scatter(hh['zpc1'], hh['zpc2'], s=20, edgecolor='green', facecolor='none', marker='o')
    p.scatter(hd['pc1'], hd['pc2'], hd['size'], color='red', marker='p')
    p.scatter(hd['zpc1'], hd['zpc2'], s=20, edgecolor='red', facecolor='none', marker='p')

    # set plot labels
    rhr = p.Rectangle((0, 0), 1, 1, fc='#90EE90')
    rdr = p.Rectangle((0, 0), 1, 1, fc='#FFA07A')
    hhr = p.Rectangle((0, 0), 1, 1, fc='green')
    hdr = p.Rectangle((0, 0), 1, 1, fc='red')
    p.legend((rhr, rdr, hhr, hdr), ('Rat Health', 'Rat Disease', 'Human Health', 'Human Disease'), loc='upper left')
    p.title(otu_name, style='italic')
    p.xlabel('PC2 ({:.2f}%)'.format(float(unifrac[-1].split('\t')[2])))
    p.ylabel('PC1 ({:.2f}%)'.format(float(unifrac[-1].split('\t')[1])))
    p.hlines(0, -0.4, 0.4)
    p.vlines(0, -0.3, 0.3)
    p.xlim(-0.4, 0.4)
    p.ylim(-0.3, 0.3)
    
    fig.savefig('./figures/' + '_'.join(otu_name.split()) + '.png', facecolor='gray', edgecolor='none')


def handle_program_options():
    parser = argparse.ArgumentParser(description="Create a series of Principle\
                                     Coordinate plots for each OTU in an \
                                     input list where the plot points are \
                                     varied in size by the relative abundance \
                                     of the OTU (relative to either Sample or\
                                     the total contribution of the OTU to the \
                                     data set.")
    parser.add_argument('-i', '--otu_table', required=True,
                        help="The biom-format file with OTU-Sample abundance \
                              data.")
    parser.add_argument('-u', '--unifrac', required=True, 
                        help='Principle coordinates analysis file. \
                              Eg. unweighted_unifrac_pc.txt')
    parser.add_argument('-d', '--names_colors_ids_fn', required=True, 
                        help='The name of an input data file containing three\
                              items: Line 1: a tab-separated list of display\
                              names for the different types in the specified \
                              mapping category (--mapping_category), Line 2:\
                              a matching tab-separated list of hexadecimal\
                              colors for each of the category types, Line 4:\
                              blank, Lines 5-end: a tab-separated pair \
                              specifying OTU ID and OTU Name. Each entry will\
                              get a separate PCoA plot under a file with the \
                              name of the OTU.')
    parser.add_argument('-m', '--mapping', required=True,
                        help="The mapping file specifying group information \
                              for each sample.")
    parser.add_argument('-o', '--output_dir', default='.',
                        help="The directory to output the PCoA plots to.")
    parser.add_argument('-c', '--map_category',
                        help="Any mapping category, such as treatment type, \
                              that will be used to group the data in the \
                              output plots. For example, one category \
                              with three types will result in three different\
                              point sets in the final output.")
    parser.add_argument('-n', '--normalization', default='sample', 
                        choices=['sample','otu'],
                        help="Specifies whether OTU abundance is \
                        normalized column-wise (per-sample) or row-wise \
                        (per-OTU).")
#    parser.add_argument('-v', '--verbose', action='store_true')
    
    return parser.parse_args()


def main():
    args = handle_program_options()
    
    with open(args.otu_table) as bF:
        biom = json.loads(bF.readline())
    
    unifrac = parse_unifrac(args.unifrac)
    
    otus = {}
    with open(args.names_colors_ids_fn, 'rU') as nciF:
        category_names = nciF.readline().split('\t')
        category_colors = nciF.readline().split('\t')
        for line in nciF.readlines()[2:]:
            line = line.split()
            otus[line[0]] = ' '.join(line[1:])
            
    imap = parse_map_file(args.mapping)
    with open(args.mapping) as mF:
        header = mF.readline().split('\t')
        try:
            category_idx = header.index(args.map_category)
        except ValueError: 
            msg = "Error: Specified mapping category '{}' not found"
            sys.exit(msg.format(args.map_category))
    category_ids = link_samples_to_categories(imap, category_idx)
            
    # plot samples based on relative abundance of some OTU ID
    for otuID in otus:
        cat_data = {cat:{'pc1':[], 'pc2':[], 'size':[], 'zpc1':[], 'zpc2':[]} 
                    for cat in category_ids}
    
        for e in pcd:
            if 'RM' in e[0]:
                category = rh
            elif 'RS' in e[0]:
                category = rd
            elif 'HH' in e[0]:
                category = hh
            elif 'HD' in e[0]:
                category = hd
    
            size = rel_abundance(otuID, e[0], biom)
            if size > 0:
                category['pc1'].append(e[1])
                category['pc2'].append(e[2])
                category['size'].append(size)
            else:
                category['zpc1'].append(e[1])
                category['zpc2'].append(e[2])
                
    
        plot_PCoA(rh, rd, hh, hd, otus[otuID])
