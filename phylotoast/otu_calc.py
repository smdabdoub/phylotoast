from __future__ import division
# python libs
import ast
from collections import (OrderedDict,defaultdict)
import json
import math
# 3rd party
#from fuzz import FuzzySet, FuzzyElement
import numpy as np
# local
import util
import biom_calc as bc


# def fuzzy_lookup(orig, keys):
#     """
#     Return the intersection of a fuzzy set and a collection of keys
#     (presumably a subset).

#     :type orig: set
#     :param orig: FuzzySet of SampleID with OTUID and relative abundances.

#     :type keys: list
#     :param keys: Genus-species taxonomic identifier.

#     :rtype: set
#     :return: Returns a new FuzzySet of genus-species identifier and relative
#              abundance for the given list of keys.
#     """
#     fset = FuzzySet()
#     for k in keys:
#         if k in orig:
#             fset.add(orig[k])
#     return fset


# def sdi(fset):
#     """
#     Calculate the Shannon Diversity Index.
#     H = -sum(p*ln(p))
#     where p is the relative abundance of a single OTU in the set.
#     Note that the Equitability Index ( $E_H = H / H_{max}$ ) could be easily
#     calculated from the returned array by:
#     >>>diversities = sdi(fset)
#     >>>equitabilities = diversities/max(diversities)

#     :type fset: FuzzySet
#     :param fset: The set of OTUs and their relative abundance values

#     :rtype: float
#     :return: The Shannon Diversity Index.
#     """
#     p = [e.mu for e in fset]
#     return -1*sum([entry*math.log(entry) for entry in p])


def otu_name_biom(biom_row):
    """
    Given an OTU row from a BIOM table, determine a Genus-species identifier
    from the taxonomic specifier (see otu_name() method)

    :type biom_row: str
    :param biom_row: Row entry of a BIOM format file containing full taxonomy.

    :rtype: str
    :return: Returns the genus-species identifier.
    """
    return otu_name(biom_row['metadata']['taxonomy'])


def otu_name(tax):
    """
    Determine a simple Genus-species identifier for an OTU, if possible.
    If OTU is not identified to the species level, name it as
    Unclassified (familly/genus/etc...).

    :type tax: list
    :param tax: QIIME-style taxonomy identifiers, e.g.
                 ['k__Bacteria', u'p__Firmicutes', u'c__Bacilli', ...

    :rtype: str
    :return: Returns genus-species identifier based on identified taxonomical
             level.
    """
    extract_name = lambda lvl: '_'.join(lvl.split('_')[2:])
    for i, lvl in enumerate(tax):
        lvl = lvl.strip()
        if i < len(tax) - 1 and len(tax[i + 1].strip()) == 3:
            if tax[i].strip()[0] == 'g':
                return extract_name(lvl) + '_spp.'
            else:
                return 'Unclassified_' + extract_name(lvl)
        elif i == len(tax) - 1:
            name = extract_name(lvl)
            if lvl[0] == 's':
                name = extract_name(tax[i-1]) + '_' + name
            return name


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
    with open(core_fp, 'rU') as in_f:
        return {otu_name(ast.literal_eval(line.split('\t')[1]))
                for line in in_f.readlines()[1:]}


# def assign_otu_membership(biom):
#     """
#     Determines the number of OTUs associated with samples using
#     fuzzy sets with membership amount determined by relative abundance.

#     :type biom: str
#     :param biom: BIOM format file converted to Python object using JSON
#                  decoder.

#     :rtype: dict
#     :return: Returns a dictionary of FuzzySet of SampleID's with OTUID and
#              relative abundance as its elements.
#     """
#     samples = defaultdict(FuzzySet)
#     otu_map = {row['id']: row['metadata']['taxonomy']
#                for row in biom['rows']}
#     rel_abd = bc.relative_abundance(biom, [col['id']
#                                     for col in biom['columns']])
#     for sid in rel_abd:
#         samples[sid].update([FuzzyElement(otu_name(otu_map[oid]), ra)
#                             for oid, ra in rel_abd[sid].items()])

#     return samples


# def print_membership(entry):
#     """
#     Given an entry from a Sample dictionary (see assign_otu_membership) of
#     fuzzy otu membership, pretty-print the members as percentages.

#     :type entry: list
#     :param entry: SampleID's from the output dict of assign_otu_membership().

#     :rtype: str
#     :return: Returns OTU name and percentage relative abundance as membership
#              for the given list of SampleID's.
#     """
#     data = [str(e).split(' \\ ') for e in entry]
#     data = [(name, float(mu)) for name, mu in data]
#     data = ['{0}: {1:4.2%}'.format(name, mu)
#             for name, mu in sorted(data, key=lambda x:x[1])]
#     for e in data:
#         print e
