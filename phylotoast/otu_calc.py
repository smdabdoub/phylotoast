from __future__ import division
import ast
from collections import defaultdict
from phylotoast import biom_calc as bc


def otu_name(tax):
    """
    Determine a simple Genus-species identifier for an OTU, if possible.
    If OTU is not identified to the species level, name it as
    Unclassified (familly/genus/etc...).

    :type tax: list
    :param tax: QIIME-style taxonomy identifiers, e.g.
                 ["k__Bacteria", u"p__Firmicutes", u"c__Bacilli", ...

    :rtype: str
    :return: Returns genus-species identifier based on identified taxonomical
             level.
    """
    extract_name = lambda lvl: "_".join(lvl.split("_")[2:])
    spname = "spp."
    for lvl in tax[::-1]:
        if len(lvl) <= 3:
            continue
        if lvl.startswith("s"):
            spname = extract_name(lvl)
        elif lvl.startswith("g"):
            return "{}_{}".format(extract_name(lvl), spname)
        else:
            return "Unclassified_{}".format(extract_name(lvl))

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
    with open(core_fp, "rU") as in_f:
        return {otu_name(ast.literal_eval(line.split("\t")[1]))
                for line in in_f.readlines()[1:]}

def assign_otu_membership(biomfile):
    """
    Determines the OTUIDs present in each sample.
    :type biomfile: biom.table.Table
    :param biomfile: BIOM table object from the biom-format library.
    :rtype: dict
    :return: Returns a dictionary keyed on Sample ID with sets containing
    the IDs of OTUIDs found in each sample.
    """
    samples = defaultdict(set)
    _ = biomfile.pa()
    for sid in biomfile.ids():
        for otuid in biomfile.ids("observation"):
            if biomfile.get_value_by_ids(otuid, sid) == 1:
                samples[sid].add(otuid)
    return samples
