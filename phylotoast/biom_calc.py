'''
Created on Feb 19, 2013

:author: Shareef Dabdoub

This module provides methods for calculating various metrics with regards to
each OTU in an input OTU abundance table. This is currently used by iTol.py
to offload the different methods.
'''
from collections import defaultdict
import math


def relative_abundance(biom):
    """
    Calculate the relative abundance of each OTUID in a Sample.

    :type biom: A BIOM file format converted to JSON string format structure.
    :param biom: OTU table format.

    :type sampleIDs: list
    :param sampleIDs: A list of column id's from BIOM format OTU table.

    :rtype: dict
    :return: Returns a keyed on SampleIDs, and the values are dictionaries
             keyed on OTUID's and their values represent the relative
             abundance of that OTUID in that SampleID.
    """
    ra = {item['id']: defaultdict(int) for item in biom['columns']}
    totals = defaultdict(float)

    for row, col, amt in biom['data']:
        otuID = biom['rows'][row]['id']
        sampleID = biom['columns'][col]['id']
        ra[sampleID][otuID] = amt
        totals[sampleID] += amt

    return {sid: {oid: ra[sid][oid] / totals[sid] for oid in ra[sid]}
            for sid in ra}


def mean_otu_pct_abundance(ra, otuIDs):
    """
    Calculate the mean OTU abundance percentage.

    :type ra: Dict
    :param ra: 'ra' refers to a dictionary keyed on SampleIDs, and the values
                are dictionaries keyed on OTUID's and their values represent
                the relative abundance of that OTUID in that SampleID. 'ra' is
                the output of relative_abundance() function.

    :type otuIDs: List
    :param otuIDs: A list of OTUID's for which the percentage abundance needs
                   to be measured.

    :rtype: dict
    :return: A dictionary of OTUID and their percent relative
                      abundance as key/value pair.
    """
    sids = ra.keys()
    otumeans = defaultdict(int)

    for oid in otuIDs:
        otumeans[oid] = sum([ra[sid][oid] for sid in sids
                            if oid in ra[sid]]) / len(sids) * 100

    return otumeans


def MRA(biom):
    """
    Calculate the mean relative abundance.

    :type biom: A BIOM file format converted to JSON string format structure.
    :param biom: OTU table format.

    :type sampleIDs: List
    :param sampleIDs: A list of column id's from BIOM format OTU table.

    :rtype: dict
    :return: A dictionary keyed on OTUID's and their mean relative abundance
             for a given number of sampleIDs.
    """
    ra = relative_abundance(biom)
    otuIDs = {row['id'] for row in biom['rows']}
    return mean_otu_pct_abundance(ra, otuIDs)


def raw_abundance(biom, sampleIDs=None, sample_abd=True):
    """
    Calculate the total number of sequences in each OTU or SampleID.

    :type biom: A BIOM file format converted to JSON string format structure.
    :param biom: OTU table format.

    :type sampleIDs: List
    :param sampleIDs: A list of column id's from BIOM format OTU table. By
                      default, the list has been set to None.

    :type sample_abd: Boolean
    :param sample_abd: A boolean operator to provide output for OTUID's or
                       SampleID's. By default, the output will be provided
                       for SampleID's.

    :rtype: dict
    :return: Returns a dictionary keyed on either OTUID's or SampleIDs
                     and their respective abundance as values.
    """
    results = defaultdict(int)
    for row, col, amt in biom['data']:
        otuID = biom['rows'][row]['id']
        sampleID = biom['columns'][col]['id']

        if sampleIDs is None or sampleID in sampleIDs:
            if sample_abd:
                results[sampleID] += amt
            else:
                results[otuID] += amt

    return results


def transform_raw_abundance(biom, fn=math.log10,
                            sampleIDs=None, sample_abd=True):
    """
    Function to transform the total abundance calculation for each sample ID
    to another format based on user given transformation function.

    :type biom: A BIOM file format converted to JSON string format structure.
    :param biom: OTU table format.

    :param fn: Mathematical function which is used to transform smax to
               another format. By default, the function has been given as
               base 10 logarithm.

    :rtype: dict
    :return: Returns a dictionary similar to output of raw_abundance function
             but with the abundance values modified by the mathematical
             operation. By default, the operation performed on the abundances
             is base 10 logarithm.
    """
    totals = raw_abundance(biom, sampleIDs, sample_abd)
    return {sid: fn(abd) for sid, abd in totals.items()}
