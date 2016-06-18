#! /usr/bin/env python
"""
Abstract: Get unique OTUIDs present in each category of your dataset.
Date: 06/09/2016
"""
import sys
import argparse
import os.path as osp
from collections import defaultdict
from phylotoast import biom_calc as bc, util
try:
    import biom
except ImportError as ie:
    sys.exit("Please install missing module: {}.".format(ie))


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
    rel_abd = bc.relative_abundance(biomfile)
    for sid in rel_abd:
        samples[sid].update([oid for oid, ra in rel_abd[sid].items() if ra > 0])
    return samples


def sample_group(sid, groups):
    """
    Iterate through all categories in an OrderedDict and return category name if SampleID
    present in that category.

    :type sid: str
    :param sid: SampleID from dataset.

    :type groups: OrderedDict
    :param groups: Returned dict from phylotoast.util.gather_categories() function.

    :return type: str
    :return: Category name used to classify `sid`.
    """
    for name in groups:
        if sid in groups[name].sids:
            return name


def combine_sets(*sets):
    """
    Combine multiple sets to create a single larger set.
    """
    combined = set()
    for s in sets:
        combined.update(s)
    return combined


def unique_otuids(groups):
    """
    Get unique OTUIDs of each category.

    :type groups: Dict
    :param groups: {Category name: OTUIDs in category}

    :return type: dict
    :return: Dict keyed on category name and unique OTUIDs as values.
    """
    uniques = {key: set() for key in groups}
    for i, group in enumerate(groups):
        to_combine = groups.values()[:i]+groups.values()[i+1:]
        combined = combine_sets(*to_combine)
        uniques[group] = groups[group].difference(combined)
    return uniques


def write_uniques(path, prefix, uniques):
    """
    Given a path, the method writes out one file for each group name in the
    uniques dictionary with the file name in the pattern

        PATH/prefix_group.txt

    with each file containing the unique OTUIDs found when comparing that group
    to all the other groups in uniques.

    :type path: str
    :param path: Output files will be saved in this PATH.

    :type prefix: str
    :param prefix: Prefix name added in front of output filename.

    :type uniques: dict
    :param uniques: Output from unique_otus() function.
    """
    for group in uniques:
        fp = osp.join(path, "{}_{}.txt".format(prefix, group))
        with open(fp, "w") as outf:
            outf.write("\n".join(uniques[group]))


def handle_program_options():
    parser = argparse.ArgumentParser(description="Parse a BIOM format file and obtain a "
                                     "list of unique OTUIDs found in each category in "
                                     "mapping file.")
    parser.add_argument("input_biom_fp", help="BIOM format file path.")
    parser.add_argument("output_dir", help="Path to save category unique OTUIDs.")
    parser.add_argument("mapping_file", help="Mapping file with category information.")
    parser.add_argument("category_column", help="Column in mapping file specifying the "
                        "category/condition of all samples.")
    parser.add_argument("-p", "--prefix", default="unique", help="Provide specific text "
                        "to prepend the output file names. By default, the 'unique' will "
                        "be added in front of output filenames.")
    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        # Load biom format file
        biomf = biom.load_table(args.input_biom_fp)
    except TypeError as te:
        sys.exit("The data in the path does not appear to be a BIOM format table. "
                 "Error: {}.".format(te))

    # Determine OTUIDs present in each sample
    sample_otus = assign_otu_membership(biomf)

    try:
        # Parse mapping file
        header, imap = util.parse_map_file(args.mapping_file)
    except ValueError as ve:
        sys.exit("Error: {}.".format(ve))

    # Get relevant category information
    group_data = util.gather_categories(imap, header, [args.category_column])

    # Initialize results dict in group_data with {"otuids": set()} for each category
    for group in group_data:
        group_data[group].results["otuids"] = set()

    # Collect all OTUIDs present in each category
    for sid in sample_otus:
        group = sample_group(sid, group_data)
        group_data[group].results["otuids"].update(sample_otus[sid])

    # Create input for unique_otus
    group_otuids = {group: group_data[group].results["otuids"] for group in group_data}

    # Write out unique OTUIDs to file
    write_uniques(args.output_dir, args.prefix, unique_otuids(group_otuids))

if __name__ == "__main__":
    sys.exit(main())
