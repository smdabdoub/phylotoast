'''
Created on Feb 2, 2013

@author: Shareef Dabdoub
'''
from collections import namedtuple, OrderedDict
import os


FASTARecord = namedtuple("FASTA_Record", "id descr data")


def storeFASTA(fastaFNH):
    """
    Parse the records in a FASTA-format file by first reading the entire file
    into memory.

    :type source: list or open file handle
    :param source: The data source from which to parse the FASTA records.
                    Expects the input to resolve to a collection that can be
                    iterated through, such as a list or an open file handle.

    :rtype: tuple
    :return: FASTA records containing entries for id, description and data.
    """
    fasta = file_handle(fastaFNH).read()
    return [
        FASTARecord(rec[0].split()[0], rec[0], ''.join(rec[1:]))
        for rec in (x.strip().split('\n') for x in fasta.split('>')[1:])
        ]


def parseFASTA(fastaFNH):
    """
    Parse the records in a FASTA-format file keeping the file open, and
    reading through one line at a time.

    :type source: list or open file handle
    :param source: The data source from which to parse the FASTA records.
                    Expects the input to resolve to a collection that can be
                    iterated through, such as a list or an open file handle.

    :rtype: tuple
    :return:FASTA records containing entries for id, description and data.
    """
    recs = []
    seq = []
    ID = ''
    descr = ''

    for line in file_handle(fastaFNH):
        line = line.strip()
        if line[0] == ';':
            continue
        if line[0] == '>':
            # conclude previous record
            if seq:
                recs.append(FASTARecord(ID, descr, ''.join(seq)))
                seq = []
            # start new record
            line = line[1:].split()
            ID, descr = line[0], ' '.join(line[1:])
        else:
            seq.append(line)

    # catch last seq in file
    if seq:
        recs.append(FASTARecord(ID, descr, ''.join(seq)))
    return recs


def parse_map_file(mapFNH):
    """
    Opens a QIIME mapping file and stores the contents in a dictionary
    keyed on SampleID (default) or a user-supplied one. The only
    required fields are SampleID, BarcodeSequence, LinkerPrimerSequence
    (in that order), and Description (which must be the final field).

    :type mapFNH: file or str
    :param mapFNH: Either the full path to the map file or an open file
                    handle

    :rtype: dict
    :return: A map associating each line of the mapping file with the
              appropriate sample ID (each value of the map also contains
              the sample ID). An OrderedDict is used so the returned map is
              guaranteed to have the same order as the input file.

    Example data:
    #SampleID BarcodeSequence LinkerPrimerSequence State   Description
    11.V13    ACGCTCGACA      GTTTGATCCTGGCTCAG    Disease Rat_Oral
    """
    m = OrderedDict()
    map_header = None

    with file_handle(mapFNH) as mapF:
        for line in mapF:
            if line.startswith('#SampleID'):
                map_header = line.split('\t')
            if line.startswith('#') or not line:
                    continue
            line = line.strip().split('\t')
            m[line[0]] = line

    return map_header, m


def write_map_file(mapFNH, items, header):
    """
    Given a list of mapping items (in the form described by the
    parse_mapping_file method) and a header line, write each row to the given
    input file with fields separated by tabs.

    :type mapFNH: file or str
    :param mapFNH: Either the full path to the map file or an open file handle

    :type items: list
    :param item: The list of row entries to be written to the mapping file

    :type header: list or str
    :param header: The descriptive column names that are required as the first
                    line of the mapping file

    :rtype: None
    """
    if isinstance(header, list):
        header = '\t'.join(header)

    with file_handle(mapFNH, 'w') as mapF:
        mapF.write(header)
        for row in items:
            mapF.write('\t'.join(row)+'\n')


def parse_taxonomy_table(idtaxFNH):
    """
    Greengenes provides a file each OTU a full taxonomic designation. This
    method parses that file into a map with (key,val) = (OTU, taxonomy).

    :type idtaxFNH: file or str
    :param idtaxFNH: Either the full path to the map file or an open file
                      handle

    :rtype: dict
    :return: A map associating each OTU ID with the taxonomic specifier.
              An OrderedDict is used so the returned map is guaranteed to have
              the same order as the input file.
    """
    idtax = OrderedDict()
    with file_handle(idtaxFNH) as idtxF:
        for line in idtxF:
            ID, tax = line.strip().split('\t')[:2]
            idtax[ID] = tax

    return idtax


def split_phylogeny(p, level='s'):
    """
    Return either the full or truncated version of a QIIME-formatted
    taxonomy string.

    :type p: str
    :param p: A QIIME-formatted taxonomy string: k__Foo; p__Bar; ...

    :type level: str
    :param level: The different level of identification are kingdom (k),
                phylum (p), class (c),order (o), family (f), genus (g)
                and species (s). If level is not provided, the default level
                of identification is species.

    :rtype: str
    :return: A QIIME-formatted taxonomy string up to the classification given
            by param level.
    """
    level = level+'__'
    result = p.split(level)
    return result[0]+level+result[1].split(';')[0]


def ensure_dir(d):
    """
    Check to make sure the supplied directory path does not exist, if so,
    create it.

    :type d: str
    :param d: It is the full path to a directory.

    :return: Does not return anything, but creates a directory path if it
             doesn't exist already.
    """
    if not os.path.exists(d):
        os.makedirs(d)


def file_handle(fnh, mode='rU'):
    """
    Takes either a file path or an open file handle, checks validity and
    returns an open file handle or raises an appropriate Exception.

    :type fnh: str
    :param fnh: It is the full path to a file, or open file handle

    :type mode: str
    :param mode: The way in which this file will be used, for example to read
                 or write or both. By default, file will be opened in rU mode.

    :return: Returns an opened file for appropriate usage.
    """
    handle = None
    if isinstance(fnh, file):
        if fnh.closed:
            raise ValueError('Input file is closed.')
        handle = fnh
    elif isinstance(fnh, str):
        handle = open(fnh, mode)

    return handle

# Meant to contain all the data necessary for calculating a single column of
# an iTol data table
DataCategory = namedtuple('DataCategory', 'sids results')


def gather_categories(imap, header, categories=None):
    """
    Find the user specified categories in the map and create a dictionary
    to contain the relevant data for each type within the categories. Multiple
    categories will have their types combined such that each possible
    combination will have its own entry in the dictionary.

    :type imap: dict
    :param imap: The input mapping file data keyed by SampleID
    :type header: list
    :param header: The header line from the input mapping file. This will be
                   searched for the user-specified categories
    :type categories: list
    :param categories: The list of user-specified categories from the mapping
                        file
    :rtype: dict
    :return: A sorted dictionary keyed on the combinations of all the types
              found within the user-specified categories. Each entry will
              contain an empty DataCategory namedtuple. If no categories are
              specified, a single entry with the key 'default' will be
              returned
    """
    if categories is None:
        return {'default': DataCategory(set(imap.keys()), {})}

    cat_ids = [header.index(cat)
               for cat in categories if cat in header and '=' not in cat]
    conditions = {}
    for cat in categories:
        if '=' in cat and cat.split('=')[0] in header:
            conditions[header.index(cat.split('=')[0])] = cat.split('=')[1]

    if not cat_ids and not conditions:
        return {'default': DataCategory(set(imap.keys()), {})}

    if not cat_ids and conditions:
        sids = {sid for sid, row in imap.iteritems()
                if all([row[c] == conditions[c] for c in conditions])}
        return {'default': DataCategory(sids, {})}

    table = {}
    for sid, row in imap.iteritems():
        if all([row[c] == conditions[c] for c in conditions]):
            key = '_'.join([row[cid] for cid in cat_ids])
            if key not in table:
                table[key] = DataCategory(set(), {})
            table[key].sids.add(sid)

    return OrderedDict(sorted(table.items(), key=lambda t: t[0]))


def parse_unifrac(unifracFN):
    """
    Parses the unifrac results file into a dictionary

    :type unifracFN: str
    :param unifracFN: The path to the unifrac results file

    :rtype: dict
    :return: A dictionary with keys: 'pcd' (principle coordinates data) which
             is a dictionary of the data keyed by sample ID,
             'eigvals' (eigenvalues), and 'varexp' (variation explained)
    """
    with open(unifracFN) as uF:
        unifrac = {'pcd': {}, 'eigvals': [], 'varexp': []}

        lines = uF.readlines()[1:]
        for line in lines:
            if line == '\n':
                break
            line = line.rstrip('\n').split('\t')
            unifrac['pcd'][line[0]] = line[1:]

        unifrac['eigvals'] = lines[-2].split('\t')
        unifrac['varexp'] = lines[-1].split('\t')
        return unifrac
