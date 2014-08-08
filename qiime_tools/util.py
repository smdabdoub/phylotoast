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
    """
    fasta = file_handle(fastaFNH).read()
    return [FASTARecord(rec[0].split()[0], rec[0], ''.join(rec[1:])) for rec in (x.strip().split('\n') for x in fasta.split('>')[1:])]

def parseFASTA(fastaFNH):
    """
    Parse the records in a FASTA-format file keeping the file open, and reading
    through one line at a time.
    
    :@type source: list or open file handle
    :@param source: The data source from which to parse the FASTA records. 
                    Expects the input to resolve to a collection that can be 
                    iterated through, such as a list or an open file handle.
                    
    :@param 
    """
    recs = []; seq = []
    ID = ''; descr = ''
    
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
    
    :@type mapFNH: file or str
    :@param mapFNH: Either the full path to the map file or an open file 
                    handle
    
    :@rtype: dict
    :@return: A map associating each line of the mapping file with the
              appropriate sample ID (each value of the map also contains 
              the sample ID). An OrderedDict is used so the returned map is
              guaranteed to have the same order as the input file.
    
    Example data:
    #SampleID BarcodeSequence LinkerPrimerSequence State   Description
    11.V13    ACGCTCGACA      GTTTGATCCTGGCTCAG    Disease Rat_Oral
    """
    m = OrderedDict()
    
    with file_handle(mapFNH) as mapF:
        for line in mapF:
            if line.startswith('#') or not line:
                    continue
            line = line.strip().split('\t')
            m[line[0]] = line
            
    return m


def write_map_file(mapFNH, items, header):
    """
    Given a list of mapping items (in the form described by the  
    parse_mapping_file method) and a header line, write each row to the given 
    input file with fields separated by tabs.
    
    :@type mapFNH: file or str
    :@param mapFNH: Either the full path to the map file or an open file handle
    :@type items: list
    :@param item: The list of row entries to be written to the mapping file
    :@type header: list or str
    :@param header: The descriptive column names that are required as the first
                    line of the mapping file
    :@rtype: None
    """
    if isinstance(header, list):
        header = '\t'.join(header)+'\n'
    
    with file_handle(mapFNH, 'w') as mapF:
        mapF.write(header)
        for row in items:
            mapF.write('\t'.join(row)+'\n')


def parse_taxonomy_table(idtaxFNH):
    """
    Greengenes provides a file each OTU a full taxonomic designation. This 
    method parses that file into a map with (key,val) = (OTU, taxonomy).
    
    :@type idtaxFNH: file or str
    :@param idtaxFNH: Either the full path to the map file or an open file 
                      handle
    :@rtype: dict
    :@return: A map associating each OTU ID with the taxonomic specifier.
              An OrderedDict is used so the returned map is guaranteed to have
              the same order as the input file.
    """
    idtax = OrderedDict()
    with file_handle(idtaxFNH) as idtxF:
        for line in idtxF:
            ID, tax = line.strip().split('\t')
            idtax[ID] = tax
    
    return idtax


def split_phylogeny(p, level='s'):
    level = level+'__'
    result = p.split(level)
    return result[0]+level+result[1].split(';')[0]


def ensure_dir(d):
    """
    Check to make sure the supplied directory path does not exist, if so, 
    create it.
    """
    if not os.path.exists(d):
        os.makedirs(d)


def file_handle(fnh, mode='rU'):
    """
    Takes either a file path or an open file handle, checks validity and 
    returns an open file handle or raises an appropriate Exception
    """
    handle = None
    if isinstance(fnh, file):
        if fnh.closed: raise ValueError('Input file is closed.')
        handle = fnh
    elif isinstance(fnh, str):
        handle = open(fnh, mode)
    
    return handle





