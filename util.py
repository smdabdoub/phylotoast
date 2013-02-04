'''
Created on Feb 2, 2013

@author: Shareef Dabdoub
'''
from collections import namedtuple
import os

FASTARecord = namedtuple("FASTA_Record", "id descr data")

def storeFASTA(fastaF):
    recs = []; seq = []
    ID = ''; descr = ''
    
    for line in fastaF:
        if line == '' or line[0] == ';':
            continue
        elif line[0] == '>':
                # conclude previous record
                if seq:
                    recs.append(FASTARecord(ID, descr, ''.join(seq)))
                    seq = []
                # start new record
                line = line[1:].split()
                ID, descr = line[0], ' '.join(line[1:])
        else:
                seq.append(line.strip())
    return recs

def parseFASTA(fastaF):
    seq = []
    ID = ''; descr = ''
    
    for line in fastaF:
        if line == '' or line[0] == ';':
            continue
        elif line[0] == '>':
                # conclude previous record
                if seq:
                    yield FASTARecord(ID, descr, ''.join(seq))
                    seq = []
                # start new record
                line = line[1:].split()
                ID, descr = line[0], ' '.join(line[1:])
        else:
                seq.append(line.strip())
                
def parse_map_file(mapFN, keyColumn=0):
    """
    Opens a QIIME mapping file and stores the contents in a dictionary
    keyed on SampleID (default) or a user-supplied one. The only 
    required fields are SampleID, BarcodeSequence, LinkerPrimerSequence 
    (in that order), and Description (which must be the final field).
    
    :@param mapFN: Full path to the map file
    :@param keyColumn: The column to key the returned dictionary on.
                       If keyColumn is larger than the number of columns,
                       it will be set to the default, 0.     
    
    Example data:
    #SampleID    BarcodeSequence    LinkerPrimerSequence    Treatment    SampleType    DOB    Description
    11.V13    ACGCTCGACA    GTTTGATCCTGGCTCAG    V13    Rat_Oral_Disease    111111    Rat_Oral
    """
    m = {}
    
    with open(mapFN) as mapfile:
        lines = mapfile.readlines()
        if keyColumn >= len(lines[0].split()): keyColumn = 0
        for line in lines[1:]:
            line = line.strip().split('\t')
            m[line[keyColumn]] = line
            
    return m


def ensure_dir(d):
    """
    Check to make sure the supplied directory path does not exist, if so, 
    create it.
    """
    if not os.path.exists(d):
        os.makedirs(d)










