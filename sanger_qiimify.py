'''
Created on Sep 19, 2012

@author: Shareef Dabdoub

Convert Sanger Sequencing data into a format useable by QIIME. 
Specifically, generate barcodes for each sample (plate), modify 
the FASTA files to include the barcode and a placeholder primer, 
concatenate all the resulting files, and generate a mapping file 
for the sample IDs and the generated barcodes (and primers). 
'''
import argparse
from collections import namedtuple
import os
import os.path as osp
import sys

from Bio import SeqIO


file_types = {'fasta': ('.fasta', '.fa', '.txt', '.seq', '.Seq'),
              'fastq': ('.fastq'),
              'abi': ('.abi','.ab1'),
              'qual': ('.qual')}

MapRecord = namedtuple('MapRecord', 'barcode primer descr')

def gather_hierarchical_sample_ids(sampleDir):
    """
    Given a directory containing directories of Sanger sequenced 
    DNA, generate sample IDs from the directory names.
    
    :@param sampleDir: The top-level directory from which to begin the search.
    :@param write: A function that will modify and write the contents of the 
                   fasta file out to the new aggregate file. 
    """
    sampleIDs = []
    dirs = []
    for item in os.listdir(sampleDir):
        iPath = os.path.join(sampleDir, item)
        if os.path.isdir(iPath):
            dirInclude = False
            for sitem in os.listdir(iPath):
                if (osp.isfile(os.path.join(iPath, sitem)) and 
                    osp.splitext(sitem)[1] in file_types['fasta']):
                    try:
                        SeqIO.read(os.path.join(iPath, sitem), 'fasta')
                        dirInclude = True
                        dirs.append(iPath)
                        break
                    except ValueError, _:
                        pass
            if dirInclude:
                sampleIDs.append(item)
    
    return sampleIDs, dirs

def generate_barcodes(nIds, codes=None, idx=-1, codeLen=12):
    """
    Given a list of sample IDs generate unique n-base barcodes for each.
    Note that only 4^n unique barcodes are possible. 
    """
    if codes is None: codes = [''.join(['A' for _ in range(codeLen)])]
    if idx < -codeLen: return
    if len(set(codes)) > nIds: codes.pop(); return
    
    def next_code(b,c,i):
        c.append(c[-1][:i] + b + (c[-1][i+1:] if i<-1 else ''))
        return c
    
    generate_barcodes(nIds, next_code('A',codes,idx), idx-1)
    generate_barcodes(nIds, next_code('C',codes,idx), idx-1)
    generate_barcodes(nIds, next_code('T',codes,idx), idx-1)
    generate_barcodes(nIds, next_code('G',codes,idx), idx-1)
    
    return list(set(codes))[:nIds]
    

def write_mapping_file(mapFN, sampleIDs, barcodes):
    """
    Given a mapping from sample IDs to barcodes/primers/other info, 
    write out a QIIME-compatible mapping file.
    
    File format described at: http://qiime.org/documentation/file_formats.html
    Note that primer column can be left blank.
    """
    header = '\t'.join(['#SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
                        'Description', '\n'])
    primer = 'GGGGGGGGGGGGGGGGGGGG'
    sampleMap = {}
    
    with open(mapFN, 'w') as mapfile:
        mapfile.write(header)
        for sid, bc in zip(sampleIDs, barcodes):
            sampleMap[sid] = MapRecord(bc, primer, sid)
            mapfile.write('\t'.join([sid, bc, primer, sid, '\n']))
    
    return sampleMap

def scrobble_fasta_dir(dataDir, sampleMap, outF, sampleID=None):
    """
    Given a sample ID and a mapping, modify a Sanger FASTA file 
    to include the barcode and 'primer' in the sequence data 
    and change the description line as needed.
    """
    for item in os.listdir(dataDir):
        if (osp.isfile(os.path.join(dataDir, item)) and 
            osp.splitext(item)[1] in file_types['fasta']):
            for record in SeqIO.parse(os.path.join(dataDir, item), 'fasta'):
                record.seq = ''.join(sampleMap[sampleID].barcode,
                                     sampleMap[sampleID].primer, record.seq)
                SeqIO.write(record, outF, 'fasta')
            

def handle_program_options():
    """
    Uses the built-in argparse module to handle command-line options for the 
    program.
     
    :@return: The gathered command-line options specified by the user 
    :@rtype: argparse.ArgumentParser 
    """
    parser = argparse.ArgumentParser(description="Convert Sanger-sequencing \
                                     derived data files for use with the \
                                     metagenomics analysis program QIIME, by \
                                     extracting Sample ID information, adding\
                                     barcodes and primers to the sequence \
                                     data, and outputting a mapping file and\
                                     single FASTA-formatted sequence file \
                                     formed by concatenating all input data.") 
                                     
    parser.add_argument('input_dir', help="The directory containing sequence \
                        data files. Assumes all data files are placed in this \
                        directory. For files organized within folders by \
                        sample, use -s in addition.")
    parser.add_argument('output_map_file', type=argparse.FileType('w'),
                        help="QIIME-formatted mapping file linking Sample IDs \
                              with barcodes and primers.")
    parser.add_argument('output_sequence_file', type=argparse.FileType('w'), 
                        help="Single file containing all sequence data found \
                              in input_dir, FASTA-formatted with barcode and \
                              primer preprended to sequence.")
    parser.add_argument('-b', '--barcode_length', type=int, default=12,
                        help="Length of the generated barcode sequences. \
                              Default is 12 (QIIME default), minimum is 8.")
    parser.add_argument('-i', '--identifier_pattern', type=tuple, default=(),
                        help="Indicates how to extract the Sample ID from the \
                              description line. Specify two things: 1. Field \
                              separator, 2. Field number of Sample ID. If the \
                              separator is a space or tab, use \s or \\t \
                              respectively. Example: >ka-SampleID-2091, use \
                              -i - 2, indicating - is the separator and the \
                              Sample ID is field #2.")
    parser.add_argument('-s', '--subdirectories', 
                        help="Specify that sequence data files are contained \
                              in subfolders named by Sample ID.",
                        action='store_true')
    return parser.parse_args(sys.argv)
    

    
def main():
    args = handle_program_options()
    pass
#    sampleIDs, dirs = gather_hierarchical_sample_ids(sys.argv[1])
#    barcodes = generate_barcodes(len(sampleIDs))
#    sampleMap = write_mapping_file(sys.argv[2], sampleIDs, barcodes)
#    outF = open(sys.argv[3], 'w')
#    
#    for dataDir, sampleID in zip(dirs,sampleIDs):
#        scrobble_fasta_dir(dataDir, sampleMap, outF, sampleID)
#    
#    outF.close()


if __name__ == '__main__':
    main()












