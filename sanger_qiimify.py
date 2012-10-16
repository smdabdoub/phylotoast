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
import codecs
import os
import os.path as osp
import sys

from Bio import SeqIO


file_types = {'fasta': ('.fasta', '.fa', '.txt', '.seq', '.Seq'),
              'fastq': ('.fastq'),
              'abi': ('.abi','.ab1'),
              'qual': ('.qual')}

MapRecord = namedtuple('MapRecord', 'barcode primer descr')
IDPattern = namedtuple('IDPattern','separator field')

class ValidateIDPattern(argparse.Action):
    def __call__(self, parser, args, values, option_string = None):
        # print '{n} {v} {o}'.format(n=args, v=values, o=option_string)
        separator, field = values
        try:
            field = int(field)
        except ValueError:
            print parser.format_usage().strip()
            print '%s: error: FIELD_NUMBER must be a number' % sys.argv[0]
            sys.exit(-1)
        
        setattr(args, self.dest, IDPattern(separator, field))


def open_enc(fn, utf16):
    """
    Open a file for reading with either the standard built-in open() or the
    open in the codecs module for reading UTF-16 encoded files.
    
    :@param fn: Name of the file to open
    :@param utf16: Boolean flag to specify a UTF-16 file
    """
    if utf16:
        return codecs.open(fn, encoding='utf-16')
    return open(fn, 'rU')


def gather_hierarchical_sample_ids(sampleDir, utf16):
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
                        with open_enc(os.path.join(iPath, sitem), utf16) as fh:
                            SeqIO.read(fh, 'fasta')
                            dirInclude = True
                            dirs.append(iPath)
                            break
                    except ValueError:
                        pass
            if dirInclude:
                sampleIDs.append(item)
    
    return sampleIDs, dirs


def gather_sample_ids(sampleDir, idPattern, utf16):
    """
    
    """
    sampleIDs = []
    separator, field = idPattern
    for item in os.listdir(sampleDir):
        iPath = os.path.join(sampleDir, item)
        if (osp.isfile(iPath) and 
            osp.splitext(item)[1] in file_types['fasta']):
            try:
                with open_enc(iPath, utf16) as fh:
                    record = SeqIO.parse(fh, 'fasta').next()
                    sampleIDs.append(record.id.split(separator)[field-1])
            except StopIteration:
                print 'Empty FASTA file: %s' % iPath
    
    return sampleIDs


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
    

def write_mapping_file(mapF, sampleIDs, barcodes):
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
    
    with mapF:
        mapF.write(header)
        for sid, bc in zip(sampleIDs, barcodes):
            sampleMap[sid] = MapRecord(bc, primer, sid)
            mapF.write('\t'.join([sid, bc, primer, sid, '\n']))
    
    return sampleMap

def scrobble_fasta_dir(dataDir, sampleMap, outF, sampleID=None, 
                       idPattern=None, utf16=False):
    """
    Given a sample ID and a mapping, modify a Sanger FASTA file 
    to include the barcode and 'primer' in the sequence data 
    and change the description line as needed.
    """
    for item in os.listdir(dataDir):
        if (osp.isfile(os.path.join(dataDir, item)) and 
            osp.splitext(item)[1] in file_types['fasta']):
            fh = open_enc(os.path.join(dataDir, item), utf16)
            records = SeqIO.parse(fh, 'fasta')
            for record in records:
                if idPattern is not None:
                    sep, field = idPattern
                    sampleID = record.id.split(sep)[field - 1]
                record.seq = (sampleMap[sampleID].barcode + 
                              sampleMap[sampleID].primer + 
                              record.seq)
                SeqIO.write(record, outF, 'fasta')
            fh.close()


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
    parser.add_argument('--version', action='version', 
                        version='Sanger-QIIMify 0.1')
    parser.add_argument('input_dir', 
                        help="The directory containing sequence data files. \
                              Assumes all data files are placed in this \
                              directory. For files organized within folders by\
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
    parser.add_argument('-i', '--identifier_pattern', action=ValidateIDPattern,
                        nargs=2, metavar=('SEPARATOR', 'FIELD_NUMBER'),
                        help="Indicates how to extract the Sample ID from the \
                              description line. Specify two things: 1. Field \
                              separator, 2. Field number of Sample ID (1 or \
                              greater). If the separator is a space or tab, \
                              use \s or \\t respectively. \
                              Example: >ka-SampleID-2091, use -i - 2, \
                              indicating - is the separator and the Sample ID \
                              is field #2.")
    parser.add_argument('-s', '--subdirectories',
                        help="Specify that sequence data files are contained \
                              in subfolders named by Sample ID.",
                        action='store_true')
    # data input file type options
    parser.add_argument('-q', '--qual', action='store_true', default=False,
                        help="Instruct the program to look for quality \
                              input files")
    parser.add_argument('-u', '--utf16', action='store_true', default=False,
                        help="UTF-16 encoded input files")
    
    return parser.parse_args(), parser
    

    
def main():
    args, parser = handle_program_options()
    
    if args.subdirectories:
        sampleIDs, dirs = gather_hierarchical_sample_ids(args.input_dir, 
                                                         args.utf16)
    else:
        if args.identifier_pattern:
            sampleIDs =  gather_sample_ids(args.input_dir, 
                                           args.identifier_pattern, 
                                           args.utf16)
        else:
            msg = ' '.join(['Error: Must specify identifier pattern:',
                            '-i separator field_number'])
            parser.exit(-1, msg)

    barcodes = generate_barcodes(len(sampleIDs), codeLen=args.barcode_length)
    sampleMap = write_mapping_file(args.output_map_file, sampleIDs, barcodes)

    with args.output_sequence_file:
        if args.subdirectories:
            for dataDir, sampleID in zip(dirs,sampleIDs):
                scrobble_fasta_dir(dataDir, sampleMap,
                                   args.output_sequence_file, 
                                   sampleID=sampleID, 
                                   utf16=args.utf16)
        else:
            scrobble_fasta_dir(args.input_dir, sampleMap, 
                               args.output_sequence_file, 
                               idPattern=args.identifier_pattern,
                                utf16=args.utf16)

if __name__ == '__main__':
    main()












