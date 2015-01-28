#!/usr/bin/env python

'''
Created on Sep 19, 2012

Author: Shareef Dabdoub

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
import random
import re
import sys
try:
    from Bio import SeqIO
except ImportError as ie:
    sys.exit('Import Error. Please install missing module: {}'.format(ie))
from Bio import SeqIO


file_types = {'fasta': ('.fasta', '.fa', '.txt', '.seq', '.Seq'),
              'fastq': ('.fastq'),
              'abi': ('.abi', '.ab1'),
              'qual': ('.qual')}

MapRecord = namedtuple('MapRecord', 'barcode primer treatment descr')
IDPattern = namedtuple('IDPattern', 'separator field')
TreatmentType = namedtuple('TreatmentType', 'name value')


class ValidateIDPattern(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
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

    :param fn: Name of the file to open
    :param utf16: Boolean flag to specify a UTF-16 file
    """
    if utf16:
        return codecs.open(fn, encoding='utf-16')  # -be')
    return open(fn, 'rU')


def gather_sample_ids(sampleDir, idopt, utf16):
    """

    """
    sampleIDs = []
    useFname = False
    if isinstance(idopt, tuple):
        separator, field = idopt
    else:
        useFname = True
    for item in os.listdir(sampleDir):
        iPath = os.path.join(sampleDir, item)
        if (osp.isfile(iPath) and
            osp.splitext(item)[1] in file_types['fasta']):
            try:
                with open_enc(iPath, utf16) as fh:
                    record = SeqIO.parse(fh, 'fasta').next()
                    if useFname:
                        sid = osp.splitext(item)[0]
                    else:
                        sid = record.id.split(separator)[field-1]
                    sampleIDs.append(sid)
            except StopIteration:
                print 'Invalid FASTA file: %s' % iPath

    return sampleIDs


def generate_barcodes(nIds, codeLen=12):
    """
    Given a list of sample IDs generate unique n-base barcodes for each.
    Note that only 4^n unique barcodes are possible.
    """
    def next_code(b, c, i):
        return c[:i] + b + (c[i+1:] if i < -1 else '')

    def rand_base():
        return random.choice(['A', 'T', 'C', 'G'])

    def rand_seq(n):
        return ''.join([rand_base() for _ in range(n)])

    # homopolymer filter regex: match if 4 identical bases in a row
    hpf = re.compile('aaaa|cccc|gggg|tttt', re.IGNORECASE)

    while True:
        codes = [rand_seq(codeLen)]
        if (hpf.search(codes[0]) is None):
            break
    idx = 0

    while len(codes) < nIds:
        idx -= 1
        if idx < -codeLen:
            idx = -1
            codes.append(rand_seq(codeLen))
        else:
            nc = next_code(rand_base(), codes[-1], idx)
            if hpf.search(nc) is None:
                codes.append(nc)
        codes = list(set(codes))

    return codes


def write_mapping_file(mapF, sampleIDs, barcodes, treatment=None):
    """
    Given a mapping from sample IDs to barcodes/primers/other info,
    write out a QIIME-compatible mapping file.

    File format described at: http://qiime.org/documentation/file_formats.html
    Note that primer column can be left blank.
    """
    header = ['#SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
              'Description', '\n']
    if treatment is not None:
        header.insert(-2, treatment.name)
    primer = 'AAACCCGGGTTTAAACTGAC'
    sampleMap = {}

    with mapF:
        mapF.write('\t'.join(header))
        for sid, bc in zip(sampleIDs, barcodes):
            sampleMap[sid] = MapRecord(bc, primer, treatment.value, sid)
            line = [sid, bc, primer, sid]
            if treatment is not None:
                line.insert(-1, treatment.value)
            mapF.write('\t'.join(line) + '\n')

    return sampleMap


def scrobble_data_dir(dataDir, sampleMap, outF, qualF=None, idopt=None,
                      utf16=False):
    """
    Given a sample ID and a mapping, modify a Sanger FASTA file
    to include the barcode and 'primer' in the sequence data
    and change the description line as needed.
    """
    seqcount = 0
    outfiles = [osp.split(outF.name)[1]]
    if qualF:
        outfiles.append(osp.split(qualF.name)[1])

    for item in os.listdir(dataDir):
        if item in outfiles or not osp.isfile(os.path.join(dataDir, item)):
            continue
        # FASTA files
        if osp.splitext(item)[1] in file_types['fasta']:
            fh = open_enc(os.path.join(dataDir, item), utf16)
            records = SeqIO.parse(fh, 'fasta')
            for record in records:
                if isinstance(idopt, tuple):
                    sep, field = idopt
                    sampleID = record.id.split(sep)[field - 1]
                else:
                    sampleID = osp.splitext(item)[0]
                record.seq = (sampleMap[sampleID].barcode +
                              sampleMap[sampleID].primer +
                              record.seq)
                SeqIO.write(record, outF, 'fasta')
                seqcount += 1
            fh.close()
        # QUAL files
        elif qualF and osp.splitext(item)[1] in file_types['qual']:
            fh = open_enc(os.path.join(dataDir, item), utf16)
            records = SeqIO.parse(fh, 'qual')
            for record in records:
                mi = sampleMap[sampleMap.keys()[0]]
                quals = [40 for _ in range(len(mi.barcode) + len(mi.primer))]
                record.letter_annotations['phred_quality'][0:0] = quals
                SeqIO.write(record, qualF, 'qual')
            fh.close()
    return seqcount


def handle_program_options():
    """
    Uses the built-in argparse module to handle command-line options for the
    program.

    :return: The gathered command-line options specified by the user
    :rtype: argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(description="Convert Sanger-sequencing \
                                     derived data files for use with the \
                                     metagenomics analysis program QIIME, by \
                                     extracting Sample ID information, adding\
                                     barcodes and primers to the sequence \
                                     data, and outputting a mapping file and\
                                     single FASTA-formatted sequence file \
                                     formed by concatenating all input data.")
    parser.add_argument('-i', '--input_dir', required=True,
                        help="The directory containing sequence data files. \
                              Assumes all data files are placed in this \
                              directory. For files organized within folders by\
                              sample, use -s in addition.")
    parser.add_argument('-m', '--map_file', default='map.txt',
                        help="QIIME-formatted mapping file linking Sample IDs \
                              with barcodes and primers.")
    parser.add_argument('-o', '--output', default='output.fasta',
                        metavar='OUTPUT_FILE',
                        help="Single file containing all sequence data found \
                              in input_dir, FASTA-formatted with barcode and \
                              primer preprended to sequence. If the -q option \
                              is passed, any quality data will also be output \
                              to a single file of the same name with a .qual \
                              extension.")
    parser.add_argument('-b', '--barcode_length', type=int, default=12,
                        help="Length of the generated barcode sequences. \
                              Default is 12 (QIIME default), minimum is 8.")

    parser.add_argument('-q', '--qual', action='store_true', default=False,
                        help="Instruct the program to look for quality \
                              input files")
    parser.add_argument('-u', '--utf16', action='store_true', default=False,
                        help="UTF-16 encoded input files")

    parser.add_argument('-t', '--treatment',
                        help="Inserts an additional column into the mapping \
                              file specifying some treatment or other variable\
                              that separates the current set of sequences \
                              from any other set of seqeunces. For example:\
                              -t DiseaseState=healthy")

    # data input options
    sidGroup = parser.add_mutually_exclusive_group(required=True)
    sidGroup.add_argument('-d', '--identifier_pattern',
                          action=ValidateIDPattern,
                          nargs=2, metavar=('SEPARATOR', 'FIELD_NUMBER'),
                          help="Indicates how to extract the Sample ID from \
                               the description line. Specify two things: \
                               1. Field separator, 2. Field number of Sample \
                               ID (1 or greater). If the separator is a space \
                               or tab, use \s or \\t respectively. \
                               Example: >ka-SampleID-2091, use -i - 2, \
                               indicating - is the separator and the Sample ID\
                               is field #2.")
    sidGroup.add_argument('-f', '--filename_sample_id', action='store_true',
                          default=False, help='Specify that the program should\
                          the name of each fasta file as the Sample ID for use\
                          in the mapping file. This is meant to be used when \
                          all sequence data for a sample is stored in a single\
                          file.')

    return parser.parse_args()


def main():
    args = handle_program_options()
    qualOutFN = ''

    try:
        with open(args.input_dir):
            pass
    except IOError as ioe:
        sys.exit('\nError with input directory:{}\n'.format(ioe))

    try:
        with open(args.map_file):
            pass
    except IOError as ioe:
        sys.exit('\nError with QIIME-formatted mapping file:{}\n'.format(ioe))

    if args.treatment is not None:
        if '=' in args.treatment and args.treatment > 2:
            name, value = args.treatment.split('=')
            args.treatment = TreatmentType(name, value)
        else:
            msg = ('Treatment column specification must be of the form: ' +
                   'Treatment=type')
            print msg

    sampleIDs = gather_sample_ids(args.input_dir,
                                  args.identifier_pattern or
                                  args.filename_sample_id,
                                  args.utf16)
    barcodes = generate_barcodes(len(sampleIDs), codeLen=args.barcode_length)

    mode = 'a' if os.path.exists(args.map_file) else 'w'
    with open(args.map_file, mode) as mapfile:
        sampleMap = write_mapping_file(mapfile, sampleIDs, barcodes,
                                       args.treatment)

    if args.qual:
        qualOutFN = osp.splitext(osp.split(args.output)[1])[0] + '.qual'
        with open(args.output, 'w') as outf, open(qualOutFN, 'w') as qoutf:
            seqcount = scrobble_data_dir(args.input_dir, sampleMap, outf,
                                         idopt=args.identifier_pattern or
                                         args.filename_sample_id,
                                         qualF=qoutf,
                                         utf16=args.utf16)
    else:
        with open(args.output, 'w') as outf:
            seqcount = scrobble_data_dir(args.input_dir, sampleMap, outf,
                                         idopt=args.identifier_pattern or
                                         args.filename_sample_id,
                                         utf16=args.utf16)

    print
    print 'Processing completed: %i samples, %i sequences' % (len(sampleMap),
                                                              seqcount)

if __name__ == '__main__':
    main()
