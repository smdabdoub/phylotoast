#!/usr/bin/env python
'''
Created on Jan 31, 2013

@author: Shareef Dabdoub

From an input FASTA file, filter out all sequences with barcodes matching those 
in an input mapping file and write to a new file.
'''
import argparse
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from qiime_tools import util


def gather_sequences(fastaFN, mapFN):
    barcodes = util.parse_map_file(mapFN, 1).keys()
    seqs = []
    bcodelen = len(barcodes[0])
    count = 0
    
    for record in SeqIO.parse(fastaFN, "fasta", generic_dna):
        count += 1
        if str(record.seq)[:bcodelen] in barcodes:
            seqs.append(record)
            
    return seqs, count

def gather_quality_data(qualityFN, seqs):
    quals = []
    seqids = [seq.id for seq in seqs]
    count = 0
    
    inQualRecs = SeqIO.to_dict(SeqIO.parse(qualityFN, "qual"))
    for recID in inQualRecs:
        count += 1
        if not count % 1000:
            print count
        if recID in seqids:
            quals.append(inQualRecs[recID])

    return quals

def handle_program_options():
    parser = argparse.ArgumentParser(description="From an input FASTA file, \
                                     filter all sequences with barcodes \
                                     matching those in an input mapping file.")
    parser.add_argument('-i','--input_fasta_fn', required=True, 
                        help="The sequence data file to be filtered.")
    parser.add_argument('-m', '--mapping_fn', required=True,
                        help="The mapping file containing the barcodes you \
                              want filtered sequenced to contain.")
    parser.add_argument('-q', '--quality_fn',
                        help="The quality data file. If you plan to use quality\
                              data with split_libraries.py, you have to filter \
                              the quality data as well.")
    parser.add_argument('-o', '--output_prefix', default='filtered',
                        help="The prefix for the output filtered data")
    parser.add_argument('-v', '--verbose', action='store_true')
    
    return parser.parse_args()
        

    
def main():
    args = handle_program_options()

    seqs, count = gather_sequences(args.input_fasta_fn, args.mapping_fn)
    SeqIO.write(seqs, args.output_prefix+'.fasta', 'fasta')

    if args.quality_fn:
        quals = gather_quality_data(args.quality_fn, seqs)
        SeqIO.write(quals, args.output_prefix+'.qual', 'qual')
    
    if args.verbose:
        print '%i sequences in input file.' % count
        print '%i sequences filtered to output file.' % len(seqs)
    
    
if __name__ == '__main__':
    main()