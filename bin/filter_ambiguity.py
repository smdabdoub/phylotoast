'''
Created on Dec 6, 2012

Author: Shareef M. Dabdoub
'''
import sys
try:
    from Bio import SeqIO
except ImportError as ie:
    sys.exit('Import Error. Please install missing module: {}'.format(ie))
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


def filter_ambiguity(records, percent=0.5):  # , repeats=6)
    """
    Filters out sequences with too much ambiguity as defined by the method
    parameters.

    :type records: list
    :param records: A list of sequences
    :type repeats: int
    :param repeats: Defines the number of repeated N that trigger truncating a
                    sequence.
    :type percent: float
    :param percent: Defines the overall percentage of N in a sequence that
                     will cause the sequence to be filtered out.
    """
    seqs = []
    # Ns = ''.join(['N' for _ in range(repeats)])
    count = 0
    for record in records:
        if record.seq.count('N')/float(len(record)) < percent:
#            pos = record.seq.find(Ns)
#            if pos >= 0:
#                record.seq = Seq(str(record.seq)[:pos])
            seqs.append(record)
        count += 1

    return seqs, count


def handle_program_options():
    """
    Uses the built-in argparse module to handle command-line options for the
    program.

    :return: The gathered command-line options specified by the user
    :rtype: argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(description="Filter an input \
                                                  FASTA-formatted file to \
                                                  remove or truncate sequences\
                                                  based on ambiguous base (N) \
                                                  content.")
    parser.add_argument('--version', action='version',
                        version='Sequence Ambiguity Filter v0.1')
    parser.add_argument('fasta_input', help="QIIME-formatted mapping file \
                                             linking Sample IDs with barcodes \
                                             and primers.")
    parser.add_argument('-o', '--output', default='output.fna',
                        help='The name of the file to output the set of \
                              filtered sequences. Default: \'output.fna\'.')
#    parser.add_argument('-r', '--repeats', type=int, default=6,
#                        help='Truncates a sequence when a string of ambiguous \
#                              bases (N) of REPEATS or longer is found. \
#                              Default: REPEATS=6.')
    parser.add_argument('-p', '--percent', type=int, default=5,
                        help='Removes any sequence containing the specified \
                              percentage (or greater) of ambiguous bases (N).\
                              Default: PERCENT=5')
    parser.add_argument('-v', '--verbose', action='store_true')

    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.fasta_input):
            pass
    except IOError as ioe:
        sys.exit('\nError with QIIME formatted mapping file:{}\n'.format(ioe))

    with open(args.fasta_input, 'rU') as inF:
        in_records = SeqIO.parse(inF, 'fasta')
        records, count = filter_ambiguity(in_records, args.percent/100.0)  # , args.repeats)

    SeqIO.write(records, args.output, "fasta")

    if args.verbose:
        print '%i sequences found.' % count
        print '%i sequences kept.' % len(records)
        print
        print 'Output written to: %s' % args.output

if __name__ == '__main__':
    main()
