#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Jan 10, 2013

Author: Shareef M Dabdoub

Generate cluster-computing job scripts for submission in order to run multiple
simultaneous runs of the QIIME parallel BLAST pick OTUs script.
"""

import argparse
import sys
import os.path as osp

PBS_JOB_NAME_SIZE = 15


def handle_program_options():
    parser = argparse.ArgumentParser(description="Generate cluster-computing \
                                    job scripts for submission in order to run \
                                    multiple simultaneous runs of the QIIME \
                                    parallel BLAST pick OTUs script.")
    parser.add_argument('-i', '--input_fna', required=True, nargs='+',
                        help="The names of the sequence files that will be \
                              have job scripts generated to process them. \
                              The expected input is from the \
                              split_sequence_data.py script (e.g. 0.fna, \
                              1.fna, ..., n.fna).")
    parser.add_argument('-s', '--similarity', default=0.97, type=float,
                        help="Sequence similarity threshold [default: 0.97]")
    parser.add_argument('-j', '--job_script_template', required=True,
                        help="A file template containing placeholders for \
                              variables that this script will fill in \
                              when creating a new job script for each input \
                              FASTA query file. An example file for PBS \
                              systems is included with qiime-tools.")
    parser.add_argument('-d', '--database', required=True,
                        help="The path to the sequence database file to run \
                              the BLAST against.")
    parser.add_argument('-t', '--walltime', type=int, default=10,
                        help="The maximum running time in hours for each \
                              script. Default is 10 hours.")
    parser.add_argument('-n', '--job_name', default='pick_otus',
                        help="A descriptive name for the job script that will\
                              appear when checking the job status. Max length\
                              is 15 characters, but '_#' will be appended to \
                              the name you provide to differentiate among all\
                              the jobs, so this parameter will be truncated if\
                              necessary to accommodate for the number of input\
                              files.")
    parser.add_argument('-v', '--verbose', action='store_true',
                        help="This will cause the program to print the full\
                        path for each output file to the command line. This \
                        can be used for informational purposes or to pipe (|)\
                        to the PBS multi-submission script to automate job\
                        submission as soon as the scripts are created.")

    return parser.parse_args()


def main():
    args = handle_program_options()

    try:
        with open(args.job_script_template):
            pass
    except IOError as ioe:
        sys.exit(
            '\nError with template file:{}\n'
            .format(ioe)
        )

    try:
        with open(args.database):
            pass
    except IOError as ioe:
        sys.exit(
            '\nError with sequence database file:{}\n'
            .format(ioe)
        )

    with open(args.job_script_template, 'rU') as tF:
        template = tF.read()

    for fname in args.input_fna:
        fnum = osp.splitext(osp.split(fname)[1])[0]
        job_id_len = len(str(len(args.input_fna))) + 1
        # for 0-based file names, this is one character too conservative
        args.job_name = args.job_name[:PBS_JOB_NAME_SIZE - job_id_len]
        outFN = '{}.pbs'.format(fnum)
        with open(outFN, 'w') as out:
            out.write(template.format(job_name=args.job_name, job_num=fnum,
                                      database_path=args.database,
                                      database_fname=osp.basename(args.database),
                                      similarity=args.similarity,
                                      walltime=args.walltime))

        if args.verbose:
            print outFN


if __name__ == '__main__':
    main()
