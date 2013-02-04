#!/usr/bin/env python

'''
Created on Jan 10, 2013

@author: Shareef M Dabdoub

Generate PBS scripts for submission to the Ohio Supercomputer Center to 
run the QIIME parallel blast pick OTUs script on multiple input sequence 
data sets.
'''
import argparse
import os.path as osp


PBS_SCRIPT = """#PBS -N {job_name}_{job_num}

# Walltime Limit: hh:nn:ss 
#PBS -l walltime={walltime}:00:00 

# Node Specification:
#PBS -l nodes=1:ppn=12

# Keep output
#PBS -j oe

set -x

module load qiime
cd $PBS_O_WORKDIR
cp $HOME/local/oakley/greengenes/gg_12_10_otus/rep_set/97_otus.fasta $TMPDIR
cp {job_num}.fna $TMPDIR
cd $TMPDIR

time parallel_pick_otus_blast.py -i {job_num}.fna -r 97_otus.fasta -O 12 -o bpo.{job_num}

cp -r bpo.{job_num} $PBS_O_WORKDIR
"""


def handle_program_options():
    parser = argparse.ArgumentParser(description="Generate PBS scripts for \
                                     submission to the OSC to run the QIIME \
                                     parallel blast pick OTUs script on \
                                     multiple input sequence data sets.")
    parser.add_argument('-i', '--input_fna', required=True, nargs='+',
                        help="The names of the sequence files that will be \
                              have PBS scripts generated to process them. \
                              The expected input is from the \
                              split_sequence_data.py script (e.g. 0.fna, \
                              1.fna, ..., n.fna).")
    parser.add_argument('-t', '--walltime', type=int, default=10,
                        help="The maximum running time to specify to the \
                              OSC queuing system for each script.")
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
    
    for fname in args.input_fna:
        fnum = osp.splitext(osp.split(fname)[1])[0]
        outFN = '{}.pbs'.format(fnum)
        with open(outFN, 'w') as outF:
            outF.write(PBS_SCRIPT.format(job_name=args.job_name, job_num=fnum,
                                         walltime=args.walltime))

        if args.verbose:
            print outFN


if __name__ == '__main__':
    main()