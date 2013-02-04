#!/usr/bin/env python

'''
Created on Feb 4, 2013

@author: Shareef Dabdoub

Submit multiple PBS job scripts to the queuing system (qsub) and store 
the output job IDs
'''
import argparse
import subprocess

def handle_program_options():
    parser = argparse.ArgumentParser(description="Submit multiple PBS job \
                                     scripts to the queuing system (qsub) and\
                                     store the output job IDs.")
    parser.add_argument('job_scripts', nargs='+',
                        help="The job script files to submit to the queuing\
                              system.")
    parser.add_argument('-t', '--test', action='store_true', 
                        help="Only print each of the qsub commands instead of\
                              actually running the commands.")
    
    return parser.parse_args()


def main():
    args = handle_program_options()
    
    for script in args.job_scripts:
        if args.test:
            print 'qsub {}'.format(script)
        else:
            print subprocess.check_output(["qsub", script])

if __name__ == '__main__':
    main()