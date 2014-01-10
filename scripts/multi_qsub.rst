Multi Qsub
============

Submit multiple PBS job scripts to the queuing system (qsub) and store the
output job IDs.

usage: multi_qsub.py [-h] [-t] job_scripts [job_scripts ...]

positional arguments:
  job_scripts  The job script files to submit to the queuing system.

optional arguments:
  -h, --help   show this help message and exit
  -t, --test   Only print each of the qsub commands instead of actually
               running the commands.