Filter Rep Set
================

Step 2 of the condensing process. Filter the representative sequence set to
include only those sequences that map to unique OTUs.

usage: filter_rep_set.py [-h] -r REP_SET_FN -u UNIQUE_OTUS_FN
                         [-o OUTPUT_FILTERED_REP_SET_FN] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -r REP_SET_FN, --rep_set_fn REP_SET_FN
                        The set of representative sequences.
  -u UNIQUE_OTUS_FN, --unique_otus_fn UNIQUE_OTUS_FN
                        The condensed assigned taxonomy file.
  -o OUTPUT_FILTERED_REP_SET_FN, --output_filtered_rep_set_fn OUTPUT_FILTERED_REP_SET_FN
                        The filtered representative set. By default outputs to
                        condensed_rep_set.fna
  -v, --verbose