Merge OTU Results
===================

Distributing sequence data across the cluster for OTU picking results in a set
of result files that need to be merged into a single pick otus result.

usage: merge_otu_results.py [-h] [-o OUTPUT_FN] [-v]
                            pick_otus_results [pick_otus_results ...]

positional arguments:
  pick_otus_results     The result files from multiple runs of a pick otus
                        script that need to be merged.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_FN, --output_fn OUTPUT_FN
                        The name of the file the merged results will be
                        written to.
  -v, --verbose
