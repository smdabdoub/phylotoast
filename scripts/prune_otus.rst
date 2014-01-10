Prune OTUs
===========

Parse the OTU-sequence data in two steps. First remove any OTUs that occur in
less than a user-defined percent of samples (default 5%). Second, remove any
OTUs that make up less than a user-defined percentage of the overall sequences
(default 0.01%)

usage: prune_otus.py [-h] -i SEQS_OTUS_FN -t ID_TO_TAXONOMY_FN
                     [-p PERCENT_OF_SAMPLES] [-s PERCENT_OF_SEQUENCES]
                     [-l {k,p,c,o,f,g,s}] [-o OUTPUT_PRUNED_OTUS_FN]
                     [--output_removed_otus_fn OUTPUT_REMOVED_OTUS_FN] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i SEQS_OTUS_FN, --seqs_otus_fn SEQS_OTUS_FN
                        The output from the pick OTUs step, e.g. seqs_otus.txt
  -t ID_TO_TAXONOMY_FN, --id_to_taxonomy_fn ID_TO_TAXONOMY_FN
                        Path to tab-delimited file mapping sequences to
                        assigned taxonomy.
  -p PERCENT_OF_SAMPLES, --percent_of_samples PERCENT_OF_SAMPLES
                        OTUs that occur in less than this percent of samples
                        will be removed. Default is 5 percent.
  -s PERCENT_OF_SEQUENCES, --percent_of_sequences PERCENT_OF_SEQUENCES
                        OTUs that occur in less than this percent of total
                        sequences will be removed. Default is 0.01 percent.
  -l {k,p,c,o,f,g,s}, --phylogenetic_level {k,p,c,o,f,g,s}
                        Set the phylogenetic level at which to join OTUs for
                        consideration in pruning. Default is 'g'(group).
  -o OUTPUT_PRUNED_OTUS_FN, --output_pruned_otus_fn OUTPUT_PRUNED_OTUS_FN
                        The main output file that will contain the remaining
                        OTUs and sequence IDs.
  --output_removed_otus_fn OUTPUT_REMOVED_OTUS_FN
                        The file to write out the set of OTUs that were
                        removed by the filter.
  -v, --verbose