Condense Workflow
===================

This workflow script will run all three steps of the OTU condensing pipeline
automatically with the default output file settings.

usage: condense_workflow.py [-h] -i ASSIGNED_TAXONOMY_FN -r REP_SET_FN -s
                            SEQS_OTUS_FN [-L {k,p,c,o,f,g,s}] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i ASSIGNED_TAXONOMY_FN, --assigned_taxonomy_fn ASSIGNED_TAXONOMY_FN
                        The taxonomy file output by the assign_taxonomy
                        script.
  -r REP_SET_FN, --rep_set_fn REP_SET_FN
                        The set of representative sequences.
  -s SEQS_OTUS_FN, --seqs_otus_fn SEQS_OTUS_FN
                        The list of OTU IDs and their associated sequence IDs.
  -L {k,p,c,o,f,g,s}, --phylogenetic_level {k,p,c,o,f,g,s}
                        Set the phylogenetic level at which to define OTUs for
                        condensing and downstream processing. Defaults to
                        species level.
  -v, --verbose