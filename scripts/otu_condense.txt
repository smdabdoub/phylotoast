===============
otu_condense.py
===============

Step 1 of the condensing process. Take a taxonomy table from the
assign_taxonomy QIIME script and prune all redundant taxonomy strings

    .. code-block:: bash
    
        usage: otu_condense.py [-h] -i INPUT_ASSIGNED_TAXONOMY [-p PRUNED_OUTPUT_FILE] [-n NON_UNIQUE_OUTPUT_FILE] [-l {k,p,c,o,f,g,s}] [-v]

Required arguments
^^^^^^^^^^^^^^^^^^

.. cmdoption:: -i INPUT_ASSIGNED_TAXONOMY, --input_assigned_taxonomy INPUT_ASSIGNED_TAXONOMY

    The taxonomy file output by the assign_taxonomy script.

Optional arguments
^^^^^^^^^^^^^^^^^^

.. cmdoption:: -p PRUNED_OUTPUT_FILE, --pruned_output_file PRUNED_OUTPUT_FILE

    The output file for the pruned taxonomy list. Defaults
    to condensed_assigned_taxonomy.txt

.. cmdoption:: -n NON_UNIQUE_OUTPUT_FILE, --non_unique_output_file NON_UNIQUE_OUTPUT_FILE

    The file will contain a list of pruned OTU IDs
    associated with the OTU IDs they replaced. Defaults to
    nonunique_otu_matrix.txt

.. cmdoption:: -l {k,p,c,o,f,g,s}, --phylogenetic_level {k,p,c,o,f,g,s}

    Set the phylogenetic level at which to define OTUs for
    condensing and downstream processing. Defaults to species level.

.. cmdoption:: -h, --help
    
    Show the help message and exit    

.. cmdoption::  -v, --verbose

    Print detailed information about script operation.