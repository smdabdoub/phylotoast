====================
condense_workflow.py
====================

This workflow script will run all three steps of the OTU condensing pipeline
automatically with the default output file settings.

    .. code-block:: bash
    
        usage: condense_workflow.py [-h] -i ASSIGNED_TAXONOMY_FN -r REP_SET_FN -s SEQS_OTUS_FN [-L {k,p,c,o,f,g,s}] [-v]

Required arguments
^^^^^^^^^^^^^^^^^^
    
.. cmdoption:: -i ASSIGNED_TAXONOMY_FN, --assigned_taxonomy_fn ASSIGNED_TAXONOMY_FN

    The taxonomy file output by the assign_taxonomy script.
    
.. cmdoption::  -r REP_SET_FN, --rep_set_fn REP_SET_FN

    The set of representative sequences.

.. cmdoption::  -s SEQS_OTUS_FN, --seqs_otus_fn SEQS_OTUS_FN

    The list of OTU IDs and their associated sequence IDs.

Optional arguments
^^^^^^^^^^^^^^^^^^
    
.. cmdoption::  -L {k,p,c,o,f,g,s}, --phylogenetic_level {k,p,c,o,f,g,s}

    Set the phylogenetic level at which to define OTUs for
    condensing and downstream processing. Defaults to species level.

.. cmdoption:: -h, --help
    
    Show the help message and exit   
    
.. cmdoption::  -v, --verbose

    Print detailed information about script operation.