==================================
assign_taxonomy_by_blast _result.py
==================================

Assign taxonomy to a rep set of OTUs that were chosen by BLAST from an
annotated database.

    .. code-block:: bash
    
        usage: assign_taxonomy_by_blast_result.py [-h] -i REP_SET_FP -t ID_TO_TAXONOMY_FP [-o ASSIGNED_TAXONOMY_FP] [-v]

Required arguments
^^^^^^^^^^^^^^^^^^

.. cmdoption:: -i REP_SET_FP, --rep_set_fp REP_SET_FP

    The set of representative sequences.
    
.. cmdoption:: -t ID_TO_TAXONOMY_FP, --id_to_taxonomy_fp ID_TO_TAXONOMY_FP

    Path to tab-delimited file mapping sequences to assigned taxonomy.
    
Optional arguments
^^^^^^^^^^^^^^^^^^

.. cmdoption:: -o ASSIGNED_TAXONOMY_FP, --assigned_taxonomy_fp ASSIGNED_TAXONOMY_FP

    The path to the result file. By default outputs to assigned_taxonomy.txt

.. cmdoption:: -h, --help
    
    Show the help message and exit    
    
.. cmdoption::  -v, --verbose

    Print detailed information about script operation.