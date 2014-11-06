=====================
pick_otus_condense.py
=====================

Step 3 of the condensing process. Condense the QIIME pick_otus.py script
output by moving the sequences associated with non-unique OTUs to OTU IDs that
were identified as unique.

    .. code-block:: bash
    
        usage: pick_otus_condense.py [-h] -s SEQS_OTUS -n NON_UNIQUE_OTU_MATRIX [-o CONDENSED_SEQS_OTUS_FILE] [-v]

Required arguments
^^^^^^^^^^^^^^^^^^

.. cmdoption:: -s SEQS_OTUS, --seqs_otus SEQS_OTUS

    The list of OTU IDs and their associated sequence IDs.

.. cmdoption:: -n NON_UNIQUE_OTU_MATRIX, --non_unique_otu_matrix NON_UNIQUE_OTU_MATRIX

    The list of unique OTU IDs associated with the OTU IDs they replaced.

.. cmdoption:: -o CONDENSED_SEQS_OTUS_FILE, --condensed_seqs_otus_file CONDENSED_SEQS_OTUS_FILE

    The condensed set of OTU IDs and the matching sequences. By default outputs to condensed_seqs_otus.txt

Optional arguments
^^^^^^^^^^^^^^^^^^

.. cmdoption:: -h, --help
    
    Show the help message and exit

.. cmdoption:: -v, --verbose

    Print detailed information about script operation.