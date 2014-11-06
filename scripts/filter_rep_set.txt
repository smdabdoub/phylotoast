=================
filter_rep_set.py
=================

Step 2 of the condensing process. Filter the representative sequence set to
include only those sequences that map to unique OTUs.

    .. code-block:: bash
    
        usage: filter_rep_set.py [-h] -r REP_SET_FN -u UNIQUE_OTUS_FN [-o OUTPUT_FILTERED_REP_SET_FN] [-v]

Required arguments
^^^^^^^^^^^^^^^^^^

.. cmdoption:: -r REP_SET_FN, --rep_set_fn REP_SET_FN

    The set of representative sequences.

.. cmdoption:: -u UNIQUE_OTUS_FN, --unique_otus_fn UNIQUE_OTUS_FN

    The condensed assigned taxonomy file.

Optional arguments
^^^^^^^^^^^^^^^^^^
    
.. cmdoption:: -o OUTPUT_FILTERED_REP_SET_FN, --output_filtered_rep_set_fn OUTPUT_FILTERED_REP_SET_FN

    The filtered representative set. By default outputs to condensed_rep_set.fna

.. cmdoption:: -h, --help
    
    Show the help message and exit
    
.. cmdoption::  -v, --verbose

    Print detailed information about script operation.