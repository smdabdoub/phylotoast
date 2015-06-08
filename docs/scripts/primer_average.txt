=================
primer_average.py
=================

Combine multi-primer pick OTUs results files into a single results file while
at the same time averaging sequence counts per sample for OTUs shared between
the primer-set results. See reference: Kumar PS et al. (2011)
doi:10.1371/journal.pone.0020956

    .. code-block:: bash
    
        usage: primer_average.py [-h] --p1 P1 --p2 P2 [-o OUTPUT_FP] [-v]

Required arguments
^^^^^^^^^^^^^^^^^^

.. cmdoption:: --p1 P1

    Primer-set 1 seqs_otus results files.

.. cmdoption:: --p2 P2

    Primer-set 2 seqs_otus results files.

Optional arguments
^^^^^^^^^^^^^^^^^^
  
.. cmdoption:: -o OUTPUT_FP, --output_fp OUTPUT_FP

    The combined seqs_otus file that has been averaged by
    shared OTU entries. Default: combined_seqs_otus.txt

.. cmdoption:: -h, --help
    
    Show the help message and exit

.. cmdoption:: -v, --verbose

    Print detailed information about script operation.