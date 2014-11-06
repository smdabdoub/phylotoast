=================
barcode_filter.py
=================

From an input FASTA file, filter all sequences with barcodes matching those in
an input mapping file.

    .. code-block:: bash
    
        usage: barcode_filter.py [-h] -i INPUT_FASTA_FN -m MAPPING_FN [-q QUALITY_FN] [-o OUTPUT_PREFIX] [-v]

Required arguments
^^^^^^^^^^^^^^^^^^
    
.. cmdoption:: -i INPUT_FASTA_FN, --input_fasta_fn INPUT_FASTA_FN

    The sequence data file to be filtered.

.. cmdoption:: -m MAPPING_FN, --mapping_fn MAPPING_FN

    The mapping file containing the barcodes you want filtered sequenced to contain.

Optional arguments
^^^^^^^^^^^^^^^^^^
    
.. cmdoption:: -q QUALITY_FN, --quality_fn QUALITY_FN

    The quality data file. If you plan to use quality data with split_libraries.py, 
    you have to filter the quality data as well.
    
.. cmdoption:: -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX

    The prefix for the output filtered data

.. cmdoption:: -h, --help
    
    Show the help message and exit    
    
.. cmdoption::  -v, --verbose

    Print detailed information about script operation.