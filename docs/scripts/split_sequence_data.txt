======================
split_sequence_data.py
======================

Split an input FASTA-formatted sequence file into a user-specified number of
smaller files such that the sequence data is evenly distributed among them.

    .. code-block:: bash
    
        usage: split_sequence_data.py [-h] -i INPUT_FASTA_FN [-n NUM_OUTPUT_FILES] [-o OUTPUT_DIR] [-v]

Required arguments
^^^^^^^^^^^^^^^^^^

.. cmdoption:: -i INPUT_FASTA_FN, --input_fasta_fn INPUT_FASTA_FN

    The sequence data file to be split up into a series of smaller files.
    
.. cmdoption:: -n NUM_OUTPUT_FILES, --num_output_files NUM_OUTPUT_FILES

    The number of files the input data should be split into.

Optional arguments
^^^^^^^^^^^^^^^^^^

.. cmdoption:: -o OUTPUT_DIR, --output_dir OUTPUT_DIR

    The location to write the split data files.

.. cmdoption:: -h, --help
    
    Show the help message and exit    
    
.. cmdoption:: -v, --verbose

    Print detailed information about script operation.