Split Sequence Data
=====================

Split an input FASTA-formatted sequence file into a user-specified number of
smaller files such that the sequence data is evenly distributed among them.

usage: split_sequence_data.py [-h] -i INPUT_FASTA_FN [-n NUM_OUTPUT_FILES]
                              [-o OUTPUT_DIR] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FASTA_FN, --input_fasta_fn INPUT_FASTA_FN
                        The sequence data file to be split up.
  -n NUM_OUTPUT_FILES, --num_output_files NUM_OUTPUT_FILES
                        The number of files the input data should be split
                        into.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        The location to write the split data files.
  -v, --verbose