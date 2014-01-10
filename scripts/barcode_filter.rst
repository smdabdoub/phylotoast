Barcode Filter
===============

From an input FASTA file, filter all sequences with barcodes matching those in
an input mapping file.

usage: barcode_filter.py [-h] -i INPUT_FASTA_FN -m MAPPING_FN [-q QUALITY_FN]
                         [-o OUTPUT_PREFIX] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FASTA_FN, --input_fasta_fn INPUT_FASTA_FN
                        The sequence data file to be filtered.
  -m MAPPING_FN, --mapping_fn MAPPING_FN
                        The mapping file containing the barcodes you want
                        filtered sequenced to contain.
  -q QUALITY_FN, --quality_fn QUALITY_FN
                        The quality data file. If you plan to use quality data
                        with split_libraries.py, you have to filter the
                        quality data as well.
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        The prefix for the output filtered data
  -v, --verbose