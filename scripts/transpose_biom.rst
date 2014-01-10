Transpose BIOM
===============

Transpose a BIOM-format file so that the matrix is sample by species.

usage: transpose_biom.py [-h] -i INPUT_BIOM_FP -m MAPPING [-c MAP_CATEGORY] -o
                         OUTPUT_BIOM_FP [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_BIOM_FP, --input_biom_fp INPUT_BIOM_FP
                        The BIOM-format file.
  -m MAPPING, --mapping MAPPING
                        The mapping file specifying group information for each
                        sample.
  -c MAP_CATEGORY, --map_category MAP_CATEGORY
                        A mapping category, such as TreatmentType, that will
                        be used to split the data into separate BIOM files;
                        one for each value found in the category.
  -o OUTPUT_BIOM_FP, --output_biom_fp OUTPUT_BIOM_FP
                        The BIOM-format file to write.
  -v, --verbose
