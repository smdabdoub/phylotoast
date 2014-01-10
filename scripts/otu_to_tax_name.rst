OTU to Tax Name
================

Convert a list of OTU IDs to a list of OTU IDs paired with Genus_species
identifiers.

usage: otu_to_tax_name.py [-h] -i OTU_ID_FP -t TAXONOMY_FP [-o OUTPUT_FP]

optional arguments:
  -h, --help            show this help message and exit
  -i OTU_ID_FP, --otu_id_fp OTU_ID_FP
                        A text file containing a list (one per line) of OTU
                        IDs.
  -t TAXONOMY_FP, --taxonomy_fp TAXONOMY_FP
                        A file associating OTU ID with a full taxonomic
                        specifier.
  -o OUTPUT_FP, --output_fp OUTPUT_FP
                        A new file containing a list OTU IDs paired with of
                        short taxonomic identifiers.