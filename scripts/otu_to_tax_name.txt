==================
otu_to_tax_name.py
==================

Convert a list of OTU IDs to a list of OTU IDs paired with Genus_species
identifiers.

    .. code-block:: bash
    
        usage: otu_to_tax_name.py [-h] -i OTU_ID_FP -t TAXONOMY_FP [-o OUTPUT_FP]


Required arguments
^^^^^^^^^^^^^^^^^^

.. cmdoption:: -i OTU_ID_FP, --otu_id_fp OTU_ID_FP

    A text file containing a list (one per line) of OTU IDs.

.. cmdoption:: -t TAXONOMY_FP, --taxonomy_fp TAXONOMY_FP

    A file associating OTU ID with a full taxonomic specifier.

Optional arguments
^^^^^^^^^^^^^^^^^^
.. cmdoption:: -o OUTPUT_FP, --output_fp OUTPUT_FP

    A new file containing a list OTU IDs paired with of short taxonomic identifiers.
    
.. cmdoption:: -h, --help
    
    Show the help message and exit