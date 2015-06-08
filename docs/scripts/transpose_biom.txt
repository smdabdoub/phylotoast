=================
transpose_biom.py
=================

Transpose a BIOM-format file so that the matrix is sample by species.

    .. code-block:: bash
    
        usage: transpose_biom.py [-h] -i INPUT_BIOM_FP -m MAPPING [-c MAP_CATEGORY] -o OUTPUT_BIOM_FP [-v]

Required arguments
^^^^^^^^^^^^^^^^^^

.. cmdoption:: -i INPUT_BIOM_FP, --input_biom_fp INPUT_BIOM_FP

    The BIOM-format file.

.. cmdoption:: -m MAPPING, --mapping MAPPING

    The mapping file specifying group information for each sample.

.. cmdoption:: -o OUTPUT_BIOM_FP, --output_biom_fp OUTPUT_BIOM_FP

    The BIOM-format file to write.

Optional arguments
^^^^^^^^^^^^^^^^^^

.. cmdoption:: -c MAP_CATEGORY, --map_category MAP_CATEGORY

    A mapping category, such as TreatmentType, that will
    be used to split the data into separate BIOM files;
    one for each value found in the category.
    
.. cmdoption:: -h, --help
    
    Show the help message and exit    
    
.. cmdoption:: -v, --verbose

    Print detailed information about script operation.