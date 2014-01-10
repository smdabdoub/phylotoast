iTol
=====

Create files appropriate for use in the iTol visualization program by using
the abundance data from a biom-format file and groups specified in a QIIME
mapping file. The program also modifies a Newick-format phylogenetic tree file
to use proper taxonomic names instead of OTU IDs for useful display in iTol.

usage: iTol.py [-h] -i OTU_TABLE -m MAPPING [-t INPUT_TREE] [-e OUTPUT_TRE]
               [-o OUTPUT_ITOL_TABLE] [-c MAP_CATEGORIES] [-a {MRA,NMRA,raw}]

optional arguments:
  -h, --help            show this help message and exit
  -i OTU_TABLE, --otu_table OTU_TABLE
                        The biom-format file with OTU-Sample abundance data.
  -m MAPPING, --mapping MAPPING
                        The mapping file specifying group information for each
                        sample.
  -t INPUT_TREE, --input_tree INPUT_TREE
                        A phylogenetic tree in Newick format to be modified by
                        exchanging the OTU ID node names for taxonomic names.
  -e OUTPUT_TRE, --output_tre OUTPUT_TRE
                        The output .tre file
  -o OUTPUT_ITOL_TABLE, --output_itol_table OUTPUT_ITOL_TABLE
                        Other than a phylogenetic tree, the main input to iTol
                        is a dataset file containing some representation of
                        the abundance of every OTU across the specified data
                        groups. This program provides multiple calculation
                        methods. See the --analysis_metric option for details.
  -c MAP_CATEGORIES, --map_categories MAP_CATEGORIES
                        Any mapping categories, such as treatment type, that
                        will be used to group the data in the output iTol
                        table. For example, one category with three types will
                        result in three data columns in the final output. Two
                        categories with three types each will result in six
                        data columns. Default is no categories and all the
                        data will be treated as a single group.
  -a {MRA,NMRA,raw}, --analysis_metric {MRA,NMRA,raw}
                        Specifies which metric is calculated on the abundance
                        data in the OTU table. Available options: MRE - mean
                        relative abundance (Abundance data is normalized by
                        total sample abundance, then averaged across OTU),
                        NMRE - normalized mean relative abundance (MRE
                        normalized by the total MRE across the groups as
                        specified in --map_categories), raw (outputs the
                        actual sequence abundance data for each OTU).