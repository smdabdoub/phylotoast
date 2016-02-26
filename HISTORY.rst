.. :changelog:

Release History
---------------

1.3.0 (2016-02-26)
++++++++++++++++++

**Addition and Improvements**

- New analysis/visualization programs:

  - LDA.py (Uses Linear Discriminant Analysis to visualize group-wise compositional differences)
  - diversity.py (Visualize alpha diversity as KDE-smoothed histograms)
- New API module

  - graph_util (Adds convenience methods for data visualization/graphing)

- biom_calc module:

  - Switched to using the official biom-format python library for parsing and handline BIOM files
  - Added arcsine square root transformation to the API and added as an option to several scripts:

    - iTol.py
    - LDA.py


- PCoA.py

  - Support for the new UniFrac format in QIIME 1.9
  - Adds the following command-line options:

    - Set font size
    - Set distance between axes and their labels
    - Set the azimuth and elevation for 3D plots
    - Set the graph styling to look like ggplot2 defaults
    - Add Sample ID annotations to graphed points
    - Replaces the --dpi option with --fig_size to set the figure size directly

- LDA.py: Added the option to save the LDA-transformed data out as a csv file.
- diversity.py

  - Added group coloring by specify a column in the mapping file
  - Support for multiple diversity metrics
  - Significance testing is now performed to compare alpha diversity between groups
    - Wilcoxon Signed Rank Test for two-group comparison
    - Kruskal-Wallis H-test for comparing three or more groups

- biom_relative_abundance.py: Support for BIOM-format input.
- biom_calc.py: Updated raw_abundance() function to handle BIOM format version 2.x input files.
- iTol.py: Updated BIOM file handling for iTol script (v2.x files).

- PhyloToAST documentation moved from separate repository into the 'docs' folder in the main code repository
- Overhauled the main page for the project documentation
- Added step-by-step instructions to the documentation for many of the executable scripts


**Bug Fixes**

- Prevents duplicate entries from being passed on into the output from restrict_repset.py.
- Fixes the --save_calculations option in diversity.py. Now prints out all metrics in a single column with the second column being group membership.
- Adds pad_taxonomy method (otu_condense.py) to ensure all taxonomy strings have all levels down to species.
- Cleans up error message for line parsing in prune_taxonomy. Makes file error message more readable.

1.2.0 (2015-01-29)
++++++++++++++++++

**Improvements**

- Renamed "QIIME-tools" to "PhyloToAST".
- Added travis.yml file to implement continuous integration using Travis IC.
- Modified the README file to conform with the name change.
- Added PCOA.py script allowing users to create 2D or 3D plots.
- Added restrict_repset.py script, which allows users to create a new OTU representative set file from a filtered BIOM table.
- Updated imports in all the scripts to use "phylotoast" instead of "qiime_tools".
- Improved error handling mechanism in all binary scripts.

  - Added try...catch blocks for ImportError and IOError for all third party imports and user required input file handles.
  - Added user friendly error message to enable quick debugging.
  
- Improved and updated all API scripts with new functions.

  - Binary scripts now import API scripts and use the functions from API, thus, reducing binary script complexities.
  
- Modified all scripts to conform to PEP8 format using flake8 module.

1.1.2 (2014-09-19)
++++++++++++++++++

**Improvements**

- Cleans up the abundance calculation code and moves the output_relative_abundance() method to the biom_relative_abundance.py script.
- Modifies name extraction from each OTU level so that underscores in the name (after the initial double underscore) are allowed.
- Adds more detail and project badges (downloads,version,format) to the README.

**Bug Fixes**

- Restores accidental reversion of commit 9413c64 introduced by commit 91d6a39.

1.1.1 (2014-09-16)
++++++++++++++++++

The most notable change is the project is now structured to be installable via setuptools and compatible with (and available on) PyPI:

$ pip install qiime-setuptools

Or if you prefer to build from source:

$ python setup.py install


**Additions and Improvements**

- SLURM template added for use with the multi_parallel_pick_otus script.
- New otu_calc module for working with abundance tables (BIOM) and providing calculations such as Shannon Diversity and relative abundance.
- otu_tax_name.py now handles classic (tsv) BIOM tables in addition to simple lists of OTU IDs.
- New script filter_keep_otus_by_sample.py as the opposite of the QIIME script filter_otus_by_sample.py

**Bug Fixes**

- In iTol.py: fixes an issue where some or all of the OTU IDs in a newick-format tree have single-quotes surrounding them ('OTU_ID'), causing the script to crash.
- otu_name_biom() in the otu_calc method now returns a value as it was meant to


1.0.0 (2014-05-17)
++++++++++++++++++

First release!
