from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

requires = ['numpy', 'scipy', 'matplotlib >= 1.5.0', 'biopython >= 1.60',
            'scikit-bio', 'scikit-learn', 'pandas', 'statsmodels',
            'palettable', 'biom-format >= 2.1.5', 'h5py']

scripts = ['bin/LDA.py',
           'bin/LDA_bubble.py',
           'bin/PCoA.py',
           'bin/PCoA_bubble.py'
           'bin/assign_taxonomy_by_blast_result.py',
           'bin/barcode_filter.py',
           'bin/biom_relative_abundance.py',
           'bin/condense_workflow.py',
           'bin/core_overlap_plot.py',
           'bin/diversity.py',
           'bin/extract_shared_or_unique_otuids.py',
           'bin/filter_biom.py',
           'bin/filter_keep_otus_by_sample.py',
           'bin/filter_rep_set.py',
           'bin/iTol.py',
           'bin/merge_otu_results.py',
           'bin/multi_parallel_pick_otus.py',
           'bin/multi_qsub.py',
           'bin/network_plots_gephi.py',
           'bin/otu_condense.py',
           'bin/otu_to_tax_name.py',
           'bin/pick_otus_condense.py',
           'bin/primer_average.py',
           'bin/prune_otus.py',
           'bin/restrict_repset.py',
           'bin/split_sequence_data.py',
           'bin/transpose_biom.py']

setup(
  name='phylotoast',
  version='1.3.0',
  description='Tools for phylogenetic data analysis including visualization and cluster-computing support. Integrates with the QIIME pipeline. http://phylotoast.org',
  author='Shareef M. Dabdoub',
  author_email='dabdoub.2@osu.edu',
  packages=['phylotoast', 'phylotoast.test'],
  include_package_data=True,
  url='https://github.com/smdabdoub/phylotoast',
  install_requires=requires,
  scripts=scripts,
  data_files=[('data', ['data/pbs_job_template.pbs', 'data/slurm_job_template.sbatch'])],
  keywords=['bioinformatics', 'QIIME', 'phylotoast', 'microbial ecology', '16S', 'microbiology'],
  license='MIT',
  classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Operating System :: Unix',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: BSD',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Microsoft :: Windows',
    ],
)
