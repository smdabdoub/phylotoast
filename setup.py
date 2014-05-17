from distutils.core import setup

requires = []

setup(
  name = 'qiime-tools',
  packages = ['qiime-tools'],
  version = '1.0',
  description = 'Useful additions to the QIIME project including visualizations and cluster-computing.',
  author = 'Shareef Dabdoub',
  author_email = 'dabdoub.2@osu.edu',
  url = 'https://github.com/smdabdoub/qiime-tools',
  download_url = 'https://github.com/smdabdoub/qiime-tools/tarball/1.0',
  install_requires=requires,
  keywords = ['bioinformatics', 'QIIME', 'microbial ecology', '16S', 'microbiology'],
  license='MIT'
  classifiers=(
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Topic :: Software Development :: Bioinformatics',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Operating System :: UNIX,
        'Operating System :: MacOS X,
        'Operating System :: POSIX :: BSD,
        'Operating System :: POSIX :: Linux',
    ),
)
