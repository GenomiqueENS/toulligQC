import os
from setuptools import setup, find_packages
from codecs import open
from os import path
import pypandoc

version_py = os.path.join(os.path.dirname(__file__), 'toulligqc', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"', '')\
    .strip()

"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""


here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = pypandoc.convert_text(f.read(), 'rst', format='md')

setup(
    name='toulligqc',

    version=version,

    description='A post sequencing QC tool for Oxford Nanopore sequencers',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/GenomicParisCentre/toulligQC',

    # Author details
    author='Genomic Paris Centre team',
    author_email='genomique_bioinfo@biologie.ens.fr',

    license='GPL V3',
    platforms = 'ALL',
    classifiers=[
        'Development Status :: 4 - Beta',

        'Environment :: Console',
        'Operating System :: POSIX',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)',

        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],

    keywords='Nanopore MinION QC report',

    packages=['toulligqc'],
    package_dir={'toulligqc': "toulligqc"},
    package_data={'toulligqc': []},
    zip_safe=False,
    include_package_data=True,

    python_requires='>=3',
    install_requires=['matplotlib>=2.0,<2.1', 'seaborn>=0.7,<0.8', 'h5py>=2.7,<2.8', \
                      'pandas>=0.19,<0.20', 'numpy>=1.12,<1.13', 'pypandoc>=1.3'],

    entry_points={
        'console_scripts': [
            'toulligqc=toulligqc.toulligqc:main',
        ],
    },
)
