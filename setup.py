import os
from setuptools import setup, find_packages
from codecs import open

version_py = os.path.join(os.path.dirname(__file__), 'toulligqc', 'version.py')
version = open(version_py).read().strip().split('=')[-1].strip().replace('\'', '').replace('\"', '')

setup(
    name='toulligqc',

    version=version,

    description='A post sequencing QC tool for Oxford Nanopore sequencers',
    long_description='See project website for more information.',

    # The project's main homepage.
    url='https://github.com/GenomicParisCentre/toulligQC',

    # Author details
    author='Genomic Paris Centre team',
    author_email='toulligqc@biologie.ens.fr',

    license='GPL V3',
    platforms='ALL',
    classifiers=[
        'Development Status :: 5 - Production/Stable',

        'Environment :: Console',
        'Operating System :: POSIX',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)',

        'Programming Language :: Python :: 3.11'
    ],

    keywords='Nanopore MinION QC report',

    packages=['toulligqc'],
    package_dir={'toulligqc': "toulligqc"},
    package_data={'toulligqc': ['resources/*.css', 'resources/*.js', 'resources/*.png']},
    zip_safe=False,
    include_package_data=True,

    python_requires='>=3.11.0',
    install_requires=['matplotlib>=3.6.3',   'plotly>=5.15.0', 'h5py>=3.7.0',
                      'pandas>=1.5.3',       'numpy>=1.24.2',  'scipy>=1.10.1',
                      'scikit-learn>=1.2.1', 'tqdm>=4.64.1',   'pysam>=0.21.0'],

    entry_points={
        'console_scripts': [
            'toulligqc=toulligqc.toulligqc:main',
        ],
    },
)
