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
    url='https://github.com/GenomiqueENS/toulligQC',

    # Author details
    author='Genomic Paris Centre team',
    author_email='toulligqc@bio.ens.psl.eu',

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

        'Programming Language :: Python :: 3.12'
    ],

    keywords='Nanopore MinION QC report',

    packages=['toulligqc'],
    package_dir={'toulligqc': "toulligqc"},
    package_data={'toulligqc': ['resources/*.css', 'resources/*.js', 'resources/*.png']},
    zip_safe=False,
    include_package_data=True,

    python_requires='>=3.11.0',
    install_requires=['matplotlib>=3.6.3',   'plotly==5.15.0', 'h5py>=3.10.0',
                      'pandas>=2.1.4',       'numpy>=1.26.4',  'scipy>=1.11.4',
                      'scikit-learn>=1.4.1', 'tqdm>=4.66.2',   'pysam>=0.22.0',
                      'pod5>=0.3.10', 'ezcharts==0.7.6'],

    entry_points={
        'console_scripts': [
            'toulligqc=toulligqc.toulligqc:main',
        ],
    },
)
