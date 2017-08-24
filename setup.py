import os
from setuptools import setup, find_packages
from codecs import open
from os import path

version_py = os.path.join(os.path.dirname(__file__), 'toulligqc', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"', '')\
    .strip()

"""A setuptools based setup module.i

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""


here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
print(version)
setup(
    name='toulligqc',

    version=version,

    description='A data analysis tool of MinION runs',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/GenomicParisCentre/toulligQC',

    # Author details
    author='Lionel Ferrato-Berberian',
    author_email='ferrato@biologie-ens.fr',

    license='GPL V3',
    platforms = 'ALL',
    classifiers=[
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        'License :: OSI Approved :: GPL V3',

        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],

    keywords='MinION analysis report',

    packages=['toulligqc'],
    package_dir={'toulligqc': "toulligqc"},
    package_data={'toulligqc': []},
    zip_safe=False,
    include_package_data=True,

    python_requires='>=3',
    install_requires=['matplotlib>=2.0,<2.1', 'seaborn>=0.7,<0.8', 'h5py>=2.7,<2.8', \
                      'pandas>=0.19,<0.20', 'numpy>=1.12,<1.13'],

    entry_points={
        'console_scripts': [
            'toulligqc=toulligqc.toulligqc:main',
        ],
    },
)
