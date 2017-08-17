import os
from setuptools import setup, find_packages
from codecs import open
from os import path

version_py = os.path.join(os.path.dirname(__file__), 'poretools', 'version.py')
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

setup(
    name='toulligqc',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='1.0.0',

    description='A data analysis tool of MinION runs',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/GenomicParisCentre/toulligQC',

    # Author details
    author='Lionel Ferrato-Berberian',
    author_email='ferrato@biologie-ens.fr',

    # Choose your license
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],

    # What does your project relate to?
    keywords='MinION analysis report',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=['toulligqc'],

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    #   py_modules=["my_module"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    python_requires='>=3',
    install_requires=['matplotlib', 'seaborn', 'h5py', 'pandas', 'numpy'],


    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    #package_data={'sample': ['package_data.dat']},


    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'toulligqc=toulligqc.toulligqc_main:main',
        ],
    },
)
