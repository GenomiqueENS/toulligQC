# -*- coding: utf-8 -*-
#
#                  ToulligQC development code
#
# This code may be freely distributed and modified under the
# terms of the GNU General Public License version 3 or later
# and CeCILL. This should be distributed with the code. If you
# do not have a copy, see:
#
#      http://www.gnu.org/licenses/gpl-3.0-standalone.html
#      http://www.cecill.info/licences/Licence_CeCILL_V2-en.html
#
# Copyright for this code is held jointly by the Genomic platform
# of the Institut de Biologie de l'École Normale Supérieure and
# the individual authors.
#
# For more information on the ToulligQC project and its aims,
# visit the home page at:
#
#      https://github.com/GenomicParisCentre/toulligQC
#
# For more information on the ToulligQC project and its aims,
# visit the home page at:
#
#      https://github.com/GenomicParisCentre/toulligQC
#
#

import tempfile
from toulligqc import version


class ToulligqcConf:
    """
    Dictionary for the storage of configuration file
    """
    def __init__(self):
        self._config_dictionary = {'app.name': "ToulligQC",
                                   'app.url': "https://github.com/GenomicParisCentre/toulligQC",
                                   'app.version':  version.__version__,
                                   'quiet': 'False',
                                   'dpi': '100',
                                   'tmpdir': tempfile.gettempdir(),
                                   'barcoding': 'False',
                                   'is_quicklaunch': 'False'}

    def __getitem__(self, item):
        return self._config_dictionary[item]

    def __setitem__(self, key, value):
        self._config_dictionary.__setitem__(key, value)

    def items(self):
        return self._config_dictionary.items()

    def __contains__(self, key):
        return key in self.keys()

    def __str__(self):
        return self._config_dictionary.__str__()

    def __len__(self):
        return len(self._config_dictionary)

    def __missing__(self, key):
        raise KeyError(key)

    def __delitem__(self, key):
        del self._config_dictionary[key]

    def __iter__(self):
        for key in self._config_dictionary:
            yield key

    def keys(self):
        return self._config_dictionary.keys()

