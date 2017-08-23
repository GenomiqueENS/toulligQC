# -*- coding: utf-8 -*-

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

class toullig_conf():
    def __init__(self):
        self._config_dictionary = {}
        self._config_dictionary['dpi'] = '100'
        self._config_dictionary['tmpdir'] = '/tmp'
        self._config_dictionary['barcoding'] = False

    def __getitem__(self, item):
        return self._config_dictionary[item]
        #return self._config_dictionary.__getitem__(item)

    def __setitem__(self, key, value):
        #self._config_dictionary['key'] = value
        self._config_dictionary.__setitem__(key, value)

    def items(self):
        return self._config_dictionary.items()
    def __contains__(self, item):
        ff
    def load(self,conf_path):

        with open(conf_path, 'r') as config_file:
            for line in config_file:
                if line.strip():
                    line = line.replace(" ", "")
                    path_list = line.strip().split('=')
                    self._config_dictionary[path_list[0]] = path_list[1]
