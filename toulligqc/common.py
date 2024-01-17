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
#
# First author: Laurent Jourdren
# Maintainer: Laurent Jourdren
# Since version 2.2

import numpy as np
from packaging import version
import glob
import os

def is_numpy_1_24():
    """
    This function checks if Numpy version is later then 1.20
    """
    return version.parse(np.__version__) >= version.parse("1.20")


def format_duration(t):
    """
    Format a time duration
    :param t: time in milliseconds
    :return: the duration in string
    """

    return "{:,d}m{:2.2f}s".format(int(t // 60), t % 60)


def set_result_dict_value(result_dict, key, tracking_id_dict, dict_key):
    """
    Set metadata values from Fast5 or pod5 dict to result_dict
    """
    value = ''
    if dict_key in tracking_id_dict:
        value = tracking_id_dict[dict_key]

    result_dict[key] = value


def find_file_in_directory(source_file, format):
    """
    Looking for a suitable Fast5 or Pod5 file in the source directory.
    :return: The path to the first suitable file in the source directory
    """
    for ext in (format, 'tar.bz2', 'tar.gz'):
        if glob.glob(source_file + '/*.' + ext):
            files_found = os.listdir(source_file)
            if len(files_found) > 0:
                return source_file + files_found[0]

    return None