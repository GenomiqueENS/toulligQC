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


def format_duration(t):
    """
    Format a time duration
    :param t: time in milliseconds
    :return: the duration in string
    """

    return "{:,d}m{:2.2f}s".format(int(t // 60), t % 60)
