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
#      https://github.com/GenomiqueENS/toulligQC
#

import numpy as np
import pandas as pd
from math import log


def occupancy_channel(dataframe):
    """
    Statistics about the channels of the flowcell
    :return: pd.Series object containing statistics about the channel occupancy without count value
    """
    total_reads_per_channel = dataframe["channel"].value_counts()
    return pd.DataFrame.describe(total_reads_per_channel)


def compute_LXX(dataframe_dict, x):
    """Compute LXX value of total sequence length"""
    data = dataframe_dict["all.reads.sequence.length"].dropna().values
    data = np.sort(data)
    half_sum = data.sum() * x / 100
    cum_sum = 0
    count = 0
    for v in data:
        cum_sum += int(v)
        count += 1
        if cum_sum >= half_sum:
            return count


def compute_NXX(dataframe_dict, x):
    """Compute NXX value of total sequence length"""
    data = dataframe_dict["all.reads.sequence.length"].dropna().values
    data = np.sort(data)
    half_sum = data.sum() * x / 100
    cum_sum = 0
    for v in data:
        cum_sum += int(v)
        if cum_sum >= half_sum:
            return int(v)


def avg_qual(quals):
    """
    Estimates mean quality Phred score
    return: float
    """
    if quals:
        qscore = -10 * log(
            sum([10 ** ((ord(q) - 33) / -10) for q in quals]) / len(quals), 10
        )
        return round(qscore, 2)
    else:
        return None
