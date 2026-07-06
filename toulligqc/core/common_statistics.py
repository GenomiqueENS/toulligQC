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

from math import log

import numpy as np
import pandas as pd


def occupancy_channel(dataframe: pd.DataFrame) -> pd.Series:
    """Compute statistics about the channels of the flowcell.

    Args:
        dataframe: DataFrame containing a ``channel`` column.

    Returns:
        A pd.Series object containing statistics about the channel
        occupancy without the count value.
    """
    total_reads_per_channel = dataframe["channel"].value_counts()
    return pd.DataFrame.describe(total_reads_per_channel)


def compute_LXX(dataframe_dict: dict, x: float) -> int | None:
    """Compute the LXX value of the total sequence length.

    Args:
        dataframe_dict: Dictionary holding an ``all.reads.sequence.length`` entry.
        x: Percentage threshold of the total sequence length.

    Returns:
        The number of reads needed to reach XX% of the total sequence
        length, or None if the threshold is never reached.
    """
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


def compute_NXX(dataframe_dict: dict, x: float) -> int | None:
    """Compute the NXX value of the total sequence length.

    Args:
        dataframe_dict: Dictionary holding an ``all.reads.sequence.length`` entry.
        x: Percentage threshold of the total sequence length.

    Returns:
        The read length at which the cumulative length reaches XX% of the
        total sequence length, or None if the threshold is never reached.
    """
    data = dataframe_dict["all.reads.sequence.length"].dropna().values
    data = np.sort(data)
    half_sum = data.sum() * x / 100
    cum_sum = 0
    for v in data:
        cum_sum += int(v)
        if cum_sum >= half_sum:
            return int(v)


def avg_qual(quals: str) -> float | None:
    """Estimate the mean quality Phred score.

    Args:
        quals: String of per-base quality characters.

    Returns:
        The mean Phred quality score as a float, or None if ``quals`` is empty.
    """
    if quals:
        qscore = -10 * log(
            sum([10 ** ((ord(q) - 33) / -10) for q in quals]) / len(quals), 10
        )
        return round(qscore, 2)
    else:
        return None
