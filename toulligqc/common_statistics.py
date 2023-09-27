import pandas as pd
from math import log

def occupancy_channel(dataframe):
    """
    Statistics about the channels of the flowcell
    :return: pd.Series object containing statistics about the channel occupancy without count value
    """
    total_reads_per_channel = pd.value_counts(dataframe["channel"])
    return pd.DataFrame.describe(total_reads_per_channel)


def compute_LXX(dataframe_dict, x):
    """Compute LXX value of total sequence length"""
    data = dataframe_dict["all.reads.sequence.length"].dropna().values
    data.sort()
    half_sum = data.sum() * x / 100
    cum_sum = 0
    count = 0
    for v in data:
        cum_sum += v
        count += 1
        if cum_sum >= half_sum:
            return count


def compute_NXX(dataframe_dict, x):
    """Compute NXX value of total sequence length"""
    data = dataframe_dict["all.reads.sequence.length"].dropna().values
    data.sort()
    half_sum = data.sum() * x / 100
    cum_sum = 0
    for v in data:
        cum_sum += v
        if cum_sum >= half_sum:
            return int(v)


def avg_qual(quals):
        """
        Estimates mean quality Phred score
        return: float
        """
        if quals:
            qscore =  -10 * log(sum([10**((ord(q)-33) / -10) for q in quals]) / len(quals), 10)
            return round(qscore, 2)
        else:
            return None

