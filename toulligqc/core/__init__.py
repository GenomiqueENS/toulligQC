"""Core building blocks shared across ToulligQC.

This subpackage holds the configuration object and the cross-cutting helpers
(generic utilities and common statistics) that the extractors, graph
generators and report modules all rely on.
"""

from toulligqc.core.common import (
    find_file_in_directory,
    format_duration,
    is_numpy_1_24,
    set_result_dict_value,
)
from toulligqc.core.common_statistics import (
    avg_qual,
    compute_LXX,
    compute_NXX,
    occupancy_channel,
)
from toulligqc.core.configuration import ToulligqcConf

__all__ = [
    "ToulligqcConf",
    "avg_qual",
    "compute_LXX",
    "compute_NXX",
    "find_file_in_directory",
    "format_duration",
    "is_numpy_1_24",
    "occupancy_channel",
    "set_result_dict_value",
]
