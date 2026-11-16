"""Report generation for ToulligQC.

Once every extractor has run, ``html_report`` renders the accumulated graphs
and statistics into the self-contained HTML report, while
``statistics_generator`` writes the plain-text ``report.data`` dump.
"""

from toulligqc.report.html_report_generator import html_report
from toulligqc.report.report_data_file_generator import (
    add_values_to_unwritten_key,
    statistics_generator,
)

__all__ = [
    "add_values_to_unwritten_key",
    "html_report",
    "statistics_generator",
]
