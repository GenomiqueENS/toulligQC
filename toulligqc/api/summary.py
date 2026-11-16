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

"""Human-readable run summary for :class:`toulligqc.api.TOULLIGQC`.

The extractors namespace every statistic (e.g.
``basecaller.sequencing.summary.1d.extractor.all.read.length.50%``). This module
turns those keys into plain-English labels ("median read length") and exposes a
curated headline summary, both as a DataFrame and as an inline Plotly table.
"""

import pandas as pd
import plotly.graph_objs as go

# Extractor id prefixes stripped from a metric before humanising its suffix.
_PREFIXES = (
    "basecaller.sequencing.summary.1dsqr.extractor.",
    "basecaller.sequencing.summary.1d.extractor.",
    "sequencing.telemetry.extractor.",
    "toulligqc.info.",
)

# Trailing statistic token -> human label.
_STAT_LABELS = {
    "50%": "median",
    "mean": "mean",
    "std": "std dev",
    "min": "min",
    "max": "max",
    "25%": "Q1",
    "75%": "Q3",
    "count": "count of",
}

# Metric base (prefix + stat token removed) -> human label.
_BASE_LABELS = {
    "all.read.length": "read length",
    "pass.reads.sequence.length": "pass read length",
    "fail.reads.sequence.length": "fail read length",
    "all.read.qscore": "read quality",
    "pass.reads.mean.qscore": "pass read quality",
    "fail.reads.mean.qscore": "fail read quality",
    "channel.occupancy.statistics": "channel occupancy",
    "read.count": "number of reads",
    "read.pass.count": "number of pass reads",
    "read.fail.count": "number of fail reads",
    "read.pass.ratio": "pass ratio",
    "read.fail.ratio": "fail ratio",
    "read.pass.frequency": "pass reads (%)",
    "read.fail.frequency": "fail reads (%)",
    "read.count.frequency": "reads (%)",
    "n50": "N50 read length",
    "l50": "L50",
    "yield": "yield (bases)",
    "run.time": "run duration (s)",
}

# Curated headline metrics: (metric suffix after the extractor prefix, label).
# Order defines the row order of the summary.
SUMMARY_METRICS = (
    ("read.count", "Total reads"),
    ("read.pass.count", "Pass reads"),
    ("read.fail.count", "Fail reads"),
    ("read.pass.frequency", "Pass reads (%)"),
    ("yield", "Total bases (yield)"),
    ("all.read.length.mean", "Mean read length (bp)"),
    ("all.read.length.50%", "Median read length (bp)"),
    ("all.read.length.max", "Max read length (bp)"),
    ("n50", "N50 read length (bp)"),
    ("all.read.qscore.mean", "Mean read quality"),
    ("all.read.qscore.50%", "Median read quality"),
    ("run.time", "Run duration (s)"),
)


def _strip_prefix(metric: str) -> str:
    """Remove the extractor-id prefix from a namespaced metric key."""
    for prefix in _PREFIXES:
        if metric.startswith(prefix):
            return metric[len(prefix) :]
    return metric


def humanize_metric(metric: str) -> str:
    """Turn a namespaced metric key into a plain-English label.

    ``...all.read.length.50%`` -> ``"median read length"``.
    Unknown keys fall back to their dotted suffix with dots replaced by spaces.
    """
    suffix = _strip_prefix(metric)

    # An exact label for the whole suffix wins (e.g. "read.pass.count").
    if suffix in _BASE_LABELS:
        return _BASE_LABELS[suffix]

    # Otherwise split a trailing statistic token: "<base>.<stat>".
    for stat_key, stat_label in _STAT_LABELS.items():
        if suffix.endswith("." + stat_key):
            base = suffix[: -len(stat_key) - 1]
            base_label = _BASE_LABELS.get(base, base.replace(".", " "))
            return f"{stat_label} {base_label}"

    return suffix.replace(".", " ")


def build_summary(stats: pd.DataFrame) -> pd.DataFrame:
    """Build the curated, human-readable summary from an :meth:`extract` DataFrame.

    Args:
        stats: The ``(metric, value)`` DataFrame returned by ``TOULLIGQC.extract``.

    Returns:
        A ``(metric, value)`` DataFrame limited to the headline metrics of
        :data:`SUMMARY_METRICS`, with plain-English labels and numeric values.
        Metrics absent from ``stats`` (e.g. for a different input type) are skipped.
    """
    # Map every metric to its extractor-agnostic suffix for lookup.
    by_suffix = {_strip_prefix(m): v for m, v in zip(stats["metric"], stats["value"])}

    rows = [
        (label, by_suffix[suffix])
        for suffix, label in SUMMARY_METRICS
        if suffix in by_suffix
    ]
    return pd.DataFrame(rows, columns=["metric", "value"])


def _format_value(value) -> str:
    """Format a summary value for display in the table figure."""
    if isinstance(value, (int, float)):
        number = float(value)
        if number.is_integer():
            return f"{int(number):,}"
        return f"{number:,.2f}"
    return str(value)


def summary_table_figure(summary: pd.DataFrame) -> go.Figure:
    """Render a summary DataFrame as an inline Plotly table.

    A table (rather than a bar chart) is used on purpose: summary metrics span
    very different scales (counts, base pairs, quality scores), which a shared
    axis would misrepresent.
    """
    fig = go.Figure(
        data=[
            go.Table(
                columnwidth=[3, 2],
                header=dict(
                    values=["<b>Metric</b>", "<b>Value</b>"],
                    fill_color="#1b3a4b",
                    font=dict(color="white", size=13),
                    align="left",
                ),
                cells=dict(
                    values=[
                        list(summary["metric"]),
                        [_format_value(v) for v in summary["value"]],
                    ],
                    fill_color=[["#f4f7f9", "#e8eef2"] * len(summary)],
                    align="left",
                    height=26,
                ),
            )
        ]
    )
    fig.update_layout(title="Run summary", margin=dict(t=40, b=10, l=10, r=10))
    return fig


class SummaryMixin:
    """``summarise`` / ``plot_summary`` methods for :class:`TOULLIGQC`.

    Relies on :meth:`TOULLIGQC.extract` being available on the same object.
    """

    def summarise(self) -> pd.DataFrame:
        """Return a curated, human-readable summary of the run as a DataFrame.

        Cryptic metric keys such as ``...all.read.length.50%`` become readable
        labels ("Median read length (bp)"). See :func:`build_summary`.
        """
        return build_summary(self.extract())

    def plot_summary(self) -> go.Figure:
        """Return the run summary rendered as an inline Plotly table."""
        return summary_table_figure(self.summarise())
