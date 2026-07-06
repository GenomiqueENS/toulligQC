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

"""Plotting methods for :class:`toulligqc.api.TOULLIGQC`.

Every ``plot_*`` method returns a :class:`plotly.graph_objs.Figure`. They are
grouped here as :class:`PlotsMixin`, which :class:`TOULLIGQC` inherits, so the
main module stays focused on the pipeline while the (long, repetitive) list of
plots lives on its own.

The mixin relies on attributes provided by :class:`TOULLIGQC`: ``_ensure_extracted``,
``_config``, ``_result_dict`` and ``_main`` (the input extractor holding the
dataframes).
"""

from toulligqc.graphs import plotly_graph_common as pgc
from toulligqc.graphs import plotly_graph_generator as pgg

# Plots always available, whatever the input type.
STANDARD_PLOTS = (
    "plot_read_count",
    "plot_read_length",
    "plot_yield",
    "plot_read_quality",
    "plot_phred_distribution",
    "plot_channel_performance",
    "plot_length_quality_density",
    "plot_length_over_time",
    "plot_phred_over_time",
    "plot_speed_over_time",
)

# Extra plots only meaningful for a barcoded run.
BARCODE_PLOTS = (
    "plot_barcode_pass_pie",
    "plot_barcode_fail_pie",
    "plot_barcode_length",
)


class PlotsMixin:
    """``plot_*`` methods for :class:`TOULLIGQC` (see module docstring)."""

    # -- helpers ---------------------------------------------------------

    def _figure(self, build):
        """Run ``build`` (a thunk producing exactly one graph) and return its Figure.

        Extraction is triggered *before* ``build`` runs, so the thunk can safely
        reference ``self._main`` and the cached dataframes.
        """
        self._ensure_extracted()
        with pgc.capture_figures() as figures:
            build()
        if not figures:
            raise RuntimeError(
                "The graph did not produce a figure for this run "
                "(the data required by this plot may be missing)."
            )
        return figures[-1]

    def _dataframe_1d(self):
        if not hasattr(self._main, "dataframe_1d"):
            raise AttributeError(
                f"{self._main.get_name()} does not expose a per-read dataframe; "
                "this plot is unavailable for the current input type."
            )
        return self._main.dataframe_1d

    def _barcode_alias(self):
        return self._config.get("barcode_alias", None)

    def _require_barcode(self):
        if str(self._config.get("barcoding", "False")).lower() != "true":
            raise ValueError(
                "Barcode plots require a barcoded run "
                "(pass barcodes=... or samplesheet=...)."
            )

    # -- standard plots --------------------------------------------------

    def plot_read_count(self):
        """Histogram of read counts (all / pass / fail)."""
        return self._figure(lambda: pgg.read_count_histogram(self._result_dict, None))

    def plot_read_length(self):
        """Distribution of read lengths."""
        return self._figure(
            lambda: pgg.read_length_scatterplot(
                self._main.dataframe_dict,
                None,
                getattr(self._main, "read_length_dist_bin_width", None),
            )
        )

    def plot_yield(self):
        """Cumulative yield through time."""
        return self._figure(lambda: pgg.yield_plot(self._dataframe_1d(), None))

    def plot_read_quality(self):
        """Read-quality box plots (all / pass / fail)."""
        return self._figure(
            lambda: pgg.read_quality_multiboxplot(self._main.dataframe_dict, None)
        )

    def plot_phred_distribution(self):
        """PHRED score frequency distribution."""
        return self._figure(
            lambda: pgg.allphred_score_frequency(self._main.dataframe_dict, None)
        )

    def plot_channel_performance(self):
        """Per-channel throughput heatmap of the flow cell."""
        return self._figure(lambda: pgg.plot_performance(self._dataframe_1d(), None))

    def plot_length_quality_density(self):
        """2D density of read length against quality."""
        return self._figure(lambda: pgg.twod_density(self._main.dataframe_dict, None))

    def plot_length_over_time(self):
        """Read length over sequencing time."""
        return self._figure(
            lambda: pgg.sequence_length_over_time(self._main.dataframe_dict, None)
        )

    def plot_phred_over_time(self):
        """PHRED score over sequencing time."""
        return self._figure(
            lambda: pgg.phred_score_over_time(
                self._main.dataframe_dict, self._result_dict, None
            )
        )

    def plot_speed_over_time(self):
        """Sequencing speed (bases/s) over time."""
        return self._figure(
            lambda: pgg.speed_over_time(self._main.dataframe_dict, None)
        )

    # -- barcode-only plots ----------------------------------------------

    def plot_barcode_pass_pie(self):
        """Pie chart of pass-read proportions per barcode."""
        self._require_barcode()
        return self._figure(
            lambda: pgg.barcode_percentage_pie_chart_pass(
                self._main.dataframe_dict,
                self._main.barcode_selection,
                None,
                self._barcode_alias(),
            )
        )

    def plot_barcode_fail_pie(self):
        """Pie chart of fail-read proportions per barcode."""
        self._require_barcode()
        return self._figure(
            lambda: pgg.barcode_percentage_pie_chart_fail(
                self._main.dataframe_dict,
                self._main.barcode_selection,
                None,
                self._barcode_alias(),
            )
        )

    def plot_barcode_length(self):
        """Read-length box plots per barcode."""
        self._require_barcode()
        return self._figure(
            lambda: pgg.barcode_length_boxplot(
                self._main.dataframe_dict,
                None,
                self._barcode_alias(),
            )
        )

    # -- convenience -----------------------------------------------------

    def plots(self) -> dict:
        """Return every applicable plot as a ``{name: Figure}`` dict.

        Barcode plots are only included for barcoded runs, and plots that need
        data unavailable for the current input type are silently skipped.
        """
        self._ensure_extracted()

        names = list(STANDARD_PLOTS)
        if str(self._config.get("barcoding", "False")).lower() == "true":
            names += list(BARCODE_PLOTS)

        figures = {}
        for name in names:
            try:
                figures[name] = getattr(self, name)()
            except (AttributeError, ValueError, RuntimeError, KeyError):
                # Plot not applicable to this run/input type.
                continue
        return figures
