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

"""Object-oriented, in-memory front end to the ToulligQC extractor pipeline.

The :class:`TOULLIGQC` class wraps the exact same extractor/graph machinery used
by the command line, but instead of writing an HTML report to disk it keeps the
computed statistics and the underlying dataframes in memory and hands them back
to the caller: :meth:`TOULLIGQC.extract` returns a :class:`pandas.DataFrame` of
run statistics, and every ``plot_*`` method (defined in :mod:`toulligqc.api.plots`)
returns a :class:`plotly.graph_objs.Figure` that renders inline in a notebook.

The heavy work (parsing the inputs) runs lazily and exactly once, the first time
:meth:`extract` or any ``plot_*`` method is called; every later call reuses the
cached dataframes.
"""

import os

import pandas as pd

from toulligqc import toulligqc as _cli
from toulligqc.api.config import build_config
from toulligqc.api.plots import PlotsMixin
from toulligqc.api.summary import SummaryMixin
from toulligqc.core import configuration


class TOULLIGQC(PlotsMixin, SummaryMixin):
    """In-memory ToulligQC run.

    Args:
        summary: ``sequencing_summary.txt`` path(s) (str or list). Mutually
            exclusive with ``fastq``/``bam``.
        fastq: FASTQ path(s), used when no sequencing summary is available.
        bam: uBAM path(s), used when no sequencing summary is available.
        telemetry: ``sequencing_telemetry.js`` path (optional, adds run metadata).
        fast5: FAST5 source (file, directory or archive) for run metadata.
        pod5: POD5 source (file, directory or archive) for run metadata.
        summary_1dsqr: 1D² sequencing summary path(s).
        barcoding: Force barcode mode. Automatically enabled when ``barcodes`` or
            ``samplesheet`` is given.
        barcodes: Comma-separated barcode list or ``start:end`` range
            (e.g. ``"barcode01:barcode12"``).
        samplesheet: MinKNOW samplesheet (.csv) used to fill sample names.
        use_alias_for_barcodes: Use the samplesheet ``alias`` column for names.
        qscore_threshold: Qscore pass/fail threshold (default 9).
        read_length_bin_width: Bin width (bp) for read-length distribution plots.
        threads: Number of worker threads for the extractors.
        batch_size: Batch size for the extractors.
        report_name: Name used for the run (defaults to a timestamp).
        quiet: Suppress the per-extractor progress messages (default ``True``).
    """

    def __init__(
        self,
        summary: str | list[str] | None = None,
        fastq: str | list[str] | None = None,
        bam: str | list[str] | None = None,
        telemetry: str | None = None,
        fast5: str | None = None,
        pod5: str | None = None,
        summary_1dsqr: str | list[str] | None = None,
        barcoding: bool = False,
        barcodes: str | None = None,
        samplesheet: str | None = None,
        use_alias_for_barcodes: bool = False,
        qscore_threshold: int | None = None,
        read_length_bin_width: float | None = None,
        threads: int = 2,
        batch_size: int = 500,
        report_name: str | None = None,
        quiet: bool = True,
    ) -> None:
        if summary is None and fastq is None and bam is None:
            raise ValueError(
                "TOULLIGQC requires one input source: 'summary', 'fastq' or 'bam'."
            )

        self._config = build_config(
            summary=summary,
            fastq=fastq,
            bam=bam,
            telemetry=telemetry,
            fast5=fast5,
            pod5=pod5,
            summary_1dsqr=summary_1dsqr,
            barcoding=barcoding,
            barcodes=barcodes,
            samplesheet=samplesheet,
            use_alias_for_barcodes=use_alias_for_barcodes,
            qscore_threshold=qscore_threshold,
            read_length_bin_width=read_length_bin_width,
            threads=threads,
            batch_size=batch_size,
            report_name=report_name,
            quiet=quiet,
        )

        # Lazily populated by _ensure_extracted().
        self._extractors: list | None = None
        self._main = None
        self._result_dict: dict = {}

    # -- pipeline ---------------------------------------------------------

    def _ensure_extracted(self) -> "TOULLIGQC":
        """Run the extractor pipeline once and cache the results.

        Runs ``check_conf`` / ``init`` / ``extract`` for every extractor (the same
        order as the CLI) but deliberately skips ``graph_generation`` and
        ``clean`` so the in-memory dataframes stay available for the ``plot_*``
        methods.
        """
        if self._extractors is not None:
            return self

        extractors = _cli._create_extractor_list(self._config)

        for extractor in extractors:
            ok, message = extractor.check_conf()
            if not ok:
                raise ValueError(
                    f"Invalid configuration for {extractor.get_name()}: {message}"
                )

        result_dict: dict = {}
        main = None
        for extractor in extractors:
            extractor.init()
            extractor.extract(result_dict)
            # The extractor exposing per-read dataframes is the input extractor
            # (sequencing summary / FASTQ / uBAM); keep the last one seen.
            if hasattr(extractor, "dataframe_dict"):
                main = extractor

        if main is None:
            raise RuntimeError("No data extractor produced a dataframe to plot.")

        self._extractors = extractors
        self._result_dict = result_dict
        self._main = main
        return self

    @property
    def config(self) -> configuration.ToulligqcConf:
        """The resolved :class:`ToulligqcConf` backing this run."""
        return self._config

    # -- statistics -------------------------------------------------------

    def extract(self) -> pd.DataFrame:
        """Run the pipeline (if needed) and return the run statistics.

        Returns:
            A tidy two-column DataFrame (``metric``, ``value``) of every scalar
            statistic collected by the extractors. Non-scalar entries (series and
            dataframes) are exposed through :attr:`stats` and :attr:`reads`.
        """
        self._ensure_extracted()
        rows = [
            (key, value)
            for key, value in sorted(self._result_dict.items())
            if isinstance(value, (int, float, str))
        ]
        return pd.DataFrame(rows, columns=["metric", "value"])

    @property
    def stats(self) -> dict:
        """The full, namespaced ``result_dict`` produced by the extractors."""
        self._ensure_extracted()
        return self._result_dict

    @property
    def reads(self) -> pd.DataFrame:
        """The per-read dataframe (``dataframe_1d``) of the input extractor.

        Raises:
            AttributeError: If the active input type does not expose per-read data.
        """
        self._ensure_extracted()
        if not hasattr(self._main, "dataframe_1d"):
            raise AttributeError(
                f"{self._main.get_name()} does not expose a per-read dataframe."
            )
        return self._main.dataframe_1d

    # -- report -----------------------------------------------------------

    def report(self, path: str, force: bool = False) -> str:
        """Write the full standalone HTML report to ``path``.

        This reuses the exact HTML generator used by the CLI.

        Args:
            path: Destination ``.html`` file.
            force: Overwrite ``path`` if it already exists.

        Returns:
            The ``path`` that was written.
        """
        from toulligqc.report import html_report_generator

        self._ensure_extracted()

        if os.path.exists(path) and not force:
            raise FileExistsError(f"{path} already exists (pass force=True).")

        graphs: list = []
        for extractor in self._extractors:
            graphs.extend(extractor.graph_generation(self._result_dict))

        self._config["html_report_path"] = path
        html_report_generator.html_report(self._config, self._result_dict, graphs)
        return path

    def __repr__(self) -> str:
        state = "extracted" if self._extractors is not None else "pending"
        return f"<TOULLIGQC report_name={self._config['report_name']!r} ({state})>"
