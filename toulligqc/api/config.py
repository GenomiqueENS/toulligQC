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

"""Build a :class:`ToulligqcConf` from the API's keyword arguments.

Kept separate from :mod:`toulligqc.api.toulligqc_api` so the (fiddly) mapping of
Python arguments onto the CLI's string-typed configuration lives in one place.
"""

import datetime
import os

from toulligqc import toulligqc as _cli
from toulligqc.core import configuration


def _join_sources(value: str | list[str] | None) -> str | None:
    """Normalise a source argument (``str`` or list of ``str``) to the tab-joined
    string form the extractors expect, or ``None`` when nothing was provided."""
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        return "\t".join(str(v) for v in value)
    return str(value)


def build_config(
    *,
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
) -> configuration.ToulligqcConf:
    """Build a :class:`ToulligqcConf` from keyword arguments.

    Mirrors what ``toulligqc._parse_args`` does for the command line: values are
    coerced to strings, optional flags are only inserted when set, and the barcode
    selection is resolved through the shared
    :func:`toulligqc.toulligqc.resolve_barcode_selection` helper so the API and the
    CLI stay in lock-step.
    """
    config = configuration.ToulligqcConf()

    barcodes = barcodes or ""
    if barcodes or samplesheet:
        barcoding = True

    # Sources (str or list -> tab-joined string) and scalar options.
    values = {
        "sequencing_summary_source": _join_sources(summary),
        "fastq": _join_sources(fastq),
        "bam": _join_sources(bam),
        "sequencing_summary_1dsqr_source": _join_sources(summary_1dsqr),
        "sequencing_telemetry_source": telemetry,
        "fast5_source": fast5,
        "pod5_source": pod5,
        "thread": threads,
        "batch_size": batch_size,
        "threshold": -1 if qscore_threshold is None else qscore_threshold,
        "barcodes": barcodes,
    }
    for key, value in values.items():
        if value is not None:
            config[key] = str(value)

    # Presence-checked options: only insert them when enabled.
    if samplesheet:
        config["samplesheet"] = str(samplesheet)
    if use_alias_for_barcodes:
        config["use_alias_for_barcodes"] = "True"
    if read_length_bin_width is not None:
        config["readlengthdist_binwidth"] = str(read_length_bin_width)

    config["barcoding"] = "True" if barcoding else "False"
    config["quiet"] = "True" if quiet else "False"

    if report_name:
        config["report_name"] = str(report_name)
    else:
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H%M%S")
        config["report_name"] = "Toulligqc-report-" + timestamp

    # The API keeps everything in memory: no image directory and no output files.
    # The graph generators treat ``images_directory=None`` as "return divs, write
    # nothing".
    config["images_directory"] = None
    config["html_report_path"] = None
    config["data_report_path"] = None
    config["report_only"] = "True"

    # Directory sources must end with '/', like the CLI does.
    for key in list(config.keys()):
        value = config.get(key, None)
        if (
            isinstance(value, str)
            and key.endswith("_source")
            and os.path.isdir(value)
            and not value.endswith("/")
        ):
            config[key] = value + "/"

    _cli.resolve_barcode_selection(config)
    return config
