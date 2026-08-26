"""ToulligQC: post-sequencing QC for Oxford Nanopore sequencers.

The package is organised into subpackages by role:

- ``core``: configuration object and cross-cutting helpers.
- ``extractors``: one extractor per input type (telemetry, FAST5, POD5, FASTQ,
  uBAM, sequencing summary, 1D²) plus shared extractor helpers.
- ``graphs``: Plotly graph generators.
- ``report``: HTML report and ``report.data`` file generation.

The CLI entry point is :func:`toulligqc.toulligqc.main`.
"""

from toulligqc.version import __version__

__all__ = ["__version__"]
