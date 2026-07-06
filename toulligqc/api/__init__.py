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

"""Programmatic (notebook-friendly) API for ToulligQC.

Example
-------
>>> from toulligqc.api import TOULLIGQC
>>> qc = TOULLIGQC(summary="path/to/sequencing_summary.txt")
>>> df = qc.extract()          # tidy DataFrame of run statistics
>>> qc.plot_yield()            # a plotly Figure, renders inline in Jupyter
"""

from toulligqc.api.toulligqc_api import TOULLIGQC

__all__ = ["TOULLIGQC"]
