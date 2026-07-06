"""Plotly graph generators for ToulligQC.

The graph functions are numerous and are meant to be called through their
modules (``plotly_graph_generator`` for standard runs,
``plotly_graph_onedsquare_generator`` for 1D²), with shared plotting helpers in
``plotly_graph_common``. Import the module you need, for example::

    from toulligqc.graphs import plotly_graph_generator as pgg
"""

from toulligqc.graphs import (
    plotly_graph_common,
    plotly_graph_generator,
    plotly_graph_onedsquare_generator,
)

__all__ = [
    "plotly_graph_common",
    "plotly_graph_generator",
    "plotly_graph_onedsquare_generator",
]
