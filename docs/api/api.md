# Python API reference

ToulligQC ships an in-memory Python API built around the `TOULLIGQC` class. It
runs the same extractor pipeline as the command line, but keeps the results in
memory: statistics are returned as `pandas` DataFrames and every `plot_*` method
returns a Plotly figure. See the [API usage tutorial](../tutorials/api_usage.ipynb)
for a walk-through.

```python
from toulligqc.api import TOULLIGQC

qc = TOULLIGQC(summary="sequencing_summary.txt")
qc.extract()      # run statistics as a DataFrame
qc.plot_yield()   # a Plotly figure
```

## `TOULLIGQC`

::: toulligqc.api.TOULLIGQC
    options:
      members_order: source
      show_source: false

## Plotting methods

::: toulligqc.api.plots.PlotsMixin
    options:
      show_source: false

## Run summary

::: toulligqc.api.summary
    options:
      show_source: false
      members:
        - humanize_metric
        - build_summary
        - SummaryMixin

## Configuration

::: toulligqc.api.config.build_config
    options:
      show_source: false
