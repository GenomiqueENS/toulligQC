# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

ToulligQC is a post-sequencing QC tool for Oxford Nanopore sequencers (MinION, etc.). It ingests basecaller output (Guppy/Dorado `sequencing_summary.txt` + `sequencing_telemetry.js`, or FASTQ/uBAM, plus optional FAST5/POD5 for run metadata) and produces a single self-contained HTML report with embedded Plotly graphs plus a plain-text `report.data` statistics file. It supports 1D², barcoded runs, and MinKNOW samplesheets.

## Commands

This project uses `uv` and `pyproject.toml` exclusively (the old `setup.py` was removed).

```bash
uv sync                       # install deps + dev group (ruff, pytest)
uv run toulligqc [options]    # run the CLI
uv run pytest                 # run the test suite
uv run pytest test/test_sequencing_summary_extractor.py::TestName   # single test
uv run ruff check .           # lint
uv run ruff format .          # format (double quotes, 4-space indent)
```

Requires Python >=3.11 (repo pins 3.12 via `.python-version`). Ruff config lives in `pyproject.toml`; `E501` (line length) is intentionally ignored.

The CLI entry point is `toulligqc.toulligqc:main`. Required input is one of: `-a` sequencing summary, `-q` FASTQ, or `-u` uBAM. Use `--force` to overwrite existing output, `--debug`/`--report-only` are hidden dev flags.

## Architecture

The pipeline is an **extractor pipeline driven by a shared mutable config dict and a shared `result_dict`**. Understanding these two dictionaries is the key to the whole codebase.

**1. Config dictionary (`configuration.py::ToulligqcConf`)** — a dict-like object seeded with app defaults, then filled by `_parse_args` from the command line. Every value is coerced to a string. It is passed to every extractor's constructor. Note the quirk: `qscore_threshold()`/`is_default_qscore_threshold()` treat a stored threshold of `"-1"` as "use default 9 / read passes_filtering from file as-is".

**2. Extractor list (`toulligqc.py::_create_extractor_list`)** — based on which inputs were provided, a list of extractor objects is built (telemetry, fast5, pod5, 1dsqr, then exactly one of fastq/bam/sequencing-summary). A `ToulligqcInfoExtractor` is always inserted at index 0. Order matters: earlier extractors' data can be consumed by later ones.

**3. Extractor lifecycle (`toulligqc.py::main`)** — for each extractor, in order, `main` calls the same five methods:
`check_conf()` → `init()` → `extract(result_dict)` → `graph_generation(result_dict)` → `clean(result_dict)`.
Every extractor implements this interface. `extract` writes stats into the shared `result_dict`; `graph_generation` returns a list of `(title, image_path/html, ...)` graph tuples that are accumulated across all extractors; `clean` frees large in-memory dataframes.

**4. result_dict namespacing (`extractor_common.py`)** — extractors never write raw keys. They go through `set_result_value`/`get_result_value`, which prefix keys with the extractor's `get_report_data_file_id()`. `_check_result_key_value` enforces that values are only int/float/str/list/pd.Series/pd.DataFrame. This namespacing is why `report_data_file_generator.py` can flatten everything into `report.data`.

**5. Report generation** — after all extractors run, `html_report_generator.py::html_report` renders the accumulated graphs + `result_dict` into the HTML report (embedding images as data URIs), and `report_data_file_generator.py::statistics_generator` writes the `report.data` text dump.

### Module groups

- **Extractors** (`*_extractor.py`): one per input type — `sequencing_summary_extractor`, `sequencing_summary_onedsquare_extractor` (1D²), `fastq_extractor`, `bam_extractor`, `fast5_extractor`, `pod5_extractor`, `sequencing_telemetry_extractor`, `toulligqc_info_extractor`. Shared helpers in `extractor_common.py`, `fastq_bam_common.py`, `common_statistics.py`.
- **Graph generators** (`plotly_graph_*.py`): `plotly_graph_generator` (standard runs), `plotly_graph_onedsquare_generator` (1D²), with shared plotting helpers in `plotly_graph_common.py` (the largest module). Extractors call these from their `graph_generation`.
- **Report** : `html_report_generator.py`, `report_data_file_generator.py`.
- **`common.py`**: small cross-cutting utilities (e.g. `format_duration`).

### Barcode handling

Barcode parsing lives in `main` (`toulligqc.py`). Barcode names from `--barcodes` (list or `start:end` range) or a samplesheet are normalized to `barcodeNN` via the regex `(BC|RB|NB|BP|BARCODE)(\d{2})`; non-matching names are kept verbatim (custom arrangements). The resolved set goes into `config_dictionary["barcode_selection"]`. Samplesheet aliases populate `barcode_alias`.

## Tests

Tests live in `test/`. `test/config.py` builds `ToulligqcConf` fixtures pointing at `test_data/sequencing_summary/`; `test_sequencing_summary_extractor.py` is the main unit test. `test/toulligqc-compare-reports.py` diffs two generated reports (useful for regression-checking output changes). Sample inputs are under `test_data/` (fast5, fastq, sequencing_summary, samplesheets, integration).

## Release

Versioning is via `release-please` (config in `release-please-config.json`, manifest in `.release-please-manifest.json`); the version single-sources from `toulligqc/version.py`. Use Conventional Commits (`fix:`, `feat:`, `chore:`) — they drive changelog and version bumps.
