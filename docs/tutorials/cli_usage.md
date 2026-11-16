# Command line interface

ToulligQC supports both RNA-Seq and DNA-Seq data and is compatible with 1D²
runs. It works with **Guppy** and **Dorado** basecalling output. Flow-cell and
kit versions are read from the telemetry file; if no telemetry file is provided,
a single FAST5 file can be used to recover the flow-cell ID and run date.

If a sequencing summary file is not available, ToulligQC also accepts **FASTQ**
or **BAM** files. Inputs may be compressed (`gz`, `tar.gz`, `bz2`, `tar.bz2`).
The tool produces a set of graphs, a plain-text statistics file and an HTML
report.

## Input layouts

A standard run needs the basecaller output `sequencing_summary.txt` (and, ideally,
`sequencing_telemetry.js`):

```
RUN_ID
├── sequencing_summary.txt
└── sequencing_telemetry.js
```

1D² analysis:

```
RUN_ID
├── sequencing_summary.txt
├── sequencing_telemetry.js
└── 1dsq_analysis
    └── sequencing_1dsq_summary.txt
```

For a barcoded run, add the Guppy/Dorado barcoding files
(`barcoding_summary_pass.txt` and `barcoding_summary_fail.txt`), or a single
combined `sequencing_summary_all.txt`. ToulligQC recognises the barcode naming
schemes `BCXX`, `RBXX`, `NBXX` and `barcodeXX` (case-insensitive), where `XX` is
the barcode number.

1D² with barcoding files:

```
RUN_ID
├── sequencing_summary.txt
├── sequencing_telemetry.js
└── 1dsq_analysis
    ├── barcoding_summary_pass.txt
    ├── barcoding_summary_fail.txt
    └── sequencing_1dsq_summary.txt
```

## Options

Required input — provide **one** of `-a`, `-q` or `-u`:

| Option | Description |
|--------|-------------|
| `-a`, `--sequencing-summary-source` | Basecaller sequencing summary source (`.gz`/`.bz2` allowed). |
| `-t`, `--telemetry-source` | Basecaller telemetry file source (`.gz`/`.bz2` allowed). |
| `-f`, `--fast5-source` | FAST5 source (needed if no telemetry file); directory or `tar.gz`/`tar.bz2` allowed. |
| `-p`, `--pod5-source` | POD5 source (needed if no telemetry file); directory or `tar.gz`/`tar.bz2` allowed. |
| `-q`, `--fastq` | FASTQ file (if no sequencing summary); `tar.gz` allowed. |
| `-u`, `--bam` | uBAM file (if no sequencing summary); SAM allowed. |

Optional arguments:

| Option | Description |
|--------|-------------|
| `-s`, `--samplesheet` | Samplesheet (`.csv`) to fill sample names in MinKNOW. |
| `--use-aliases-for-barcodes` | Use the samplesheet `alias` column for barcode names instead of `barcode`. |
| `--thread` | Number of threads. |
| `--batch-size` | Batch size. |
| `--qscore-threshold` | Qscore pass/fail threshold. |
| `-n`, `--report-name` | Report name. |
| `--output-directory` | Output directory. |
| `-o`, `--html-report-path` | Output HTML report path. |
| `--data-report-path` | Output `report.data` path. |
| `--images-directory` | Images directory. |
| `-d`, `--sequencing-summary-1dsqr-source` | Basecaller 1D² summary source. |
| `-b`, `--barcoding` | Enable barcode usage. |
| `-l`, `--barcodes` | Comma-separated barcode list (`BC05,RB09,NB01,barcode10`) or a range (`barcode01:barcode19`). |
| `--quiet` | Quiet mode. |
| `--force` | Overwrite existing output files. |
| `-h`, `--help` | Show help and exit. |
| `--version` | Show the version and exit. |

## Examples

**Sequencing summary alone** — note that the flow-cell ID and run date will be
missing (they come from the telemetry or a FAST5 file):

```bash
toulligqc --report-name summary_only \
          --sequencing-summary-source /path/to/sequencing_summary.txt \
          --html-report-path /path/to/output/report.html
```

**Sequencing summary + telemetry file:**

```bash
toulligqc --report-name summary_plus_telemetry \
          --telemetry-source /path/to/sequencing_telemetry.js \
          --sequencing-summary-source /path/to/sequencing_summary.txt \
          --html-report-path /path/to/output/report.html
```

**Telemetry file + FAST5 files:**

```bash
toulligqc --report-name telemetry_plus_fast5 \
          --telemetry-source /path/to/sequencing_telemetry.js \
          --fast5-source /path/to/fast5_files.fast5.gz \
          --html-report-path /path/to/output/report.html
```

**FASTQ / BAM only:**

```bash
toulligqc --report-name FAF0256 \
          --fastq /path/to/fastq_files.fq.gz \
          --html-report-path /path/to/output/report.html
```

**1D² analysis:**

```bash
toulligqc --report-name FAF0256 \
          --telemetry-source /path/to/sequencing_telemetry.js \
          --sequencing-summary-source /path/to/sequencing_summary.txt \
          --sequencing-summary-1dsqr-source /path/to/sequencing_1dsqr_summary.txt \
          --html-report-path /path/to/output/report.html
```

**Barcoded samples:**

```bash
toulligqc --report-name FAF0256 \
          --barcoding \
          --telemetry-source /path/to/sequencing_telemetry.js \
          --sequencing-summary-source /path/to/sequencing_summary.txt \
          --sequencing-summary-source /path/to/barcoding_summary_pass.txt \
          --sequencing-summary-source /path/to/barcoding_summary_fail.txt \
          --html-report-path /path/to/output/report.html \
          --barcodes BC01,BC02,BC03
```

## Sample data

[Sample raw data](http://outils.genomique.biologie.ens.fr/leburon/downloads/toulligqc-example/toulligqc_demo_data.tar.bz2)
is provided to try the software. It was generated on a MinION MkIb with an
R9.4.1 flow cell (FLO-MIN106) in 1D (SQK-LSK108) mode with barcoded samples
(BC01–BC05 and BC07); acquisition with MinKNOW 1.11.5, basecalling with Guppy
3.2.4.

Download the demo scripts:

```bash
wget https://raw.githubusercontent.com/GenomiqueENS/toulligQC/refs/heads/master/help/demo/run-toulligqc-demo-with-docker.sh
wget https://raw.githubusercontent.com/GenomiqueENS/toulligQC/refs/heads/master/help/demo/run-toulligqc-demo.sh
chmod +x run-toulligqc-demo*.sh
```

Run the demo in Docker:

```bash
./run-toulligqc-demo-with-docker.sh
```

…or with a local install:

```bash
./run-toulligqc-demo.sh
```

Or launch ToulligQC manually on the sample data:

```bash
toulligqc \
    --report-name               'ToulligQC Demo Data' \
    --barcoding \
    --telemetry-source          sequencing_telemetry.js \
    --sequencing-summary-source sequencing_summary.txt \
    --sequencing-summary-source barcoding_summary_pass.txt \
    --sequencing-summary-source barcoding_summary_fail.txt \
    --barcodes                  BC01:BC07 \
    --output-directory          output
```

## Output

If `--output-directory` / `--html-report-path` are not given, ToulligQC writes
its output to the current directory; without `--report-name` a default name is
used.

**HTML report** (path set with `--html-report-path`), containing:

- useful information about the sequencing run;
- read-count and read-length histograms per read type;
- a graph checking that sequencing was homogeneous over the run;
- a graph to locate potential flow-cell spatial biases;
- PHRED-score and density distributions across read types;
- length / speed / quality / read count over sequencing time;
- per-barcode quality, length and read counts.

**`report.data`** log file (path set with `--data-report-path`), containing:

- information about the ToulligQC execution and environment;
- full per-module statistics in key-value form (the key prefix is the module's
  report-data-file id), for complementary analyses;
- the nucleotide rate per read.

With a directory output (the default), the layout is:

```
RUN_ID
├── report.html
├── report.data
└── images
    ├── plots.html
    └── plot.png
```
