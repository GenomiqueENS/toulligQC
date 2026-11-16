# Getting started

## Installing ToulligQC

`ToulligQC` can be installed from `PyPI` on all OS, for any Python version `>=3.11`.

!!! note "Advice (optional)"
    Install [uv](https://docs.astral.sh/uv/) if you don't have it. 
    uv is a fast Python package and project manager.

Pick the installation method that fits your setup ([uv](#using-uv-recommended) is recommended):

=== "From PyPI"

    ```sh
    pip install toulligqc
    ```

=== "uv (specific version)"

    ``` bash
    git clone https://github.com/GenomiqueENS/toulligQC.git
    cd toulligQC
    git checkout vX.X
    uv sync
    ```

=== "with conda"

    ```bash
    git clone https://github.com/GenomiqueENS/toulligQC.git
    cd toulligQC
    conda env create -f environment.yml
    conda activate toulligqc
    pip install .
    ```

=== "Editable mode"

    ``` bash
    git clone https://github.com/GenomiqueENS/toulligQC.git
    cd toulligQC
    git checkout vX.X
    pip install .
    ```

## Run toulligQC

Run toulligQC with uv:

```bash
uv run toulligqc [options]
```

Or activate the virtual environment first:

```bash
source .venv/bin/activate
toulligqc [options]
```

## Docker



!!! note "Advice (optional)"
    ToulligQC and its dependencies are available as a Docker image hosted on
    [Docker Hub](https://hub.docker.com/r/genomicpariscentre/toulligqc). Even though
    Docker can run on Windows or macOS VMs, running ToulligQC on a Linux host is
    recommended.

Pull the image:

```bash
docker pull genomicpariscentre/toulligqc:latest
```

Run it (mount your input and output paths):

```bash
docker run -ti \
    -u $(id -u):$(id -g) \
    --rm \
    -v /path/to/sequencing_summary:/path/to/sequencing_summary \
    -v /path/to/sequencing_telemetry:/path/to/sequencing_telemetry \
    -v /path/to/result/directory:/path/to/result/directory \
    genomicpariscentre/toulligqc:latest
```

## nf-core module

ToulligQC is available as an [nf-core](https://nf-co.re/docs/usage/introduction)
Nextflow module:

```bash
nf-core modules install toulligqc
```

## Your first report

With a basecaller sequencing summary:

```bash
toulligqc \
    --report-name my_first_run \
    --sequencing-summary-source /path/to/sequencing_summary.txt \
    --output-directory output/
```

ToulligQC writes an `output/my_first_run/` directory containing `report.html`,
`report.data` and an `images/` folder.

Next steps:

- The [command-line guide](tutorials/cli_usage.md) covers every option, input
  layout and the report contents.
- The [Python API tutorial](tutorials/api_usage.ipynb) shows how to run the
  pipeline in memory and get statistics and figures directly.

