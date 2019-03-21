<a href="https://raw.githubusercontent.com/GenomicParisCentre/toulligQC/master/Docs/Toulligqc.png"><img src="https://raw.githubusercontent.com/GenomicParisCentre/toulligQC/master/Docs/Toulligqc.png" align="middle" height="50" width="190" > </a> 
[![PyPI version](https://badge.fury.io/py/toulligqc.svg)](https://badge.fury.io/py/toulligqc) [![Downloads](https://pepy.tech/badge/toulligqc)](https://pepy.tech/project/toulligqc) [![Python 3.6](https://img.shields.io/badge/python-3.5-blue.svg)](https://www.python.org/downloads/release/python-360/)  [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) 


ToulligQC is a program written in Python and developped by the [Genomic facility](https://genomique.biologie.ens.fr/) of the [Institute of Biology of the Ecole Normale Superieure (IBENS)](http://www.ibens.ens.fr/).

This program is dedicated to the QC analyses of Oxford Nanopore runs.
Moreover it is adapted to RNA-Seq along with DNA-Seq and it is compatible with 1Dsquare runs.
It can work with Albacore or Guppy outputs. 
It also needs a single FAST5 file (to catch the flowcell Id and the run date) if a telemetry file is not provided and the Albacore outputted FASTQ file (to compute the sequence statistics).
Flow cells and kits version are retrieved using the telemetry file (Albacore and Guppy) or using the pipeline.log file (Albacore only).
ToulligQC can take barcoding samples into account with a samplesheet.csv describing the barcodes used.

To do so, ToulligQC deals with different file formats: gz, tar.gz, bz2, tar.bz2, FASTQ and FAST5.
This tool will produce a set of graphs, statistic files in txt format and a HTML report.

<a href="https://htmlpreview.github.com/?https://github.com/GenomicParisCentre/toulligQC/blob/master/Docs/report.html" rel="some text">![Report preview](https://raw.githubusercontent.com/GenomicParisCentre/toulligQC/master/Docs/images.png)</a>

Click on the [image](https://htmlpreview.github.com/?https://github.com/GenomicParisCentre/toulligQC/blob/master/Docs/report.html) to see an report example! 

## Authors / Support

Bérengère Laffay, Lionel Ferrato-Berberian, Laurent Jourdren, Sophie Lemoine and Stéphane Le Crom.

Support is availlable on [GitHub issue page](https://github.com/GenomicParisCentre/toulligQC/issues) and at **toulligqc** **at** **biologie.ens.fr**.

## Table of Contents

* 1.[Get ToulligQC](#get-toulligqc)
  * 1.1 [Local installation](#local-installation)
  * 1.2 [PyPi package installation](#pypi-installation) 
  * 1.3 [Docker](#docker)
     *  [Docker image recovery](#docker-image-recovery)    
     *  [Launching Docker image with docker run](#launching-Docker-image-with-docker-run)
     
* 2.[Usage](#usage)
  * 2.1 [Command line](#command-line)
      * [Options](#options)
      * [Examples](#examples)
  * 2.2 [Configuration file](#configuration-file)
  * 2.3 [Sample sheet for barcoded samples](#sample-sheet-for-barcoded-samples)
  
* 3.[Output](#output) 


<a name="get-toulligqc"></a>
## 1. Get ToulligQC 
<a name="local-installation"></a>
### 1.1 Local
This option is also suitable if you are interested in further developments of the package, but requires a little bit more hands-on. Install the dependencies required and clone the repository locally.

```bash
$ git clone https://github.com/GenomicParisCentre/toulligQC.git
# X.X here is the version of ToulligQC to install
$ git checkout vX.X
$ cd toulligqc && python3 setup.py build install
```

* **Requirements**

ToulligQC is written with Python 3.
To run ToulligQC without Docker, you need to install the following Python modules:

* matplotlib
* h5py
* pandas
* seaborn
* numpy
* plotly

<a name="pypi-installation"></a>
### 1.3 Using a PyPi package

ToulligQC can be more easlily installed with a pip package availlable on the PyPi repository. The following command line  will install the latest version of ToulligQC: 
```bash
$ pip3 install toulligqc
```

<a name="docker"></a>
### 1.3 Using Docker
ToulligQC and its dependencies are available through a Docker image. To install docker on your system, go to the Docker website(https://docs.docker.com/engine/installation/). 
Even if Docker can run on Windows or macOS virtual machines, we recommend to run ToulligQC on a Linux host. 
<a name="docker-image-recovery"></a>
* ####  Docker image recovery
An image of ToulligQC is hosted on the Docker hub on the genomicpariscentre repository(genomicpariscentre/toulligqc).

```bash
$ docker pull genomicpariscentre/toulligqc:latest
```


<a name="launching-docker-image-with-a-shell-script"></a>
* ####  Launching Docker image with docker run

```
$ docker run -ti \
             -u $(id -u):$(ig -g) \
             --rm \  
             -v /path/to/fast5/directory:/path/to/telemetry/directory \
             -v /path/to/fastq/directory:/path/to/fastq/directory \
             -v /path/to/basecaller/sequencing/summary/file:/path/to/basecaller/sequencing/summary/file \
             -v /path/to/basecaller/sequencing/telemetry/file:/path/to/basecaller/telemetry/summary/file \
             -v /path/to/samplesheet/file:/path/to/samplesheet/file \
             -v /path/to/configuration/file:/path/to/configuration/file \
             -v /path/to/result/directory:/path/to/result/directory \
             toulligqc:latest 
```
<a name="usage"></a>
## 2. Usage
<a name="command-line"></a>

To run ToulligQC, you need one of your initial Fast5 ONT file (or the telemetry file) and you may also need an Albacore output directory to get the ``` Fastq files ```  and ``` the sequencing_summary.txt ```.
ToulligQC can perform analyses on your data if the directory is organise as following:

```
RUN_ID
├── sequencing_summary.txt
├── pipeline.log
├── configuration.cfg
├── sequencing_telemetry.js
└── workspace
    └── run_id.fastq
```

or 

```
RUN_ID
├── sequencing_summary.txt
├── pipeline.log
├── configuration.cfg
├── sequencing_telemetry.js
└── workspace
    └── pass
        └── run_id.fastq
```
for 1Dsquare analysis

```
RUN_ID
├── sequencing_summary.txt
├── pipeline.log
├── configuration.cfg
├── workspace
│   └── pass
│       └── run_id.fastq
└── 1dsq_analysis
    └── sequencing_1dsq_summary.txt
 ```

### 2.1 Command line

<a name="options"></a>
* #### Options

General Options:
```
toulligqc  [-h] [-c FILE] [-n REPORT_NAME] [-f FAST5_SOURCE]
                [-a SEQUENCING_SUMMARY_SOURCE] [-d SUMMARY_SUMMARY_1DSQR_SOURCE]
                [-p ALBACORE_PIPELINE_SOURCE] [-q FASTQ_SOURCE] [-t TELEMETRY_SOURCE]
                [-o OUTPUT] [-s SAMPLE_SHEET_FILE] [-b] [--quiet] [--version]

optional arguments:
  -h, --help                                                       show this help message and exit
  -c FILE, 
  --conf-file FILE                                                 Specify config file
  -n REPORT_NAME, 
  --report-name REPORT_NAME                                        Report name
  -f FAST5_SOURCE, 
  --fast5-source FAST5_SOURCE                                      Fast5 file source
  -a SEQUENCING_SUMMARY_SOURCE,
  --sequencing-summary-source SEQUENCING_SUMMARY_SOURCE, 
  --albacore-summary-source SEQUENCING_SUMMARY_SOURCE              Basecaller sequencing summary source
  -d SUMMARY_SUMMARY_1DSQR_SOURCE, 
  --sequencing-summary-1dsqr-source SUMMARY_SUMMARY_1DSQR_SOURCE, 
  --albacore-1dsqr-summary-source SUMMARY_SUMMARY_1DSQR_SOURCE     Basecaller 1dsq summary source
  -p ALBACORE_PIPELINE_SOURCE,
  --albacore-pipeline-source ALBACORE_PIPELINE_SOURCE              Albacore pipeline log source
  -q FASTQ_SOURCE,
  --fastq-source FASTQ_SOURCE                                      Fastq file source
  -t TELEMETRY_SOURCE,
   --telemetry-source TELEMETRY_SOURCE                             Telemetry file source
  -o OUTPUT,
   --output OUTPUT                                                 Output directory
  -s SAMPLE_SHEET_FILE,
   --samplesheet-file SAMPLE_SHEET_FILE                            Path to sample sheet file
  -b, --barcoding                                                  Barcode usage
  --quiet                                                          Quiet mode
  --version                                                        show program's version number and exit

```

 <a name="example"></a>
 * #### Examples
 

Example with optional arguments:

```bash
$ toulligqc --report-name FAF0256 \
            --telemetry-source /path/to/basecaller/output/sequencing_telemetry.js \
            --sequencing-summary-source /path/to/basecaller/output/sequencing_summary.txt \
            --sequencing-summary-1dsqr-source /path/to/basecaller/output/sequencing_1dsqr_summary.txt \ (optional)
            --fastq-source /path/to/fastq/source \
            --output /path/to/output/directory \
```


Example with optional arguments to deal with barcoded samples:

```bash
$ toulligqc --report-name FAF0256 \
            --barcoding \
            --telemetry-source /path/to/basecaller/output/sequencing_telemetry.js \
            --sequencing-summary-source /path/to/basecaller/output/sequencing_summary.txt \
            --sequencing-summary-source /path/to/basecaller/output/barcoding_summary_pass.txt \         (optional)
            --sequencing-summary-source /path/to/basecaller/output/barcoding_summary_fail.txt \         (optional)
            --sequencing-summary-1dsqr-source /path/to/basecaller/output/sequencing_1dsqr_summary.txt \ (optional)
            --sequencing-summary-1dsqr-source /path/to/basecaller/output/barcoding_summary_pass.txt \   (optional)
            --sequencing-summary-1dsqr-source /path/to/basecaller/output/barcoding_summary_fail.txt \   (optional)
            --fastq-source /path/to/fastq/source \
            --output /path/to/output/directory \
            --sample-sheet-source /path/to/sample/sheet
``` 
 
Example with optional arguments with a configuration file:

```bash
$ python3 toulligqc.py --run-name FAF0256 \
                       --barcoding \
                       --albacore-pipeline-source /path/to/albacore/pipeline.log \
                       --config-file /path/to/configuration/file/
```


<a name="configuration-file"></a>
### 2.2 Configuration file

A configuration file can be used, the required informations have to be defined as following in the same order :

```ini
[config]
;(path to the directory or sequencing_telemetry.js file)
telemetry_source=/path/to/basecaller/output/directory/or/telemetry/file
;(containing either FAST5, FAST5.tar.gz or FAST5.tar.bz2 files)
fast5_source=/path/to/fast5/directory/or/file
;(Sequencing summary path, if there are multiple sequencing summaries for bacode, a tab character must be used as file separator)
sequencing_summary_source=/path/to/basecaller/sequencing/summary/directory/or/file
;(Sequencing summary path, if there are multiple sequencing summaries for bacode, a tab character must be used as file separator)
sequencing_summary_1dsqr_source=/path/to/basecaller/sequencing/1dsqr/summary/directory/or/file
;(directory where the 1dsqr results are stored)
result_directory =/path/to/result/directory/
fastq_source=/path/to/fastq/directory/or/file 
sample_sheet_file=/path/to/sample/sheet/file
```

The samplesheet directory can be omitted if barcodes were not used in the run.

<a name="sample-sheet-for-barcoded-samples"></a>
### 2.3 Samplesheet
 
The sample sheet file describes the different barcodes and their corresponding samples.
The **Index column is mandatory** and  must contain the Oxford Nanopore Technology **barcode number** like **BC01**.
The other columns are optional but can be useful to define your samples for the following analyses. They can be modified at your convenience.

samplesheet.csv example:

index | Reads | 
------- | ------- 
 BC01 | ```dnacpc14_20170328_FNFAF04250_MN17734_mux_scan_1D_validation_test1_45344_barcode01_template.fastq.bz2``` 

## 3.Output



ToulligQC gives different informations such as:

Found in the HTML report:
- A graph allowing to locate potential flowcell spatial biaises
- A read length histogram adapted to transcripts
- A graph checking that the sequencing was homogeneous during a run
- Graphs representing the phred score frequency
- A set of graphs providing quality, length informations and read counts for each barcode

Found in the report.data log file: 
- The information about ToulligQC execution
- The environment variables
- Full statistics are provided for complementary analyses if needed : The information by modules is retained in a key-value form, the prefix of a key being the report data file id of the module
- The nucleotide rate per read or per read and per barcode if FastQ files have been processed 

Organised in a output directory  like this : 
   
```
RUN_ID
├── report.html
├── statistics
│   └── report.data                                                                                                                                                                                                                                              
└── images
    └── graphes.png
```
