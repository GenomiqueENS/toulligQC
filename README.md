# ToulligQC

ToulligQC is a program written in Python and developped by the [Genomic facility](https://genomique.biologie.ens.fr/) of the [Institute of Biologie of the Ecole Normale Superieure (IBENS)](http://www.ibens.ens.fr/).

This program is dedicated to the QC analyses of Oxford Nanopore runs.
Moreover it is adapted to RNA-Seq along with DNA-Seq.
It partly relies on the summary.txt file produced during the basecalling process by the Oxford Nanopore basecaller, Albacore.
It also needs a single FAST5 file (to catch the flowcell Id and the run date) and the Albacore outputted FASTQ file (to compute the sequence statistics).
ToulligQC can take barcoding samples into account with a samplesheet.csv describing the barcodes used.

To do so, ToulligQC deals with different file formats: gz, tar.gz, bz2, tar.bz2, FASTQ and FAST5.
This tool will produce a set of graphs, statistic files in txt format and a HTML report.

<a href="https://htmlpreview.github.com/?https://github.com/GenomicParisCentre/toulligQC/blob/master/report.html" rel="some text">![Report preview](https://raw.githubusercontent.com/GenomicParisCentre/toulligQC/master/report.png)</a>

Click on the [image](https://htmlpreview.github.com/?https://github.com/GenomicParisCentre/toulligQC/blob/master/report.html) to see an report example! 

## Table of Contents

* 1.[Get ToulligQC](#get-toulligqc)
  * 1.1 [Local installation](#local-installation)
  * 1.2 [PyPi installation](#pypi-installation) 
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

<a name="pypi-installation"></a>
### 1.3 Using PyPi

ToulligQC can be more easlilly installed with a pip package availlable on the PyPi repository:
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
             -v /path/to/fast5/directory:/path/to/fast5/file \
             -v /path/to/fastq/directory:/path/to/fastq/file \
             -v /path/to/albacore/sequencing/summary/file:/path/to/albacore/sequencing/summary/file \ 
             -v /path/to/samplesheet/file/:/path/to/samplesheet/file/ \
             -v /path/to/configuration/file:/path/to/configuration/file \
             -v /path/to/result/directory:/path/to/result/directory \
             toulligqc:latest 
```
<a name="usage"></a>
## 2. Usage
<a name="command-line"></a>

To run ToulligQC, you need one of your initial Fast5 ONT file and you may also need an Albacore output directory to get the ``` Fastq files ```  and ``` the sequencing_summary.txt ```.
ToulligQC can perform analyses on your data if the directory is organise as following :

```
RUN_ID
├── sequencing_summary.txt
├── pipeline.log
├── configuration.cfg
└── workspace
    └── run_id.fastq
```

or 

```
RUN_ID
├── sequencing_summary.txt
├── pipeline.log
├── configuration.cfg
└── workspace
    └── pass
        └── run_id.fastq
```
 
### 2.1 Command line

<a name="options"></a>
* #### Options

General Options:
```
toulligqc [-h] [-n RUN_NAME] [-b] [-c CONFIG_FILE] [-f FAST5_SOURCE]
               [-a ALBACORE_SUMMARY_SOURCE] [-q FASTQ_SOURCE]
               [-o OUTPUT_DIRECTORY] [-s SAMPLE_SHEET_SOURCE]

               
optional arguments:

  -h, --help                       Show this help message and exit
  -n RUN_NAME, --run-name          Run name                   
  -b, --barcoding                  Search for barcodes to demultiplex sequencing data
  -c CONFIG_FILE, --config-file    Optional configuration file to use
  -f FAST5_SOURCE, --fast5-source  Fast5 file source (.fast5 or fastq.tar.bz2 format) 
                       
  -a ALBACORE_SUMMARY_SOURCE, --albacore-summary-source      Albacore summary source (.txt format)
  -q FASTQ_SOURCE, --fastq-source                            Fastq file source (.fastq or fastq.tar.gz format)                  
  -o OUTPUT_DIRECTORY, --output                              Path to save the output directory
  -s SAMPLE_SHEET_SOURCE, --sample-sheet-source              Path to samplesheet to take barcodes into account 
```

 <a name="example"></a>
 * #### Examples
 

Example with optional arguments:

```bash
$ python3 toulligqc.py --run-name FAF0256 \
                       --fast5-source /path/to/fast5/source \
                       --albacore-summary-source /path/to/albacore/sequencing_summary.txt \
                       --fastq-source /path/to/fastq/source \
                       --output /path/to/output/directory \
```


Example with optional arguments to deal with barcoded samples:

```bash
$ python3 toulligqc.py --run-name FAF0256 \
                       --barcoding \
                       --fast5-source /path/to/fast5/source \
                       --albacore-summary-source /path/to/albacore/sequencing_summary.txt \
                       --fastq-source /path/to/fastq/source \
                       --output /path/to/output/directory \
                       --sample-sheet-source /path/to/sample/sheet
``` 
 
Example with optional arguments with a configuration file:

```bash
$ python3 toulligqc.py --run-name FAF0256 \
                       --barcoding \
                       --config-file /path/to/configuration/file/
```


<a name="configuration-file"></a>
### 2.2 Configuration file

A configuration file can be used, the required informations have to be defined as following in the same order :

```ini
[config]
;(containing either FAST5, FAST5.tar.gz or FAST5.tar.bz2 files)
fast5_source=/path/to/fast5/directory/or/file
albacore_summary_source=/path/to/albacore/sequencing/summary/directory/or/file
;(directory where the results are stored)
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

index | ``` Reads ``` | 
------- | ------- 
 BC01| ```dnacpc14_20170328_FNFAF04250_MN17734_mux_scan_1D_validation_test1_45344_barcode01_template.fastq.bz2 ``` 

## 3.Output



ToulligQC gives different informations such as:

- A graph allowing to locate potential flowcell spatial biaises
- A read length histogram adapted to transcripts
- A graph checking that the sequencing was homogeneous during a run
- Graphs representing the phred score frequency
- A set of graphs providing quality, length informations and read counts for each barcode
- Full statistics are provided in a text file for complementary analyses if needed

in a output directory organised like this : 
   
```
RUN_ID
├── report.html
├── statistics
│   ├── run_statistics_file.txt                                                                                                                                                                                                                                                               
│   └── save_result_statistics.txt 
│   └── run_id.fastq
└── images
    └── graphes.png
```


