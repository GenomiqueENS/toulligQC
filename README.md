# ToulligQC
This program is dedicated to the QC analyses of Oxford Nanopore runs, barcoded or not. It requires a design file describing the barcodes used if the run was barcoded. It partly relies on the log file produced during the basecalling process by the Oxford Nanopore basecaller, Albacore. It also needs a single FAST5 file (to catch the flowcell Id and the run date) and the the Albacore outputted FASTQ file (to compute the sequence statistics). To do so, ToulligQC deals with different file formats: gz, tar.gz, bz2, tar.bz2, fastq and fast5.
This program will produce a set of graphs and statistic files and a report in pdf, html and docx formats.

## Table of Contents

* 1.[Get ToulligQC](#get-toulligqc)
  * 1.1 [Docker](#docker)
     *  [Docker image recovery](#docker-image-recovery)
  
     *  [Launching docker image with a shell script](#launching-docker-image-with-a-shell-script)
     
     *  [Launching Docker image with docker run](#launching-Docker-image-with docker-run)
     
  * 1.2 [Local installation](#local-installation)
* 2.[Usage](#usage)
    * 2.1 [Command line](#command-line)

      * [Options](#options)
  
      * [Example](#example)

     *  2.2 [Configuration file](#configuration-file)
  
     * 2.3 [Sample sheet for barcoded samples](#sample-sheet-for-barcoded-samples)
* 3.[Output](#output) 

<a name="get-toulligqc"></a>
## 1. Get ToulligQC 
<a name="docker"></a>
### 1.1 Docker
ToulligQC and its dependencies are available through a Docker image. To install docker on your system, go to the Docker website. Even if Docker can run on Windows or macOS virtual machines, we recommend to run ToulligQC on a Linux host. 
<a name="docker-image-recovery"></a>
* ####  Docker image recovery
An image of ToulligQC is hosted on the Docker hub on the genomicpariscentre repository(genomicpariscentre/toulligqc).

```$ docker pull genomicpariscentre/toulligqc:latest ```

<a name="launching-docker-image-with-a-shell-script"></a>
   * #### Launching docker image with a shell script
A shell script called read_file.sh is provided to launch the image and mount  automatically the directories contained in the configuration file. 

Example:
```$ ./read_file.sh /path/to/configuration/file ```

<a name="launching-docker-image-with-a-shell-script"></a>
* ####  Launching Docker image with docker run
```$ docker run -ti --rm  -v /path/to/result/directory/

-v /path/to/fast5/directory\

-v /path/to/fastq/directory\

-v /path/to/design/file/\

-v /path/to/configuration/file/
 (-v /path/to/sequencing/summary/file if not include in fastq file directory)\

toulligqc:latest ```
 
 <a name="local-installation"></a>
#### 1.2 Local
This option is also suitable if you are interested in further developing the package, but requires a little bit more hands-on. Install the dependencies required and clone the repository locally.

```$ git clone https://github.com/GenomicParisCentre/toulligQC.git```

* **Requirements**

To run ToulligQC without Docker, you need to install the following softwares:
* matplotlib
* h5py
* pandas
* seaborn
* numpy

On Debian/Ubuntu, you can install requirements using the 'apt-get' command, here is an example: 

```$ sudo apt-get install matplotlib```

<a name="usage"></a>
## 2. Usage
<a name="command-line"></a>
### 2.1 Command line

<a name="options"></a>
* #### Options

The run name is indicated before the file extension in the FAST5 and FASTQ files.
Example:

usage: ```main.py [-h] [-n RUN_NAME] [-b] [-c CONFIG_FILE] [-f FAST5_SOURCE]
               [-a ALBACORE_SUMMARY_SOURCE] [-q FASTQ_SOURCE]
               [-o OUTPUT_DIRECTORY] [-s SAMPLE_SHEET_SOURCE]`

               
optional arguments:


  -h, --help            show this help message and exit
 
  -n RUN_NAME, --run_name RUN_NAME
                        Run name
                        
  -b, --barcoding         Barcode usage
  
  -c CONFIG_FILE, --config_file CONFIG_FILE
                        Path to the configuration file
                        
  -f FAST5_SOURCE, --fast5-source FAST5_SOURCE
                        Fast5 file source
                        
  -a ALBACORE_SUMMARY_SOURCE, --albacore-summary-source ALBACORE_SUMMARY_SOURCE Albacore summary source
  
  -q FASTQ_SOURCE, --fastq-source FASTQ_SOURCE
                        fastq file source
                        
  -o OUTPUT_DIRECTORY, --output OUTPUT_DIRECTORY
                        output directory
                        
  -s SAMPLE_SHEET_SOURCE, --sample-sheet-source SAMPLE_SHEET_SOURCE
                        Sample sheet source```
                       
 <a name="example"></a>
 * #### Example
 >>>
Example with optional arguments:

```$ python3 toulligqc.py --run_name FAF0256 --barcoding -config_file /path/to/configuration/file/```

Example with optional arguments but no config file:

```$ python3 toulligqc.py --run_name FAF0256 --barcoding --fast5_source /path/to/fast5/source \
-albacore_summary_source /path/to/albacore/summary/source --fastq_source /path/to/fastq/source\
 --output /path/to/output/directory --sample-sheet-source /path/to/sample/sheet```

<a name="configuration-file"></a>
### 2.2 Configuration file

A configuration file can be used, the required informations has to be defined as following in the same order :

```[config]

fast5.directory=/path/to/fast5/directory/ (containing either FAST5, FAST5.tar.gz or FAST5.tar.bz2 files)

albacore.summary.directory=/path/to/albacore/sequencing/summary/directory/or/file

result.directory =/path/to/result/directory/(directory where the results are stored)

fastq.directory=/path/to/fastq/directory/ (containing either FASTQ or FASTQ.bz2 files)

design.file=/path/to/sample/sheet

[extension]

fast5.file.extension=tar.bz2

fastq.file.extension=bz2```

In the config part, the sample.sheet directory can be omitted if barcodes were not used in the run.
In the extension part, you must indicate the extensions used for your files. The file types currently supported are:

for FAST5 file : tar.gz, tar.bz2, fast5

for FASTQ file : bz2, gz (?), fastq

<a name="sample-sheet-for-barcoded-samples"></a>
### 2.3 Sample sheet
 
A sample sheet is required if barcodes were used. The sample sheet file describes the different samples and their corresponding barcodes. The **Index column is mandatory** and  must contain the Oxford Nanopore Technology **barcode number**. For example **01, 02, 11**. The **other columns are optional** but can be useful to define your samples for the following analyses. They can be modified at your convenience.

design.csv example:

index | Reads | 
------- | ------- 
 2015341_BC01 | dnacpc14_20170328_FNFAF04250_MN17734_mux_scan_1D_validation_test1_45344_barcode01_template.fastq.bz2 

## 3.Output
We can see a report example in the git repository.


