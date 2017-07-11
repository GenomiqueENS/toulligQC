ToulligQC
=========
This program is dedicated to the QC analyses of Oxford Nanopore runs, barcoded or not. It requires a design file describing the barcodes used if the run was barcoded. It partly relies on log file produced during the basecalling process by the Oxford Nanopore basecaller, Albacore. This program will produce a set of graphs and statistic files in the form of pdf, html and docx report.
ToulligQC accept different formats: bz2, tar.bz2, fastq and fast5
First of all a set of files are required before the programm runs:
-a configuration file that must be in the following form :
[config] 
fast5.directory=path to our fast5 directory 
log.file= path to our log file
fastq.directory=path to fastq files in the form of bz2 files
design_file_directory=path to design file directory

[extension]
fast5.file.extension=tar.bz2
fastq.file.extension=bz2 

Warnings : The fast5 and fastq files must be present just after the directory indicated in the configuration file. The directory indicated must not contained another directory containing the files.

All that is before the equal sign must be conserved along with [config] and [extension]. For the rest you can indicate the path toward your files. Only the design_file_directory can be ommited if you didn’t use the barcodes.
In extension you must indicate the extension used for your file. For the moment the type of extension are supported are:
for fast5 file : tar.gz, tar.bz2, fast5
for fastq file : bz2, fastq
- a design_file if you used the barcodes. This one must be named design.csv and describes the different sample barcoded. Only the first column is important because it must contain the barcodes number that you used in the form of BC followed by two digits. For example BC01, BC11. The rest of files may be modified at your convenience.

An example is provided thereafter.

Installation
=============
## Option1 : Installation using Docker

ToulligQC and its dependancies are available throw Docker images. 
You can use a Docker image with Aozan and all its optional dependencies 
To see how install docker on your system, go to the Docker website. Even if Docker can run in virtual machines in Windows or macOS, we recommand to only run ToulligQC on a Linux host.
You can use a Docker image with ToulligQC and all its optional dependencies  instead of installating manually ToulligQC. This image is named genomicpariscentre/aozan:2.0. When you use this Docker image you need to mount all the required directories by Aozan in the Docker container.
A shell script called read_file.sh is provided which make mounting automatic with the configuration file.
Ce script prend les noms de fichiers indiqués dans le fichier de configuration et les traduit en leur vrai nom dans le cas où l’on utilise des liens symboliques puis il monte ces mêmes fichiers dans le contenair de Docker.
 
## Option2 : Local installation 
This option is also suitable if you are interested in further developping the package, but requires a little bit more hands-on.
Clone the repository locally 
git clone https://github.com/GenomicParisCentre/toulligQC.git
Install all dependencies indicated below.

Requirements:
=============
To run ToulligQC, you need to install the following software:
matplotlib 
h5py 
pandas 
seaborn 
numpy 
PyPDF2 
csv 
python-docx

On Debian/Ubuntu, you can install requirements  using the 'apt-get' command, here is an example:
$sudo apt-get install matplotlib

Organisation of your directory
===============================
The directory where the files are presented must be in the following form :
for fast5 directory the fast5 file must be named with the run name given in the argument line. For example FAF042450.fast5 or FAF04250.tar.gz for the run name argument FAF04250
for fastq file we must have a directory after the fastq directory named with the same run name that in the fast5 file above. This one is essential for the using of barcode because we have a fastq file for each barcode. We can have for example ten files in the fastq directory. 

Launching ToulligQC
=========================

A set of option are available :
python3 main.py -h
optional arguments:
  -h, --help          	Show this help message and exit
  -n RUN_NAME 	Run name
  -d                          	Docker usage
  -b                        	Barcode usage
  -s ARG [ARG ...], --arg ARG [ARG …]   Selected files

The run name correspond to this is indicated before the extension file for fastq and fast5 files.
For example if you have FAF2056.tar.bz2 the run name correspond to FAF0256 and not FAF0256.tar.bz2.
You run ToulligQC as follows :
python3 main.py -n run name with n argument mandatory.

Usage example
========================

Here I will provide a complete example.
Our run calls FAF04250.
The fast5 files and fastq files calls FAF04250.tar.bz2 and FAF04250.fastq.bz2.
First of all, we need a configuration file :
[config] 
fast5.directory=path to our fast5 directory containing the fast5 files
log.file= path to our log file containing the fastq files
fastq.directory=path to fastq files in the form of bz2 files
design_file_directory=path to design file directory

[extension]
fast5.file.extension=tar.bz2
fastq.file.extension=bz2 

The design file directory being optional if you don’t use the barcodes.
Then a design file if we use the barcodes named imperatively design.csv which describes the different sample barcoded. It's only the first column which is important. The rest of files may be modified at your convenience. An example might be:

index | Reads | Description | Date | FastqFormat | RepTechGrou
------- | ------- | ------------- | -------- | -------------- | ---------------
 2015341_BC01 | dnacpc14_20170328_FNFAF04250_MN17734_mux_scan_1D_validation_test1_45344_barcode01_template.fastq.bz2 |  WT1_BC01 | 2017-01-24 | fastq-sanger | WT1_BC01


After you must modify the read_file.sh script with the path toward your config file for the two first line.
 sed '/^$/d' /import/config.txt > /import/conf.txt
config_file=/import/conf.txt
For these two lines you must indicate the path toward your config file for the first line before the >     
token. Then you must indicate the path where your configuration file is with another name than your configuration file. 

After that, we launch the script.
Without barcode and  without docker
python3 main.py -n 20170104FAF04250
With barcodes :
python3 main.py -n 20170104FAF04250 -b
With docker :
python3 main.py -n  20170104FAF04250 -d
With docker and barcodes :
python3 main.py -n  20170104FAF04250 -b  -d

Here the run name is  20170104FAF04250. So just after the fastq directory another directory called  20170104FAF04250 must be created which will contain the fastq files with extension indicated in the config file without the point at the beginning(bz2 and not .bz2)
For the fast5 file ……..(to see).
The program generates a set of graphs and statistics. More precisely we have got 8 graphs or 7 graphs without barcode. Moreover either a global statistics file is yielded or a statistic file by barcode if barcodes are used.
ToulligQC yield a report in the form of a html, pdf or docx file.
