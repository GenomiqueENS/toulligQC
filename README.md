toulligQC
==========
This program is dedicated to the QC analyses of Oxford Nanopore runs, barcoded or not.
It requires a design file describing the barcodes used if the run was barcoded.
It partly relies on log file produced during the basecalling process by the Oxford Nanopore basecaller, Albacore.
This program will produce a set of graphs and statistic files


# Table of content

 1.Python module

  * [basecalling_stat_plotter1D](#basecalling_stat_plotter1D)

  * [getter1D](#getter1D)

  * [fast5_data_extractor](#fast5_data_extractor)

2.[example](#example)

3.[requirements](#requirements)

Python module
==============

## 1)basecalling_stat_plotter1D


This module provides support for create a report as a word document and a pdf file in the images directory. Files containing statistics for each barcode is provided in the statistics directory.

Several python files are provided:

 #### class basecalling_stat_plotter1D
 
   returns a set of graphs and statistics for the creation of report as well as statistics for the statistics files.

     meanqscore_barcode
     Write the mean qscore extracted from the log file provided by albacore in the statistics files according to the barcodes

     date
      Returns the date of a Minion run from the log file provided by albacore in researching the date in the run id field

     stat_generation
        Generates a dictionary of statistics such as quartile, std for the creation of a log file like aozan from the          summary log provided by albacore

    barcode_pie_chart
        Plots the barcode pie chart from the selection of barcodes

    reads_size_selection_barcode
        Plots the histogram of reads size by bins of 100 for the barcode choosed from design file

    statistics_read_size
        Gets statistics on the fastq from the file containing the fastq bz2 files decompressed

    histogram_count_reads
        Plots the histogram of count of different types of reads: template, complement, full_2D from albacore log file

    quality_reads_boxplot
        Plots a boxplot of reads quality

    channel_count
        Plots an histogram of channel count

    read_time
        Plots an histogram of reads length

    minion_flowcell_layout
        Represents the layout of a minion flowcell

    plot_performance(pore_measure)
        Plots the pore performance in terms of reads per pore
        parameters: pore_measure: reads number per pore

    get_selection
       Returns the selection of barcode from the design file

    statistics_dataframe
        Presents statitstics retrieved from statistics files in the statistics directory for each barcode as a dataframe to make
        the reading easier.

    read_size_total
       Plots an histogram of reads size by bins of 100 for all barcodes


## 2)getter1D

This module provided informations about the minion runs and the fastq sequences

### get_MinknowVersion(h5py_file)
      parameter: fast5 file open with h5py
      Gets the Minknow version from fast5 file

### getFlowcellId(h5py_file)
      parameter: fast5 file open with h5py
      Gets the flowcell id from fast5 file

### get_Hostname(h5py_file)
      parameter: fast5 file open with h5py
      Gets the hostname from fast5 file

### getNumMinION(h5py_file)
      parameter: fast5 file open with h5py
      Gets the number of Minion run

### getProtocolRunId(h5py_file)
      parameter: fast5 file open with h5py
      Gets the run id protocol from fast5 file

### get_barcode()
      parameter: fast5 file open with h5py
      Gets the barcode from a file given in input

### get_fastq(selection)
      parameter: fast5 file open with h5py
      Gets the fastq sequence

## 3)fast5_data_extractor

### fast5_data_extractor(fast5_file_directory)
      Creates a dataframe from collections of fast5 files
      param fast5_file_directory: directory where fast5 files are stored
      return : tuple with different informations about fast5 files

### write_data(tuple_array)
      parameter: tuple array containing informations about the fast5 files
      Writes data in a tsv file

### read_data(data_file)
      Reads the file created previously

example
==========

First of all a set of files are required before the python scripts run:

        * a design file named design.csv which describes the different sample barcoded. It's only the first column which is important. The rest of files may be modified at your convenience. An example might be:

  index | Reads | Description | Date | FastqFormat | RepTechGrou
------- | ------- | ------------- | -------- | -------------- | ---------------
 2015341_BC01 | dnacpc14_20170328_FNFAF04250_MN17734_mux_scan_1D_validation_test1_45344_barcode01_template.fastq.bz2 |  WT1_BC01 | 2017-01-24 | fastq-sanger | WT1_BC01

First Header | Second Header
------------ | -------------
Content from cell 1 | Content from cell 2
Content in the first column | Content in the second column

        An important thing is that the barcodes must be written with BC followed by two digits (BC01, BC02,....,BC80)

        * a configuration file: this file includes the path at your different files. These files numbered four in the following order:

                * the directory where the fast5 are holded
                * the file where the Albacore log is placed
                * the directory where the fastq are located in the form of bz2 files
                * the design file mentioned above

This file must be in the following form with these names:

[config]
fast5.directory=path to our fast5 directory
log.file= path to our log file
bz2.fastq.directory=path to fastq files in the form of bz2 files
design_file=path to design file

The first line is important for the python script because they use a parser in order to parse the configuration file and mustn't be modified

Once the files described above are ready an placed in the dcripts directory, we can launch the programm by means of the following command :
        * python3 main.py name_sample * where name sample represents the Minion run's sample name.


Afterwards two directories are created:

   * statistics containing statistics for each barcode
   * images containing the images generated by the script. It will serve to create the report.

Then two files are also created:

   *  a report .docx containing a set of graphs generated by the script
                                                                                                                                1,1          Haut

requirements
===============

Some modules needed be installed in order to the script is executed:

* matplotlib
* h5py
* pandas
* seaborn
* numpy
* PyPDF2
* csv
* python-docx
* biopython

This script is written in python3 and not in python2.

A docker file was created containing these modules with the correct version of python.

A config file was created which must be modified according to the configuration of its system.
