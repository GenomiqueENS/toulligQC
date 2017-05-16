toulligQC
==========


# Table of content

1.[basecalling_stat_plotter1D](#basecalling_stat_plotter1D)

2.[getter1D](#getter1D)

3.[fast5_data_extractor](#fast5_data_extractor)

4.[example](#example)

4.[requirements](#requirements)

basecalling_stat_plotter1D
===========================

This module provides support for create a report as a word document and a pdf file in the images directory. Files containing statistics for each barcode is provided in the statistics directory.

Several python files are provided:

 class basecalling_stat_plotter1D
   return a set of graphs and statistics for the creation of report as well as statistics for the statistics files.

    meanqscore_barcode
     Write the mean qscore extracted from the log file provided by albacore in the statistics files according to the barcodes

    date
      Return the date of a Minion run from the log file provided by albacore in researching the date in the run id field

    stat_generation
        Generate a dictionary of statistics such as quartile, std for the creation of a log file like aozan from the          summary log provided by albacore

    barcode_pie_chart
        Plot the barcode pie chart from the selection of barcodes

    reads_size_selection_barcode
        Plot the histogram of reads size by bins of 100 for the barcode choosed from design file

    statistics_read_size
        Get statistics on the fastq from the file containing the fastq bz2 files decompressed

    histogram_count_reads
        Plot the histogram of count of different types of reads: template, complement, full_2D from albacore log file

    quality_reads_boxplot
        Plot a boxplot of reads quality

    channel_count
        Plot an histogram of channel count

    read_time
        Plot an histogram of reads length

    minion_flowcell_layout
        Represent the layout of a minion flowcell

    plot_performance(pore_measure)
        Plot the pore performance in terms of reads per pore
        parameters: pore_measure: reads number per pore

    get_selection
       Return the selection of barcode from the design file

    statistics_dataframe
        Presents statitstics retrieved from statistics files in the statistics directory for each barcode as a dataframe to make
        the reading easier.

    read_size_total
       Plot an histogram of reads size by bins of 100 for all barcodes


getter1D
=========

This module provided informations about the minion runs and the fastq sequences

### get_MinknowVersion(h5py_file)
      parameter: fast5 file open with h5py
      Get the Minknow version from fast5 file

### getFlowcellId(h5py_file)
      parameter: fast5 file open with h5py
      Get the flowcell id from fast5 file

### get_Hostname(h5py_file)
      parameter: fast5 file open with h5py
      Get the hostname from fast5 file

### getNumMinION(h5py_file)
      parameter: fast5 file open with h5py
      Get the number of Minion run

### getProtocolRunId(h5py_file)
      parameter: fast5 file open with h5py
      Get the run id protocol from fast5 file

### get_barcode()
      parameter: fast5 file open with h5py
      Get the barcode from a file given in input

### get_fastq(selection)
      parameter: fast5 file open with h5py
      Get the fastq sequence

fast5_data_extractor
=====================

### fast5_data_extractor(fast5_file_directory)
      Create a dataframe from collections of fast5 files
      param fast5_file_directory: directory where fast5 files are stored
      return : tuple with different informations about fast5 files

### write_data(tuple_array)
      parameter: tuple array containing informations about the fast5 files
      Write data in a tsv file

### read_data(data_file)
      Read the file created previously

example
==========

A configuration file must be created named config2.txt(different for platform, see below) and not an other name or modify the python scripts.

The config file must be modify according to your system. It's very important.
This file must be organized as following:

*[config]*

*fast5.directory=/home/ferrato/shares-net/sequencages/nanopore/fast5/raw*

*log.file=/home/ferrato/shares-net/sequencages/nanopore/albacore-logs*

*bz2.fastq.directory=/home/ferrato/shares-net/sequencages/nanopore/fastq*

*design_file=/home/ferrato/ownCloud/fast5_1D*

*script_directory=$HOME/ownCloud/fast5_1D*

The [config] must not be modified because it's used by the scripts python. Always at the top of config file.

For the platform:
* the config file created must be specified in the shell script 
* go to the shell scripts and modify the file and file 2 variables at your liking but the file2 must always end by config2.txt
* don't put / to the end of each line for the config file even if it's a directory
* a file named config2.txt is created. Don't modify it.

The program runs with the following command:
  python3 main.py name_sample

name_sample: sample name of run (ex.20170328_FAF04250)

Then a path at fast5 directory is asked: if you use Docker on the platform ibens you must put /fast5_directory which corresponds to /home/ferrato/shares-net/sequencages/nanopore/fast5/raw path for my config file then the remainder.

Afterwards two directories are created:

* statistics containing statistics for each barcode
* images containing the images generated by the script. It will serve to create the report.

Then three files are also created:

* a report.docx containing a set of graphs generated by the script
* a pdf file which contains the different graph in the case of some figures will be unreadable or not present in the docx

A file named layout.pdf is necessary. That file represent the layout of a minion flowcell.
It must be in the directory where the script is executed. This latter is present in the git directory.


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
