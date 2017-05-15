===========================
basecalling_stat_plotter1D
===========================

This module provides support for create a report as a word document and a pdf file in the images directory. Files containing statistics for each barcode is provided in the statistics directory.

Several python files are provided:

 **class basecalling_stat_plotter1D**
   return a set of graphs and statistics for the creation of report as well as statistics for the statistics files.

    *meanqscore_barcode*
     Write the mean qscore extracted from the log file provided by albacore in the statistics files according to the barcodes

    *date*
      Return the date of a Minion run from the log file provided by albacore in researching the date in the run id field

    *stat_generation*
        Generate a dictionary of statistics such as quartile, std for the creation of a log file like aozan from the summary log
        provided by albacore

    *barcode_pie_chart*
        Plot the barcode pie chart from the selection of barcodes

    *reads_size_selection_barcode*
        Plot the histogram of reads size by bins of 100 for the barcode choosed from design file

    *statistics_read_size*
        Get statistics on the fastq from the file containing the fastq bz2 files decompressed

    *histogram_count_reads*
        Plot the histogram of count of different types of reads: template, complement, full_2D from albacore log file

    *quality_reads_boxplot*
        Plot a boxplot of reads quality

    *channel_count*
        Plot an histogram of channel count

    *read_time*
        Plot an histogram of reads length


    *minion_flowcell_layout*
        Represent the layout of a minion flowcell

    *plot_performance(self, pore_measure)*
        Plot the pore performance in terms of reads per pore
        parameters: pore_measure: reads number per pore

    *get_selection*:
      Return the selection of barcode from the design file


    *statistics_dataframe*
        Presents statitstics retrieved from statistics files in the statistics directory for each barcode as a dataframe to make
        the reading easier.

    *read_size_total*
       Plot an histogram of reads size by bins of 100 for all barcodes
