# -*- coding: utf-8 -*-

#                  ToulligQC development code
#
# This code may be freely distributed and modified under the
# terms of the GNU General Public License version 3 or later
# and CeCILL. This should be distributed with the code. If you
# do not have a copy, see:
#
#      http://www.gnu.org/licenses/gpl-3.0-standalone.html
#      http://www.cecill.info/licences/Licence_CeCILL_V2-en.html
#
# Copyright for this code is held jointly by the Genomic platform
# of the Institut de Biologie de l'École Normale Supérieure and
# the individual authors.


#Graph generation

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pandas.tools.plotting import table
import re

def _make_table(value, ax, metric_suppression=''):
    '''
    Creation of a statistics table printed with the graph
    :param value: information measured
    :param ax: axes used
    :param metric_suppression: suppression of a metric when we use the describe pandas function
    '''

    if metric_suppression:
        the_table = table(ax, np.round(value.describe().drop(metric_suppression), 2), loc='center', )
    else:
        the_table = table(ax, np.round(value.describe(), 2), loc='center', )

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.axis('off')

    the_table.set_fontsize(12)
    the_table.scale(1, 1.2)


def _safe_log(x):
    '''
    Verification that we haven't a null value
    :param x: tested value
    :return: log2 value or 0
    '''
    if x <= 0:
        return 0
    return np.log2(x)

def read_length_histogram(albacore_log, my_dpi,result_directory):
    """
    Plots an histogram of the reads length by bins of 100 for each of the barcodes described in the design file or without barcode
    """
    output_file = result_directory + '/read_length_histogram.png'
    sequence_length_template = albacore_log['sequence_length_template']
    minimum, maximum = min(sequence_length_template), max(sequence_length_template)

    fig = plt.figure(figsize=(14, 8), dpi=my_dpi)
    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])

    n, bins, patches = ax.hist(sequence_length_template, edgecolor='black',
                               bins=2 ** np.linspace(_safe_log(minimum), _safe_log(maximum), 30))
    ax.set_xscale('log', basex=2)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
    ax.set_xticks(bins)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_rotation('vertical')
    ax.set_xlabel('Read size(in pb)')
    ax.set_ylabel('Read number')
    ax.set_title('Read size histogram')

    ax2 = plt.subplot(gs[1])

    the_table = table(ax2, np.round(sequence_length_template.describe(), 2), loc='center')

    ax2.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)
    ax2.axis('off')

    the_table.set_fontsize(12)
    the_table.scale(1, 1)

    plt.savefig(output_file)
    plt.close()

    return 'Read size histogram', output_file

def phred_score_frequency(albacore_log, my_dpi, result_directory):
    '''
    Plot the distribution of the phred score
    '''
    output_file = result_directory + '/phred_score_frequency.png'
    plt.figure(figsize=(800 / my_dpi, 8), dpi=my_dpi)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])

    sns.distplot(albacore_log['mean_qscore_template'], bins=15, color='green',
                 hist_kws=dict(edgecolor="k", linewidth=1))
    plt.xlabel("mean_qscore")
    plt.ylabel("Frequency")
    plt.title("Phred score frequency")
    rd = albacore_log['mean_qscore_template'].describe().drop('count').round(2).reset_index()

    plt.axvline(x=albacore_log['mean_qscore_template'].describe()['50%'], color='green')

    ax2 = plt.subplot(gs[1])
    _make_table(rd, ax2, 'count')
    plt.savefig(output_file)
    plt.close()

    return "Phred score frequency", output_file

def scatterplot(albacore_log, my_dpi, result_directory):
    '''
    Plot the scatter plot representing the relation between the phred score and the sequence length
    '''
    output_file = result_directory + '/scatter_plot.png'
    plt.figure(figsize=(1200 / my_dpi, 8), dpi=my_dpi)
    plt.scatter(x=albacore_log['sequence_length_template'], y=albacore_log['mean_qscore_template'])
    plt.xlim(0, 100000)
    plt.xlabel("sequence_length_template")
    plt.ylabel("mean_qscore_template")
    plt.title("mean template qscore function of template read length")
    plt.savefig(output_file)
    plt.close()

    return "mean template qscore function of template read length", output_file

def read_count_histogram(albacore_log, my_dpi, result_directory):
    """
    Plots the count histograms of count  of the different types of reads eventually available in a Minion run: template, complement, full_2D.
    """
    output_file = result_directory + '/read_count_histogram.png'
    plt.figure(figsize=(1200 / my_dpi, 8), dpi=my_dpi)
    fast5_raw = len(albacore_log['num_events'])
    fast5_template_basecalled = \
        len(albacore_log['num_events']) - len(albacore_log[albacore_log['num_called_template'] == 0])

    read_type = [fast5_raw, fast5_template_basecalled]
    label = ("fast5_raw", "fast5_template_basecalled")
    nd = np.arange(len(read_type))

    bars = plt.bar(nd, read_type, align='center', color=["lightblue", "salmon"])
    plt.xticks(nd, label)
    plt.xlabel("read type")
    plt.ylabel("Counts")
    plt.title("Counts of read template")

    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2., 1 * height, '%d' % int(height), ha='center', va='bottom')

    plt.savefig(output_file)
    plt.close()

    return "Counts of read template", output_file

def read_quality_boxplot(albacore_log, my_dpi, result_directory):
    """
    Plots a boxplot of reads quality
    """
    output_file = result_directory + '/read_quality_boxplot.png'
    plt.figure(figsize=(1200 / my_dpi, 8), dpi=my_dpi)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])

    dataframe = albacore_log.loc[:, ["mean_qscore_template"]]
    sns.boxplot(data=dataframe)
    plt.title('Boxplot of read quality')
    plt.ylabel('Phred score')

    ax2 = plt.subplot(gs[1])
    _make_table(dataframe, ax2, 'count')

    plt.savefig(output_file)
    plt.close()

    return 'Boxplot of read quality', output_file

def channel_count_histogram(albacore_log, my_dpi, result_directory):
    """
    Plots an histogram of the channel count according to the channel number
    """
    output_file = result_directory + '/channel_count_histogram.png'
    plt.figure(figsize=(1200 / my_dpi, 8), dpi=my_dpi)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])
    ax.hist(albacore_log['channel'], edgecolor='black',
            bins=range(min(albacore_log['channel']), max(albacore_log['channel']) + 64, 64))
    ax.set_xlabel("Channel number")
    ax.set_ylabel("Count")
    ax.set_title("Channel counts")

    channel_count = albacore_log['channel']
    total_number_reads_per_channel = pd.value_counts(channel_count)
    ax2 = plt.subplot(gs[1])
    _make_table(total_number_reads_per_channel, ax2, metric_suppression=['mean', 'std', '50%', '75%', '25%'])

    plt.savefig(output_file)
    plt.close()

    return "Channel counts", output_file


def read_number_run(albacore_log, my_dpi, result_directory):
    """
    Plots the reads produced along the run against the time(in hour)
    """
    output_file = result_directory + '/read_number_run.png'
    plt.figure(figsize=(1200 / my_dpi, 8), dpi=my_dpi)
    start_time = albacore_log["start_time"] / 3600
    start_time_sorted = sorted(start_time)
    plt.scatter(start_time_sorted, np.arange(len(start_time_sorted)))
    plt.ylabel("produced reads")
    plt.xlabel("hour")
    plt.title("Read produced along the run")

    plt.savefig(output_file)
    plt.close()

    return "Read produced along the run", output_file

def _minion_flowcell_layout():
    """
    Represents the layout of a minion flowcell
    """
    seeds = [125, 121, 117, 113, 109, 105, 101, 97,
             93, 89, 85, 81, 77, 73, 69, 65,
             61, 57, 53, 49, 45, 41, 37, 33,
             29, 25, 21, 17, 13, 9, 5, 1]

    flowcell_layout = []
    for s in seeds:
        for block in range(4):
            for row in range(4):
                flowcell_layout.append(s + 128 * block + row)
    return flowcell_layout


def plot_performance(pore_measure, my_dpi, result_directory):
    """
    Plots the channels occupancy by the reads
    @:param pore_measure: reads number per pore
    """
    output_file = result_directory + '/channel_occupancy.png'
    flowcell_layout = _minion_flowcell_layout()

    pore_values = []
    for pore in flowcell_layout:
        if pore in pore_measure:
            pore_values.append(pore_measure[pore])
        else:
            pore_values.append(0)

    d = {'rownum': list(range(1, 17)) * 32,
         'colnum': sorted(list(range(1, 33)) * 16),
         'tot_reads': pore_values,
         'labels': flowcell_layout}

    df = pd.DataFrame(d)

    d = df.pivot("rownum", "colnum", "tot_reads")
    d2 = df.pivot("rownum", "colnum", "labels")
    plt.figure(figsize=(1200 / my_dpi, 8), dpi=my_dpi)
    sns.heatmap(d, fmt="", annot=d2, linewidths=.5, cmap="YlGnBu", annot_kws={"size": 7})
    plt.title('Channel occupancy')

    plt.savefig(output_file)
    plt.close()

    return 'Channel occupancy', output_file

def barcode_percentage_pie_chart(albacore_log, barcode_selection, my_dpi, result_directory):
    """
    Plots a pie chart of the barcode percentage of a run. Needs the design file describing the barcodes to run
    """
    output_file = result_directory + '/barcode_percentage_pie_chart.png'
    plt.figure(figsize=(800 / my_dpi, 800 / my_dpi), dpi=my_dpi)
    for element in barcode_selection:

        if all(albacore_log['barcode_arrangement'] != element):
            print("The barcode {} doesn't exist".format(element))
            return False

    barcode = albacore_log['barcode_arrangement']
    barcode_count = barcode.value_counts()
    count_sorted = barcode_count.sort_index()[barcode_selection]
    total = sum(count_sorted)

    cs = cm.Paired(np.arange(len(barcode_selection)) / len(barcode_selection))

    sizes = [(100 * chiffre) / total for chiffre in count_sorted.values]
    if len(barcode_selection) <= 10:
        fig1, ax1 = plt.subplots()
        ax1.pie(sizes, labels=barcode_selection, autopct='%.2f%%', startangle=90, colors=cs)
        ax1.axis('equal')

    else:
        fig = plt.figure(figsize=(20, 10))
        ax1 = fig.add_subplot(111)
        length = np.arange(0, len(barcode_count))
        ax1.set_title('Percentage of different barcodes')
        ax1.bar(length, barcode_count, color=cs)
        ax1.set_xticks(length)
        ax1.set_xticklabels(barcode_selection)

    plt.savefig(output_file)
    plt.close()

    return 'Percentage of different barcodes', output_file

def barcode_length_boxplot(albacore_log, barcode_selection, my_dpi, result_directory):
    '''
    Plot the length boxplot for each barcode indicated in the sample sheet
    '''
    output_file = result_directory + '/barcode_length_boxplot.png'
    pattern = '(\d{2})'
    dico = {}


    fig = plt.figure(figsize=(1200 / my_dpi, 8), dpi=my_dpi)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])


    for barcode in barcode_selection:
        barcode_selected_dataframe = albacore_log[albacore_log['barcode_arrangement'] == barcode]
        match = re.search(pattern, barcode)
        if match:
            dico[match.group(0)] = barcode_selected_dataframe['sequence_length_template']
        else:
            dico['Unclassified'] = barcode_selected_dataframe['sequence_length_template']
    barcode_selection_sequence_length_dataframe = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dico.items()]))
    sns.boxplot(data=barcode_selection_sequence_length_dataframe, showfliers=False)
    plt.xlabel('Barcodes')
    plt.ylabel('Read size(in pb)')
    plt.title('Read size distribution for each barcode')

    ax2 = plt.subplot(gs[1])
    _make_table(barcode_selection_sequence_length_dataframe,ax2)
    plt.savefig(output_file)
    plt.close()

    return 'Read size distribution for each barcode', output_file

def barcoded_phred_score_frequency(albacore_log, barcode_selection, my_dpi, result_directory):
    '''
    Plot the phred score distribution boxplot for each barcode indicated in the sample sheet
    '''
    output_file = result_directory + '/barcode_phred_score_boxplot.png'
    dico = {}
    pattern = '(\d{2})'

    fig = plt.figure(figsize=(1200 / my_dpi, 8), dpi=my_dpi)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])

    for barcode in barcode_selection:
        barcode_selected_dataframe = albacore_log[albacore_log['barcode_arrangement'] == barcode]
        match = re.search(pattern, barcode)
        if match:
            dico[match.group(0)] = barcode_selected_dataframe['mean_qscore_template']
        else:
            dico['Unclassified'] = barcode_selected_dataframe['mean_qscore_template']
    barcode_selection_phred_scrore_dataframe = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dico.items()]))
    sns.boxplot(data= barcode_selection_phred_scrore_dataframe, showfliers=False)
    plt.xlabel('Barcodes')
    plt.ylabel('Phred score')
    plt.title('Phred score distribution for each barcode')

    ax2 = plt.subplot(gs[1])
    _make_table(barcode_selection_phred_scrore_dataframe,ax2)
    plt.savefig(output_file)
    plt.close()

    return 'Phred score distribution for each barcode', output_file
