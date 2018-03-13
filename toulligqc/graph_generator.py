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


# Graph generation

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

pd.options.display.float_format = '{:,}'.format


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

def _make_table_html(value):
    '''
    Creation of a statistics table printed with the graph in report.html
    :param value: information measured
    :param metric_suppression: suppression of a metric when we use the describe pandas function
    '''

    table=np.round(value.describe(),2)
    pd.options.display.float_format = '{:,}'.format
    table_html = pd.DataFrame.to_html(table)
    return table_html

def _safe_log(x):
    '''
    Verification that we haven't a null value
    :param x: tested value
    :return: log2 value or 0
    '''
    if x <= 0:
        return 0
    return np.log2(x)


#  1D plots

def read_count_histogram(albacore_log, main, my_dpi, result_directory, desc):
    """
    Plots the count histograms of count  of the different types of reads eventually available in a Minion run: template, complement, full_2D.
    """
    output_file = result_directory + '/' + main + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    #ax = plt.subplot(gs[0])

    fast5_raw = len(albacore_log['num_events'])
    fast5_template_basecalled = len(albacore_log[albacore_log['num_called_template'] != 0])
    read_pass = len(albacore_log[albacore_log['passes_filtering'] == True])
    read_fail = len(albacore_log[albacore_log['passes_filtering'] == False])
    read_type = [fast5_raw, fast5_template_basecalled, read_pass, read_fail]
    label = ("fast5", "1D", "1D pass", "1D fail")
    nd = np.arange(len(read_type))

    bars = plt.bar(nd, read_type, align='center', color=["lightblue", "salmon", "yellowgreen", "orangered"])
    plt.xticks(nd, label)
    plt.xlabel("Read type")
    plt.ylabel("Counts")
    plt.title(main)

    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2., 1 * height, '%d' % int(height), ha='center', va='bottom')
    table_html=None

    plt.savefig(output_file)
    plt.close()

    return main, output_file, table_html, desc

def read_length_multihistogram(albacore_log, main, my_dpi, result_directory, desc):
    """
    Plots an histogram of the reads length by bins of 100 for each of the barcodes described in the design file or without barcode
    """
    output_file = result_directory + '/' + main + '.png'
    read_pass = albacore_log.sequence_length_template.loc[True == albacore_log['passes_filtering']]
    read_fail = albacore_log.sequence_length_template.loc[False == albacore_log['passes_filtering']]
    sequence_length_template = albacore_log.sequence_length_template[albacore_log['num_called_template'] != 0]
    minimum, maximum = min(sequence_length_template), max(sequence_length_template)
    read_type=['1D','1D pass','1D fail']

    fig = plt.figure(figsize=(12, 7), dpi=my_dpi)
    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])
    plt.subplots_adjust(bottom=0.015, top=1.0)

    n, bins, patches = ax.hist([sequence_length_template,read_pass,read_fail],color=["salmon", "yellowgreen", "orangered"], edgecolor='black', label=read_type,
                               bins=2 ** np.linspace(_safe_log(minimum), _safe_log(maximum), 30))
    plt.legend()

    ax.set_xscale('log', basex=2)
    ax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.1f}'))
    ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))
    ax.set_xticks(bins)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_rotation('vertical')
    ax.set_xlabel('Read length(bp)')
    ax.set_ylabel('Read number')
    ax.set_title(main)

    # ax2 = plt.subplot(gs[1])

    dataframe = pd.DataFrame({"1D": sequence_length_template, "1D pass": read_pass, "1D fail": read_fail})
    dataframe = dataframe[["1D","1D pass","1D fail"]]

    plt.savefig(output_file)
    plt.close()

    table_html = _make_table_html(dataframe)

    return main, output_file, table_html, desc

def allread_number_run(albacore_log, main, my_dpi, result_directory, desc):
    """
    Plots the reads produced along the run against the time(in hour)
    """
    output_file = result_directory + '/' + main + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi)

    start_time_sorted = sorted(sorted(albacore_log['start_time']/3600))
    read_pass = sorted(albacore_log.start_time.loc[True == albacore_log['passes_filtering']]/3600)
    read_fail = sorted(albacore_log.start_time.loc[False == albacore_log['passes_filtering']]/3600)
    #read_type = ['1D','1D pass', '1D fail']

    plt_read_1d=plt.scatter(start_time_sorted, np.arange(len(start_time_sorted)),color='salmon')
    plt_read_pass=plt.scatter(read_pass, np.arange(len(read_pass)),color='yellowgreen')
    plt_read_fail=plt.scatter(read_fail, np.arange(len(read_fail)),color='orangered')
    plt.ylabel("Read number")
    plt.xlabel("Time (Hour)")
    plt.legend((plt_read_1d,plt_read_pass,plt_read_fail),('1D','1D pass','1D fail'))
    plt.title(main)

    plt.savefig(output_file)
    plt.close()

    table_html = None

    return main, output_file, table_html, desc

def read_quality_multiboxplot(albacore_log, main, my_dpi, result_directory, desc):
    """
    Plots a boxplot of reads quality
    """
    output_file = result_directory + '/' + main + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    mean_qscore = albacore_log.loc[:, "mean_qscore_template"]
    read_pass = albacore_log.mean_qscore_template.loc[True == albacore_log['passes_filtering']]
    read_fail = albacore_log.mean_qscore_template.loc[False == albacore_log['passes_filtering']]

    my_pal = {"1D": "salmon", "1D pass": "yellowgreen", "1D fail": "orangered"}
    order=["1D","1D pass","1D fail"]
    dataframe = pd.DataFrame({"1D": mean_qscore, "1D pass": read_pass, "1D fail": read_fail })
    sns.boxplot(data=dataframe, ax=plt.subplot(gs[0]),palette=my_pal,order=order)

    plt.subplots_adjust(bottom=0.015, top=1.0)
    plt.ylabel('Mean Phred score')
    plt.title(main)
    plt.legend()

    #arr = dataframe.as_matrix()
    #print(arr)

    #ax2 = plt.subplot(gs[1])
    #_make_table(dataframe, ax2, 'count')

    dataframe = dataframe[["1D","1D pass","1D fail"]]
    table_html = _make_table_html(dataframe)


    plt.savefig(output_file)
    plt.close()
    return main, output_file, table_html, desc

def phred_score_frequency(albacore_log, main, my_dpi, result_directory, desc):
    '''
    Plot the distribution of the phred score
    '''
    output_file = result_directory + '/' + main + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    plt.subplots_adjust(bottom=0.015, top=1.0)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])
    plt.subplots_adjust(bottom=0.015, top=1.0)

    mean_qscore = albacore_log['mean_qscore_template']

    sns.distplot(mean_qscore, bins=15, color='salmon',
                 hist_kws=dict(edgecolor="k", linewidth=1),hist=True,label="1D")

    plt.legend()
    plt.xlabel("Mean Phred score")
    plt.ylabel("Frequency")
    plt.title(main)

    dataframe = pd.DataFrame({"1D": mean_qscore})
    rd = dataframe.describe().drop('count').round(2).reset_index()

    plt.axvline(x=mean_qscore.describe()['50%'], color='salmon')

    # ax2 = plt.subplot(gs[1])
    # _make_table(rd, ax2, 'count')
    plt.savefig(output_file)
    plt.close()

    table_html = _make_table_html(rd)

    return main, output_file, table_html, desc

def allphred_score_frequency(albacore_log, main, my_dpi, result_directory, desc):
    '''
    Plot the distribution of the phred score
    '''
    output_file = result_directory + '/' + main + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])
    plt.subplots_adjust(bottom=0.015, top=1.0)
    mean_qscore = albacore_log['mean_qscore_template']
    mean_qscore_read_pass = albacore_log.mean_qscore_template.loc[True == albacore_log['passes_filtering']]
    mean_qscore_read_fail = albacore_log.mean_qscore_template.loc[False == albacore_log['passes_filtering']]

    sns.distplot(mean_qscore_read_pass, bins=15,hist_kws=dict(edgecolor="k", linewidth=1) ,color='yellowgreen',
                 hist=True, label='1D pass')
    sns.distplot(mean_qscore_read_fail, bins=15,hist_kws=dict(edgecolor="k", linewidth=1) ,color='orangered', label='1D fail',
                 hist=True)
    plt.legend()
    plt.xlabel("Mean Phred score")
    plt.ylabel("Frequency")
    plt.title(main)

    dataframe = pd.DataFrame({"1D": mean_qscore, "1D pass": mean_qscore_read_pass, "1D fail": mean_qscore_read_fail})

    rd = dataframe.describe().drop('count').round(2).reset_index()

    plt.axvline(x=mean_qscore_read_pass.describe()['50%'], color='yellowgreen')
    plt.axvline(x=mean_qscore_read_fail.describe()['50%'], color='orangered')

    plt.savefig(output_file)
    plt.close()

    dataframe = dataframe[["1D","1D pass","1D fail"]]
    table_html = _make_table_html(dataframe)

    return main, output_file, table_html, desc

def all_scatterplot(albacore_log, main, my_dpi, result_directory, desc):
    '''
    Plot the scatter plot representing the relation between the phred score and the sequence length in log
    '''
    output_file = result_directory + '/' + main + '.png'
    #plt.figure(figsize=(1100/my_dpi, 6), dpi=my_dpi)
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    #plt.subplots_adjust(bottom=0.015, top=1.0)

    length_read_pass = albacore_log.sequence_length_template.loc[True == albacore_log['passes_filtering']]
    length_read_fail = albacore_log.sequence_length_template.loc[False == albacore_log['passes_filtering']]

    qscore_read_pass = albacore_log.mean_qscore_template.loc[True == albacore_log['passes_filtering']]
    qscore_read_fail = albacore_log.mean_qscore_template.loc[False == albacore_log['passes_filtering']]

    read_pass=plt.scatter(x=length_read_pass, y=qscore_read_pass,color="yellowgreen")
    read_fail=plt.scatter(x=length_read_fail, y=qscore_read_fail,color="orangered")
    plt.legend((read_pass,read_fail),("1D pass","1D fail"))
    plt.xlim(np.min(albacore_log.sequence_length_template.loc[albacore_log.sequence_length_template>0]),np.max(albacore_log.sequence_length_template))
    plt.xscale('log')
    plt.xlabel("Sequence length")
    plt.ylabel("Mean Phred score")
    plt.title(main)
    plt.savefig(output_file)
    plt.close()

    table_html = None

    return main, output_file, table_html, desc

def channel_count_histogram(albacore_log, main, my_dpi, result_directory, desc):
    """
    Plots an histogram of the channel count according to the channel number
    """
    output_file = result_directory + '/' + main + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])
    ax.hist(albacore_log['channel'], edgecolor='black',
            bins=range(min(albacore_log['channel']), max(albacore_log['channel']) + 64, 64))
    ax.set_xlabel("Channel number")
    ax.set_ylabel("Count")
    ax.set_title(main)

    channel_count = albacore_log['channel']
    total_number_reads_per_channel = pd.value_counts(channel_count)
    ax2 = plt.subplot(gs[1])
    _make_table(total_number_reads_per_channel, ax2, metric_suppression=['mean', 'std', '50%', '75%', '25%'])

    plt.savefig(output_file)
    plt.close()
    table = total_number_reads_per_channel.describe()
    table_html = pd.DataFrame.to_html(table)

    return main, output_file, table_html, desc

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


def plot_performance(pore_measure, main, my_dpi, result_directory, desc):
    """
    Plots the channels occupancy by the reads
    @:param pore_measure: reads number per pore
    """
    output_file = result_directory + '/' + main + '.png'
    flowcell_layout = _minion_flowcell_layout()

    pore_values = []
    for pore in flowcell_layout:
        if pore in pore_measure:
            pore_values.append(pore_measure[pore])
        else:
            pore_values.append(0)

    d = {'Row number': list(range(1, 17)) * 32,
         'Column number': sorted(list(range(1, 33)) * 16),
         'tot_reads': pore_values,
         'labels': flowcell_layout}

    df = pd.DataFrame(d)

    d = df.pivot("Row number", "Column number", "tot_reads")
    d2 = df.pivot("Row number", "Column number", "labels")
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    sns.heatmap(d, fmt="", linewidths=.5, cmap="YlGnBu", annot_kws={"size": 7}, cbar_kws={'label': 'Read number per pore channel',"orientation": "horizontal"})
    plt.title(main)

    plt.savefig(output_file)
    plt.close()

    table_html = None

    return main, output_file, table_html, desc

# For eache barcode 1D

def barcode_percentage_pie_chart_pass(albacore_log, main, barcode_selection, my_dpi, result_directory, desc):
    """
    Plots a pie chart of the barcode percentage of a run. Needs the design file describing the barcodes to run
    """
    output_file = result_directory + '/' + main + '.png'
    plt.figure(figsize=(800 / my_dpi, 800 / my_dpi), dpi=my_dpi)
    for element in barcode_selection:

        if all(albacore_log['barcode_arrangement'] != element):
            print("The barcode {} doesn't exist".format(element))
            return False

    barcode = albacore_log.barcode_arrangement.loc[True == albacore_log['passes_filtering']]
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
        ax1.set_title(main)
        ax1.bar(length, barcode_count, color=cs)
        ax1.set_xticks(length)
        ax1.set_xticklabels(barcode_selection)

    plt.savefig(output_file)
    plt.close()

    table_html = None

    return main, output_file, table_html, desc

def barcode_percentage_pie_chart_fail(albacore_log, main, barcode_selection, my_dpi, result_directory, desc):
    """
    Plots a pie chart of the barcode percentage of a run. Needs the design file describing the barcodes to run
    """
    output_file = result_directory + '/' + main + '.png'
    plt.figure(figsize=(800 / my_dpi, 800 / my_dpi), dpi=my_dpi)
    for element in barcode_selection:

        if all(albacore_log['barcode_arrangement'] != element):
            print("The barcode {} doesn't exist".format(element))
            return False

    barcode = albacore_log.barcode_arrangement.loc[False == albacore_log['passes_filtering']]
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
        ax1.set_title(main)
        ax1.bar(length, barcode_count, color=cs)
        ax1.set_xticks(length)
        ax1.set_xticklabels(barcode_selection)

    plt.savefig(output_file)
    plt.close()

    table_html = None

    return main, output_file, table_html, desc

def barcode_length_boxplot(albacore_log, main, barcode_selection, my_dpi, result_directory, desc):
    '''
    Plot the length boxplot for each barcode indicated in the sample sheet
    '''
    output_file = result_directory + '/' + main + '.png'
    pattern = '(\d{2})'
    dico = {}

    fig = plt.figure(figsize=(12, 7), dpi=my_dpi)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])
    plt.subplots_adjust(bottom=0.015, top=1.0)


    dico["passes_filtering"] = albacore_log['passes_filtering']
    for barcode in barcode_selection:
        barcode_selected_dataframe = albacore_log[albacore_log['barcode_arrangement'] == barcode]
        match = re.search(pattern, barcode)
        if match:
            dico[match.group(0)] = barcode_selected_dataframe['sequence_length_template']
        else:
            dico['Unclassified'] = barcode_selected_dataframe['sequence_length_template']
    barcode_selection_sequence_length_dataframe = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dico.items()]))
    dico.clear()
    melted = pd.melt(barcode_selection_sequence_length_dataframe, id_vars=["passes_filtering"],
                     var_name="barcodes", value_name="length")
    ax = sns.boxplot(data=melted,x='barcodes',y='length',hue='passes_filtering', showfliers=False,palette={True : "yellowgreen", False : "orangered"},hue_order=[True,False])
    handles, _ = ax.get_legend_handles_labels()
    plt.legend(bbox_to_anchor=(0.92, 1), loc=2, borderaxespad=0., labels=["Pass", "Fail"], handles=handles)
    plt.xlabel('Barcodes')
    plt.ylabel('Read length(bp)')
    plt.title(main)

    # ax2 = plt.subplot(gs[1])
    # _make_table(barcode_selection_sequence_length_dataframe, ax2)
    plt.savefig(output_file)
    plt.close()
    # table=np.round(barcode_selection_sequence_length_dataframe.describe(),2)
    # pd.options.display.float_format = '{:,}'.format
    # table_html = pd.DataFrame.to_html(table)
    table_html=_make_table_html(barcode_selection_sequence_length_dataframe)

    return main, output_file, table_html, desc

def barcoded_phred_score_frequency(albacore_log, main, barcode_selection, my_dpi, result_directory, desc):
    '''
    Plot the phred score distribution boxplot for each barcode indicated in the sample sheet
    '''
    output_file = result_directory + '/' + main + '.png'
    dico = {}
    pattern = '(\d{2})'

    fig = plt.figure(figsize=(12, 7), dpi=my_dpi)
    plt.subplots_adjust(bottom=0.015, top=1.0)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])

    dico["passes_filtering"] = albacore_log['passes_filtering']
    for barcode in barcode_selection:
        barcode_selected_dataframe = albacore_log[albacore_log['barcode_arrangement'] == barcode]
        match = re.search(pattern, barcode)
        if match:
            dico[match.group(0)] = barcode_selected_dataframe['mean_qscore_template']
        else:
            dico['Unclassified'] = barcode_selected_dataframe['mean_qscore_template']

    barcode_selection_phred_scrore_dataframe = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dico.items()]))
    dico.clear()

    melted = pd.melt(barcode_selection_phred_scrore_dataframe, id_vars=["passes_filtering"],
                     var_name="barcodes", value_name="qscore")

    ax = sns.boxplot(data=melted, x='barcodes', y='qscore', hue='passes_filtering', showfliers=False,
                palette={True: "yellowgreen", False: "orangered"}, hue_order=[True, False])
    handles, _ = ax.get_legend_handles_labels()
    plt.legend(bbox_to_anchor=(0.92, 1), loc=2, borderaxespad=0., labels=["Pass", "Fail"], handles=handles)
    plt.xlabel('Barcodes')
    plt.ylabel('Mean Phred score')
    plt.title(main)

    table_html=_make_table_html(barcode_selection_phred_scrore_dataframe)

    plt.savefig(output_file)
    plt.close()

    return main, output_file, table_html, desc


#  1Dsqr plots

def dsqr_read_count_histogram(albacore_log_1d, albacore_log_1dsqr, main, my_dpi, result_directory, desc):
    """
    Plots the count histograms of count  of the different types of reads eventually available in a Minion run: 1D, full_1Dsquare.
    """
    output_file = result_directory + '/' + main + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    fast5_template_basecalled = len(albacore_log_1d[albacore_log_1d['num_called_template'] != 0])
    fast5_1dsqr = len(albacore_log_1dsqr['passes_filtering'])
    read_pass_1dsqr = len(albacore_log_1dsqr[albacore_log_1dsqr['passes_filtering'] == True])
    read_fail_1dsqr = len(albacore_log_1dsqr[albacore_log_1dsqr['passes_filtering'] == False])

    read_type = [fast5_template_basecalled, fast5_1dsqr, read_pass_1dsqr, read_fail_1dsqr]
    label = ("1D", "1Dsquare", "1Dsquare pass", "1Dsquare fail")
    nd = np.arange(len(read_type))

    bars = plt.bar(nd, read_type, align='center', color=["salmon", "orange", "yellowgreen", "orangered"])
    plt.xticks(nd, label)
    plt.xlabel("Read type")
    plt.ylabel("Counts")
    plt.title(main)

    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2., 1 * height, '%d' % int(height), ha='center', va='bottom')
    table_html = None

    plt.savefig(output_file)
    plt.close()

    return main, output_file, table_html, desc

def dsqr_read_length_multihistogram(albacore_log_1d, albacore_log_1dsqr, main, my_dpi, result_directory, desc):
    """
    Plots an histogram of the reads length by bins of 100 for each of the barcodes described in the design file or without barcode
    """
    output_file = result_directory + '/' + main + '.png'
    read_1d = albacore_log_1d.sequence_length_template.loc[albacore_log_1d['num_called_template'] != 0]
    read_1dsqr = albacore_log_1dsqr.loc[:,"sequence_length_2d"]
    read_pass_1dsqr = albacore_log_1dsqr.sequence_length_2d.loc[True == albacore_log_1dsqr['passes_filtering']]
    read_fail_1dsqr = albacore_log_1dsqr.sequence_length_2d.loc[False == albacore_log_1dsqr['passes_filtering']]
    minimum, maximum = min(read_1d), max(read_1d)
    read_type=['1D','1Dsquare','1Dsquare pass','1Dsquare fail']

    plt.figure(figsize=(12, 7), dpi=my_dpi)
    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])
    plt.subplots_adjust(bottom=0.015, top=1)

    n, bins, patches = ax.hist([read_1d,read_1dsqr,read_pass_1dsqr,read_fail_1dsqr],color=["salmon","goldenrod","yellowgreen", "orangered"], edgecolor='black', label=read_type,
                               bins=2 ** np.linspace(_safe_log(minimum), _safe_log(maximum), 30))
    plt.legend()

    ax.set_xscale('log', basex=2)
    ax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.1f}'))
    ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))
    ax.set_xticks(bins)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_rotation('vertical')
    ax.set_xlabel('Read length(bp)')
    ax.set_ylabel('Read number')
    ax.set_title(main)

    dataframe = pd.DataFrame({"1D": read_1d,'1Dsquare': read_1dsqr ,"1Dsquare pass": read_pass_1dsqr, "1Dsquare fail": read_fail_1dsqr})

    dataframe = dataframe[["1D","1Dsquare","1Dsquare pass","1Dsquare fail"]]
    table_html = _make_table_html(dataframe)

    # ax2.xaxis.set_visible(False)
    # ax2.yaxis.set_visible(False)
    # ax2.axis('off')

    #the_table.set_fontsize(12)
    #the_table.scale(1, 1)

    plt.savefig(output_file)
    plt.close()

    return main, output_file, table_html, desc

def dsqr_read_quality_multiboxplot(albacore_log_1d, albacore_log_1dsqr, main, my_dpi, result_directory, desc):
    """
    Plots a boxplot of reads quality
    """
    output_file = result_directory + '/' + main + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    plt.subplots_adjust(bottom=0.015, top=1.0)
    mean_qscore_1d = albacore_log_1d.loc[:,"mean_qscore_template"]
    mean_qscore_1dsqr = albacore_log_1dsqr.loc[:,"mean_qscore_2d"]
    read_pass_1dsqr = albacore_log_1dsqr.mean_qscore_2d.loc[True == albacore_log_1dsqr['passes_filtering']]
    read_fail_1dsqr = albacore_log_1dsqr.mean_qscore_2d.loc[False == albacore_log_1dsqr['passes_filtering']]

    my_pal = {"1D": "salmon","1Dsquare": "goldenrod", "1Dsquare pass": "yellowgreen", "1Dsquare fail": "orangered"}
    order=["1D","1Dsquare","1Dsquare pass","1Dsquare fail"]
    dataframe = pd.DataFrame({"1D": mean_qscore_1d,"1Dsquare": mean_qscore_1dsqr, "1Dsquare pass": read_pass_1dsqr, "1Dsquare fail": read_fail_1dsqr})
    sns.boxplot(data=dataframe, ax=plt.subplot(gs[0]),palette=my_pal,order=order)


    plt.ylabel('Mean Phred score')
    plt.title(main)
    plt.legend()

    dataframe = dataframe[["1D","1Dsquare","1Dsquare pass","1Dsquare fail"]]
    table_html = _make_table_html(dataframe)

    plt.savefig(output_file)
    plt.close()

    return main, output_file, table_html, desc

def dsqr_phred_score_frequency(albacore_log_1dsqr, main, my_dpi, result_directory, desc):
    '''
    Plot the distribution of the phred score
    '''
    output_file = result_directory + '/' + main + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    plt.subplots_adjust(bottom=0.015, top=1.0)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])
    qscore_1dsqr = albacore_log_1dsqr['mean_qscore_2d']

    sns.distplot(qscore_1dsqr, bins=15, color='goldenrod',
                 hist_kws=dict(edgecolor="k", linewidth=1),hist=True,label="1Dsquare")

    plt.legend()
    plt.xlabel("Mean Phred score")
    plt.ylabel("Frequency")
    plt.title(main)

    dataframe = pd.DataFrame({"1Dsquare": qscore_1dsqr})
    rd = dataframe.describe().drop('count').round(2).reset_index()

    plt.axvline(x=qscore_1dsqr.describe()['50%'], color='goldenrod')

    plt.savefig(output_file)
    plt.close()

    table_html = _make_table_html(rd)

    return main, output_file, table_html, desc

def dsqr_allphred_score_frequency(albacore_log_1d, albacore_log_1dsqr, main, my_dpi, result_directory, desc):
    '''
    Plot the distribution of the phred score
    '''
    output_file = result_directory + '/' + main + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    plt.subplots_adjust(bottom=0.015, top=1.0)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])
    qscore_1d = albacore_log_1d['mean_qscore_template']
    qscore_1dsqr = albacore_log_1dsqr['mean_qscore_2d']
    mean_qscore_read_pass = albacore_log_1dsqr.mean_qscore_2d.loc[True == albacore_log_1dsqr['passes_filtering']]
    mean_qscore_read_fail = albacore_log_1dsqr.mean_qscore_2d.loc[False == albacore_log_1dsqr['passes_filtering']]

    sns.distplot(mean_qscore_read_pass, bins=15,hist_kws=dict(edgecolor="k", linewidth=1) ,color='yellowgreen',
                 hist=True, label='1Dsquare pass')
    sns.distplot(mean_qscore_read_fail, bins=15,hist_kws=dict(edgecolor="k", linewidth=1) ,color='orangered', label='1Dsquare fail',
                 hist=True)
    plt.legend()
    plt.xlabel("Mean Phred score")
    plt.ylabel("Frequency")
    plt.title(main)

    dataframe = pd.DataFrame({"1D": qscore_1d,"1Dsquare": qscore_1dsqr, "1Dsquare pass": mean_qscore_read_pass, "1Dsquare fail": mean_qscore_read_fail})

    rd = dataframe.describe().drop('count').round(2).reset_index()

    plt.axvline(x=mean_qscore_read_pass.describe()['50%'], color='yellowgreen')
    plt.axvline(x=mean_qscore_read_fail.describe()['50%'], color='orangered')

    rd = rd[["1D","1Dsquare","1Dsquare pass","1Dsquare fail"]]
    table_html = _make_table_html(rd)

    plt.savefig(output_file)
    plt.close()

    return main, output_file, table_html, desc

def scatterplot_1dsqr(albacore_log_1d,albacore_log_1dsqr, main, my_dpi, result_directory, desc):
    '''
    Plot the scatter plot representing the relation between the phred score and the sequence length
    '''
    output_file = result_directory + '/' + main + '.png'
    #plt.figure(figsize=(1100/my_dpi, 6), dpi=my_dpi)
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    #plt.subplots_adjust(bottom=0.015, top=1.0)

    length_1dsqr_read_pass = albacore_log_1dsqr.sequence_length_2d.loc[True == albacore_log_1dsqr['passes_filtering']]
    length_1dsqr_read_fail = albacore_log_1dsqr.sequence_length_2d.loc[False == albacore_log_1dsqr['passes_filtering']]

    qscore_1dsqr_read_pass = albacore_log_1dsqr.mean_qscore_2d.loc[True == albacore_log_1dsqr['passes_filtering']]
    qscore_1dsqr_read_fail = albacore_log_1dsqr.mean_qscore_2d.loc[False == albacore_log_1dsqr['passes_filtering']]

    read_pass=plt.scatter(x=length_1dsqr_read_pass, y=qscore_1dsqr_read_pass,color="yellowgreen")
    read_fail=plt.scatter(x=length_1dsqr_read_fail, y=qscore_1dsqr_read_fail,color="orangered")
    plt.legend((read_pass,read_fail),("1Dsquare pass","1Dsquare fail"))
    plt.xlim(np.min(albacore_log_1dsqr.sequence_length_2d.loc[albacore_log_1dsqr.sequence_length_2d>0]),np.max(albacore_log_1dsqr.sequence_length_2d))
    plt.xscale('log')
    plt.xlabel("Sequence length")
    plt.ylabel("Mean Phred score")
    plt.title(main)
    plt.savefig(output_file)
    plt.close()

    table_html = None
    return main, output_file, table_html, desc

# For eache barcode 1Dsqr

def barcode_percentage_pie_chart_1dsqr_pass(albacore_log, main, barcode_selection, my_dpi, result_directory, desc):
    """
    Plots a pie chart of the barcode percentage of a run. Needs the design file describing the barcodes to run
    """
    output_file = result_directory + '/' + main + '.png'
    plt.figure(figsize=(800 / my_dpi, 800 / my_dpi), dpi=my_dpi)
    for element in barcode_selection:

        if all(albacore_log['barcode_arrangement'] != element):
            print("The barcode {} doesn't exist".format(element))
            return False

    barcode = albacore_log.barcode_arrangement.loc[True == albacore_log['passes_filtering']]
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
        ax1.set_title(main)
        ax1.bar(length, barcode_count, color=cs)
        ax1.set_xticks(length)
        ax1.set_xticklabels(barcode_selection)
    table_html = None

    plt.savefig(output_file)
    plt.close()
    return main, output_file, table_html, desc

def barcode_percentage_pie_chart_1dsqr_fail(albacore_log, main, barcode_selection, my_dpi, result_directory, desc):
    """
    Plots a pie chart of the barcode percentage of a run. Needs the design file describing the barcodes to run
    """
    output_file = result_directory + '/' + main + '.png'
    plt.figure(figsize=(800 / my_dpi, 800 / my_dpi), dpi=my_dpi)
    for element in barcode_selection:

        if all(albacore_log['barcode_arrangement'] != element):
            print("The barcode {} doesn't exist".format(element))
            return False

    barcode = albacore_log.barcode_arrangement.loc[False == albacore_log['passes_filtering']]
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
        ax1.set_title(main)
        ax1.bar(length, barcode_count, color=cs)
        ax1.set_xticks(length)
        ax1.set_xticklabels(barcode_selection)
    table_html = None

    plt.savefig(output_file)
    plt.close()

    return main, output_file, table_html, desc

def barcode_length_boxplot_1dsqr(albacore_log, main, barcode_selection, my_dpi, result_directory, desc):
    '''
    Plot the length boxplot for each barcode indicated in the sample sheet
    '''
    output_file = result_directory + '/' + main + '.png'
    pattern = '(\d{2})'
    dico = {}

    fig = plt.figure(figsize=(12, 7), dpi=my_dpi)
    plt.subplots_adjust(bottom=0.015, top=1.0)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])


    dico["passes_filtering"] = albacore_log['passes_filtering']

    for barcode in barcode_selection:
        barcode_selected_dataframe = albacore_log[albacore_log['barcode_arrangement'] == barcode]

        match = re.search(pattern, barcode)
        if match:
            dico[match.group(0)] = barcode_selected_dataframe['sequence_length_2d']
        else:
            dico['Unclassified'] = barcode_selected_dataframe['sequence_length_2d']


    barcode_selection_sequence_length_dataframe = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dico.items()]))
    dico.clear()
    melted = pd.melt(barcode_selection_sequence_length_dataframe, id_vars=["passes_filtering"],
                     var_name="barcodes", value_name="length")
    ax = sns.boxplot(data=melted,x='barcodes',y='length',hue='passes_filtering', showfliers=False,palette={True : "yellowgreen", False : "orangered"},hue_order=[True,False])
    handles, _ = ax.get_legend_handles_labels()
    plt.legend(bbox_to_anchor=(0.92, 1), loc=2, borderaxespad=0.,labels=["Pass","Fail"],handles=handles)
    plt.xlabel('Barcodes')
    plt.ylabel('Read length(bp)')
    plt.title(main)

    table_html = _make_table_html(barcode_selection_sequence_length_dataframe)
    plt.savefig(output_file)
    plt.close()
    return main, output_file, table_html, desc

def barcoded_phred_score_frequency_1dsqr(albacore_log, main, barcode_selection, my_dpi, result_directory, desc):
    '''
    Plot the 1Dsquare phred score distribution boxplot for each barcode indicated in the sample sheet
    '''
    output_file = result_directory + '/' + main + '.png'
    dico = {}
    pattern = '(\d{2})'

    fig = plt.figure(figsize=(12, 7), dpi=my_dpi)
    plt.subplots_adjust(bottom=0.015, top=1.0)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])

    dico["passes_filtering"] = albacore_log['passes_filtering']

    for barcode in barcode_selection:
        barcode_selected_dataframe = albacore_log[albacore_log['barcode_arrangement'] == barcode]
        match = re.search(pattern, barcode)
        if match:
            dico[match.group(0)] = barcode_selected_dataframe['mean_qscore_2d']
        else:
            dico['Unclassified'] = barcode_selected_dataframe['mean_qscore_2d']
    barcode_selection_phred_scrore_dataframe = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dico.items()]))
    dico.clear()
    melted = pd.melt(barcode_selection_phred_scrore_dataframe, id_vars=["passes_filtering"],
                     var_name="barcodes", value_name="qscore")
    ax = sns.boxplot(data=melted,x='barcodes',y='qscore',hue='passes_filtering', showfliers=False,palette={True : "yellowgreen", False : "orangered"},hue_order=[True,False])
    handles, _ = ax.get_legend_handles_labels()
    plt.legend(bbox_to_anchor=(0.92, 1), loc=2, borderaxespad=0., labels=["Pass", "Fail"], handles=handles)
    plt.xlabel('Barcodes')
    plt.ylabel('Mean Phred score')
    plt.title(main)

    table_html = _make_table_html(barcode_selection_phred_scrore_dataframe)
    plt.savefig(output_file)
    plt.close()

    return main, output_file, table_html, desc