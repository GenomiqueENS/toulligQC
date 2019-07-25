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
#
# First author: Lionel Ferrato-Berberian
# Maintainer: Bérengère Laffay
# Since version 0.1

# Functions to generate graphs and statistics tables in HTML format, they use the result_dict dictionary.

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib.pyplot import table


def _is_in_result_dict(dict, dict_key, default_value):
    """
    Global function to check for the presence of an entry in a dictionary
    and give it a default value.
    :param result_dict: result_dict dictionary
    :param dict_key: entry (string)
    :param default_value:
    :return:
    """
    if dict_key not in dict or not dict[dict_key]:
        dict[dict_key] = default_value
    return dict[dict_key]


def _make_desribe_dataframe(value):
    """
    Creation of a statistics table printed with the graph in report.html
    :param value: information measured (series)
    """

    desc = value.describe()
    desc.loc['count'] = desc.loc['count'].astype(int).astype(str)
    desc.iloc[1:] = desc.iloc[1:].applymap(lambda x: '%.2f' % x)
    desc.rename({'50%': 'median'}, axis='index', inplace=True)
    # desc.iloc[1:] = desc.iloc[1:].applymap('{:.2f}'.format) (pandas > 0.23)

    return desc


def _safe_log(x):
    """
    Verification that we haven't a null value
    :param x: tested value
    :return: log2 value or 0
    """
    if x <= 0:
        return 0
    return np.log2(x)

#
#  1D plots
#


def read_count_histogram(result_dict, main, my_dpi, result_directory, desc):
    """
    Plots the histogram of count of the different types of reads:
    Fast5 submitted to MinKNOW
    FastQ return by MinKNOW
    1D read return by Albacore
    1D pass read return by Albacore (Qscore > 7.5)
    1D fail read return by Albacore (Qscore < 7.5)
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'

    # Histogram completed with the number of basecalling errors found in the pipeline.log file

    if 'albacore.log.extractor.fast5.files.submitted' in result_dict or result_dict.keys().__contains__('albacore.log.extractor.fast5.files.submitted'):

        # Histogram completed with the number of basecalling reads (is.barcode == True)
        if 'basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count' in result_dict or result_dict.keys().__contains__('basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count'):

            plt.figure(figsize=(16, 7), dpi=my_dpi)
            gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
            ax = plt.subplot(gs[0])
            plt.mec = 'black'
            plt.mfc = 'white'
            plt.subplots_adjust(bottom=0.015, top=1.0)

            read_type = [result_dict['albacore.log.extractor.fast5.files.submitted'],
                         result_dict['albacore.log.extractor.fast5.files.basecalled.error.count'],
                         result_dict['basecaller.sequencing.summary.1d.extractor.fastq.entries'],
                         result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
                         result_dict['basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.count'],
                         result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
                         result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count"],
                         result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"],
                         result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.count"]
                         ]
            label = ("Raw Fast5", "Raw Fast5 with error", "FastQ entries", "1D",
                     'Null sequence length', "1D pass","1D pass barcoded" ,"1D fail","1D fail barcoded")
            nd = np.arange(len(read_type))
            bars = ax.bar(nd, read_type, align='center', color=["Green", "yellow", "lightblue", "salmon",
                                                                'purple', "yellowgreen", "lightgoldenrodyellow", "orangered","darksalmon"],
                          edgecolor="black", linewidth=1)

            array = np.array([[result_dict['albacore.log.extractor.fast5.files.submitted'],
                               result_dict['albacore.log.extractor.fast5.files.basecalled.error.count'],
                               result_dict["basecaller.sequencing.summary.1d.extractor.fastq.entries"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
                               result_dict['basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.count'],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.count"]],
                              [result_dict['albacore.log.extractor.fast5.files.frequency'],
                               result_dict['albacore.log.extractor.fast5.files.basecalled.error.frequency'],
                               result_dict['basecaller.sequencing.summary.1d.extractor.fastq.entries.frequency'],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.count.frequency"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.frequency"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.frequency"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.frequency"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.frequency"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.frequency"]]])
            dataframe = pd.DataFrame(array, index=['count', 'frequency'],
                                     columns=["Raw Fast5", "Raw Fast5 with error", "FastQ_entries", "1D",
                                              'Null sequence length', "1D pass","1D pass barcoded" ,"1D fail","1D fail barcoded"])

        else:

            plt.figure(figsize=(12, 7), dpi=my_dpi)
            gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
            ax = plt.subplot(gs[0])
            plt.mec = 'black'
            plt.mfc = 'white'
            plt.subplots_adjust(bottom=0.015, top=1.0)

            read_type = [result_dict['albacore.log.extractor.fast5.files.submitted'],
                             result_dict['albacore.log.extractor.fast5.files.basecalled.error.count'],
                             result_dict['basecaller.sequencing.summary.1d.extractor.fastq.entries'],
                             result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
                             result_dict['basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.count'],
                             result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
                             result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"]
                             ]
            label = ("Raw Fast5", "Raw Fast5 with error", "FastQ entries", "1D",
                     'Null sequence length', "1D pass","1D fail")
            nd = np.arange(len(read_type))
            bars = ax.bar(nd, read_type, align='center', color=["Green", "yellow", "lightblue", "salmon",
                                                                'purple', "yellowgreen", "orangered"],
                          edgecolor="black", linewidth=1)

            array = np.array([[result_dict['albacore.log.extractor.fast5.files.submitted'],
                               result_dict['albacore.log.extractor.fast5.files.basecalled.error.count'],
                               result_dict["basecaller.sequencing.summary.1d.extractor.fastq.entries"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
                               result_dict['basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.count'],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"]],
                              [result_dict['albacore.log.extractor.fast5.files.frequency'],
                               result_dict['albacore.log.extractor.fast5.files.basecalled.error.frequency'],
                               result_dict['basecaller.sequencing.summary.1d.extractor.fastq.entries.frequency'],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.count.frequency"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.frequency"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.frequency"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.frequency"]]])
            dataframe = pd.DataFrame(array, index=['count', 'frequency'],
                                     columns=["Raw Fast5", "Raw Fast5 with error", "FastQ_entries", "1D",
                                              'Null sequence length', "1D pass","1D fail"])

    else:
        plt.figure(figsize=(12, 7), dpi=my_dpi)
        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
        ax = plt.subplot(gs[0])
        plt.mec = 'black'
        plt.mfc = 'white'
        plt.subplots_adjust(bottom=0.015, top=1.0)

        # Histogram completed with the number of basecalling reads (is.barcode == True)
        if 'basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count' in result_dict or result_dict.keys().__contains__('basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count'):

            read_type = [result_dict['basecaller.sequencing.summary.1d.extractor.fastq.entries'],
                         result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
                         result_dict['basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.count'],
                         result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
                         result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count"],
                         result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"],
                         result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.count"]
                         ]
            label = ("FastQ entries", "1D",'Null sequence length', "1D pass","1D pass barcoded" ,"1D fail","1D fail barcoded")
            nd = np.arange(len(read_type))
            bars = ax.bar(nd, read_type, align='center', color=["lightblue", "salmon",'purple', "yellowgreen", "lightgoldenrodyellow", "orangered","darksalmon"],
                          edgecolor="black", linewidth=1)

            array = np.array([[result_dict["basecaller.sequencing.summary.1d.extractor.fastq.entries"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
                               result_dict['basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.count'],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.count"]],
                              [result_dict['basecaller.sequencing.summary.1d.extractor.fastq.entries.frequency'],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.count.frequency"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.frequency"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.frequency"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.frequency"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.frequency"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.frequency"]]])
            dataframe = pd.DataFrame(array, index=['count', 'frequency'],
                                     columns=["FastQ_entries", "1D",'Null sequence length', "1D pass","1D pass barcoded" ,"1D fail","1D fail barcoded"])
        else:

            read_type = [result_dict['basecaller.sequencing.summary.1d.extractor.fastq.entries'],
                         result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
                         result_dict['basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.count'],
                         result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
                         result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"]
                         ]
            label = ("FastQ entries", "1D",'Null sequence length', "1D pass","1D pass barcoded" ,"1D fail","1D fail barcoded")
            nd = np.arange(len(read_type))
            bars = ax.bar(nd, read_type, align='center', color=["lightblue", "salmon",'purple', "yellowgreen","orangered"],
                          edgecolor="black", linewidth=1)

            array = np.array([[result_dict["basecaller.sequencing.summary.1d.extractor.fastq.entries"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
                               result_dict['basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.count'],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"]],
                              [result_dict['basecaller.sequencing.summary.1d.extractor.fastq.entries.frequency'],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.count.frequency"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.frequency"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.frequency"],
                               result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.frequency"]]])
            dataframe = pd.DataFrame(array, index=['count', 'frequency'],
                                     columns=["FastQ_entries", "1D", 'Null sequence length', "1D pass", "1D fail"])

    plt.xticks(nd, label)
    plt.xlabel("Read type")
    plt.ylabel("Counts")
    plt.title(main)

    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2., 1 * height, '%d' % int(height), ha='center', va='bottom')

    plt.savefig(output_file)
    plt.close()

    dataframe.iloc[0] = dataframe.iloc[0].astype(int).astype(str)
    dataframe.iloc[1:] = dataframe.iloc[1:].applymap('{:.2f}'.format)
    table_html = pd.DataFrame.to_html(dataframe)

    return main, output_file, table_html, desc


def read_length_multihistogram(result_dict, main, my_dpi, result_directory, desc):
    """
    Plots an histogram of the read length for the different types of read:
    1D, 1Dpass, 1D fail
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'

    minimum, maximum = \
        min(result_dict["basecaller.sequencing.summary.1d.extractor.sequence.length"]), \
        max(result_dict["basecaller.sequencing.summary.1d.extractor.sequence.length"])
    read_type = ['1D', '1D pass', '1D fail']

    plt.figure(figsize=(12, 7), dpi=my_dpi)
    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])
    plt.subplots_adjust(bottom=0.015, top=1.0)

    data = [result_dict["basecaller.sequencing.summary.1d.extractor.sequence.length"],
            result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.length"],
            result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.length"]]
    ls = 2 ** np.linspace(_safe_log(minimum), _safe_log(maximum), 30)
    n, bins, patches = ax.hist(data, color=["salmon", "yellowgreen", "orangered"],
                               edgecolor='black', label=read_type,
                               bins=ls)
    plt.legend()

    ax.set_xscale('log', basex=2)
    ax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.0f}'))
    ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.0f}'))
    ax.set_xticks(bins)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_rotation('vertical')
    ax.set_xlabel('Read length(bp)')
    ax.set_ylabel('Read number')
    ax.set_title(main)
    dataframe = \
        pd.DataFrame({"1D": result_dict["basecaller.sequencing.summary.1d.extractor.sequence.length"],
                      "1D pass": result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.length"],
                      "1D fail": result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.length"]})
    dataframe = dataframe[["1D", "1D pass", "1D fail"]]

    plt.savefig(output_file)
    plt.close()

    table_html = pd.DataFrame.to_html(_make_desribe_dataframe(dataframe))

    return main, output_file, table_html, desc


def allread_number_run(result_dict, main, my_dpi, result_directory, desc):
    """
    Plots the different reads (1D, 1D pass, 1D fail) produced along the run against the time(in hour)
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi)

    plt.plot(result_dict["basecaller.sequencing.summary.1d.extractor.start.time.sorted"],
             np.arange(len(result_dict["basecaller.sequencing.summary.1d.extractor.start.time.sorted"])),
             color='salmon', linewidth=1, label="1D")

    plt.plot(result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.sorted"],
             np.arange(len(result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.sorted"])),
             color='yellowgreen', linewidth=1, label="1D pass")

    plt.plot(result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.sorted"],
             np.arange(len(result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.sorted"])),
             color='orangered', linewidth=1, label="1D fail")

    plt.xticks(np.arange(0, max(result_dict["basecaller.sequencing.summary.1d.extractor.start.time.sorted"]), 8))

    plt.ylabel("Read number")
    plt.xlabel("Time (Hour)")
    plt.legend()

    plt.savefig(output_file)
    plt.close()

    table_html = None

    return main, output_file, table_html, desc


def read_quality_multiboxplot(result_dict, main, my_dpi, result_directory, desc):
    """
    Plots a boxplot of reads quality per read type (1D, 1D pass, 1D fail)
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi, facecolor='white', edgecolor='black')
    gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[2, 1])

    my_pal = {"1D": "salmon", "1D pass": "yellowgreen", "1D fail": "orangered"}
    order = ["1D", "1D pass", "1D fail"]
    dataframe = \
        pd.DataFrame({"1D": result_dict["basecaller.sequencing.summary.1d.extractor.mean.qscore"],
                      "1D pass": result_dict['basecaller.sequencing.summary.1d.extractor.read.pass.qscore'],
                      "1D fail": result_dict['basecaller.sequencing.summary.1d.extractor.read.fail.qscore']})
    ax = plt.subplot(gs[0])
    sns.boxplot(data=dataframe, ax=ax, palette=my_pal, order=order, linewidth=1)
    plt.ylabel('Mean Phred score')
    ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.0f}'))

    ax2 = plt.subplot(gs[1])
    sns.violinplot(data=dataframe, ax=ax2, palette=my_pal, inner=None, cut=0, order=order, linewidth=1)
    ax2.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.0f}'))
    plt.subplots_adjust(bottom=0.015, top=1.0)
    plt.ylabel('Mean Phred score')
    plt.title(main)

    dataframe = dataframe[["1D", "1D pass", "1D fail"]]
    table_html = pd.DataFrame.to_html(_make_desribe_dataframe(dataframe))

    plt.savefig(output_file)
    plt.close()
    return main, output_file, table_html, desc


def phred_score_frequency(result_dict, main, my_dpi, result_directory, desc):
    """
    Plot the distribution of the phred score (not use anymore)
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    plt.subplots_adjust(bottom=0.015, top=1.0)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    plt.subplot(gs[0])
    plt.subplots_adjust(bottom=0.015, top=1.0)

    sns.distplot(result_dict["basecaller.sequencing.summary.1d.extractor.mean.qscore"], bins=15, color='salmon',
                 hist_kws=dict(edgecolor="k", linewidth=1), hist=True, label="1D")

    plt.legend()
    plt.xlabel("Mean Phred score")
    plt.ylabel("Frequency")
    plt.title(main)

    dataframe = pd.DataFrame({"1D": result_dict["basecaller.sequencing.summary.1d.extractor.mean.qscore"]})
    rd = dataframe.describe().drop('count').round(2).reset_index()

    plt.axvline(x=result_dict["basecaller.sequencing.summary.1d.extractor.mean.qscore"].describe()['50%'], color='salmon')

    plt.savefig(output_file)
    plt.close()

    table_html = pd.DataFrame.to_html(rd)

    return main, output_file, table_html, desc


def allphred_score_frequency(result_dict, main, my_dpi, result_directory, desc):
    """
    Plot the distribution of the phred score per read type (1D , 1D pass, 1D fail)
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])
    plt.subplots_adjust(bottom=0.015, top=1.0)

    sns.distplot(result_dict["basecaller.sequencing.summary.1d.extractor.mean.qscore"], bins=15, ax=ax, color='salmon',
                 hist_kws=dict(edgecolor="k", linewidth=1), hist=True, label="1D")
    plt.axvline(x=result_dict["basecaller.sequencing.summary.1d.extractor.mean.qscore"].describe()['50%'], color='salmon')

    ax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.0f}'))
    ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.2f}'))

    plt.legend()
    plt.xlabel("Mean Phred score")
    plt.ylabel("Frequency")
    plt.title(main)

    ax2 = plt.subplot(gs[1])

    max_frequency_pass = max(result_dict['basecaller.sequencing.summary.1d.extractor.read.pass.qscore'].value_counts(normalize=True))

    max_frequency_fail = max(result_dict['basecaller.sequencing.summary.1d.extractor.read.fail.qscore'].value_counts(normalize=True))

    if max_frequency_pass > max_frequency_fail:
        sns.distplot(result_dict['basecaller.sequencing.summary.1d.extractor.read.pass.qscore'],
                     ax=ax2, bins=15, hist_kws=dict(edgecolor="k", linewidth=1),
                     color='yellowgreen', hist=True, label='1D pass')

        sns.distplot(result_dict['basecaller.sequencing.summary.1d.extractor.read.fail.qscore'],
                     ax=ax2, bins=15, hist_kws=dict(edgecolor="k", linewidth=1),
                     color='orangered', label='1D fail', hist=True)

    else:
        sns.distplot(result_dict['basecaller.sequencing.summary.1d.extractor.read.fail.qscore'],
                     ax=ax2, bins=15, hist_kws=dict(edgecolor="k", linewidth=1),
                     color='orangered', label='1D fail', hist=True)

        sns.distplot(result_dict['basecaller.sequencing.summary.1d.extractor.read.pass.qscore'],
                     ax=ax2, bins=15, hist_kws=dict(edgecolor="k", linewidth=1),
                     color='yellowgreen', hist=True, label='1D pass')

    ax2.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.0f}'))
    ax2.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.2f}'))

    plt.legend()
    plt.xlabel("Mean Phred score")
    plt.ylabel("Frequency")
    plt.title(main)
    plt.axvline(x=result_dict['basecaller.sequencing.summary.1d.extractor.read.pass.qscore'].describe()['50%'], color='yellowgreen')
    plt.axvline(x=result_dict['basecaller.sequencing.summary.1d.extractor.read.fail.qscore'].describe()['50%'], color='orangered')

    dataframe = \
        pd.DataFrame({"1D": result_dict["basecaller.sequencing.summary.1d.extractor.mean.qscore"],
                      "1D pass": result_dict['basecaller.sequencing.summary.1d.extractor.read.pass.qscore'],
                      "1D fail": result_dict['basecaller.sequencing.summary.1d.extractor.read.fail.qscore']})

    plt.savefig(output_file)
    plt.close()

    dataframe = _make_desribe_dataframe(dataframe).drop('count')

    table_html = pd.DataFrame.to_html(dataframe)

    return main, output_file, table_html, desc


def all_scatterplot(result_dict, main, my_dpi, result_directory, desc):
    """
    Plot the scatter plot representing the relation between the phred score and the sequence length in log
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'
    plt.figure(figsize=(12, 6), dpi=my_dpi)
    ax = plt.gca()

    read_pass = plt.scatter(x=result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.length"],
                            y=result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.qscore"],
                            color="yellowgreen")

    read_fail = plt.scatter(x=result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.length"],
                            y=result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.qscore"],
                            color="orangered")

    plt.legend((read_pass, read_fail), ("1D pass", "1D fail"))
    plt.xlim(np.min(result_dict["basecaller.sequencing.summary.1d.extractor.sequence.length"]
                    .loc[result_dict["basecaller.sequencing.summary.1d.extractor.sequence.length"] > 0]),
             np.max(result_dict["basecaller.sequencing.summary.1d.extractor.sequence.length"]))

    plt.yticks()
    plt.xscale('log')
    plt.xlabel("Sequence length")
    plt.ylabel("Mean Phred score")
    ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.0f}'))
    plt.title(main)
    plt.savefig(output_file)
    plt.close()

    table_html = None

    return main, output_file, table_html, desc


def channel_count_histogram(albacore_log, main, my_dpi, result_directory, desc):
    """
    Plots an histogram of the channel count according to the channel number (not use anymore)
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'
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
    plt.subplot(gs[1])

    dataframe = table(ax, np.round(total_number_reads_per_channel
                                   .describe().drop(['mean', 'std', '50%', '75%', '25%']), 2), loc='center')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.axis('off')

    dataframe.set_fontsize(12)
    dataframe.scale(1, 1.2)

    plt.savefig(output_file)
    plt.close()
    table_html = pd.DataFrame.to_html(total_number_reads_per_channel.describe())

    return main, output_file, table_html, desc


def _minion_flowcell_layout():
    """
    Represents the layout of a minion flowcell (not use anymore)
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
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'
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
    df.pivot("Row number", "Column number", "labels")
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    sns.heatmap(d, fmt="", linewidths=.5, cmap="YlGnBu", annot_kws={"size": 7},
                cbar_kws={'label': 'Read number per pore channel', "orientation": "horizontal"})
    plt.title(main)

    plt.savefig(output_file)
    plt.close()

    table_html = None

    return main, output_file, table_html, desc

#
# For each barcode 1D
#


def barcode_percentage_pie_chart_pass(result_dict, main, barcode_selection, my_dpi, result_directory, desc):
    """
    Plots a pie chart of 1D read pass percentage per barcode of a run.
    Needs the samplesheet file describing the barcodes to run
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'
    plt.figure(figsize=(800 / my_dpi, 800 / my_dpi), dpi=my_dpi)
    for element in barcode_selection:

        if all(result_dict['basecaller.sequencing.summary.1d.extractor.barcode.arrangement'] != element):
            print("The barcode {} doesn't exist".format(element))

    count_sorted = result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded"]
    barcodes = count_sorted.index.values.tolist()

    cs = plt.get_cmap('Spectral')(np.arange(len(barcodes)) / len(barcodes))

    sizes = [(100 * chiffre) / sum(count_sorted) for chiffre in count_sorted.values]
    if len(barcode_selection) <= 10:
        fig1, ax1 = plt.subplots()
        ax1.pie(sizes, labels=None, startangle=90, colors=cs, wedgeprops={'linewidth': 1, 'edgecolor': 'k'})
        ax1.axis('equal')
        ax1.legend(labels=['%s, %1.1f %%' % (l, s) for l, s in zip(barcodes, sizes)],
                   loc="upper right", bbox_to_anchor=(1.1, 1.175), edgecolor="black")

    else:
        fig = plt.figure(figsize=(20, 10))
        ax1 = fig.add_subplot(111)
        length = np.arange(0, len(count_sorted))
        ax1.set_title(main)
        ax1.bar(length, count_sorted, color=cs)
        ax1.set_xticks(length)
        ax1.set_xticklabels(barcodes)
        plt.legend(labels=['%s, %1.1f %%' % (l, s) for l, s in zip(barcodes, sizes)],
                   loc="upper right", bbox_to_anchor=(1.1, 1.175))

    plt.savefig(output_file)
    plt.close()

    barcode_table = pd.DataFrame({"barcode arrangement": count_sorted/sum(count_sorted)*100,
                                 "read count": count_sorted})
    barcode_table.sort_index(inplace=True)
    pd.options.display.float_format = '{:.2f}%'.format
    table_html = pd.DataFrame.to_html(barcode_table)

    return main, output_file, table_html, desc


def barcode_percentage_pie_chart_fail(result_dict, main, barcode_selection, my_dpi, result_directory, desc):
    """
    Plots a pie chart of 1D read fail percentage per barcode of a run.
    Needs the samplesheet file describing the barcodes to run
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'
    plt.figure(figsize=(800 / my_dpi, 800 / my_dpi), dpi=my_dpi)
    for element in barcode_selection:

        if all(result_dict['basecaller.sequencing.summary.1d.extractor.barcode.arrangement'] != element):
            print("The barcode {} doesn't exist".format(element))
            return False

    count_sorted = result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded"]
    barcodes = count_sorted.index.values.tolist()

    cs = plt.get_cmap('Spectral')(np.arange(len(barcodes)) / len(barcodes))

    sizes = [(100 * chiffre) / sum(count_sorted) for chiffre in count_sorted.values]
    if len(barcode_selection) <= 10:
        fig1, ax1 = plt.subplots()
        ax1.pie(sizes, labels=None, startangle=90, colors=cs, wedgeprops={'linewidth': 1, 'edgecolor': 'k'})
        ax1.axis('equal')
        ax1.legend(labels=['%s, %1.1f %%' % (l, s) for l, s in zip(barcodes, sizes)],
                   loc="upper right", bbox_to_anchor=(1.1, 1.175), edgecolor='black')

    else:
        fig = plt.figure(figsize=(20, 10))
        ax1 = fig.add_subplot(111)
        length = np.arange(0, len(count_sorted))
        ax1.set_title(main)
        ax1.bar(length, count_sorted, color=cs)
        ax1.set_xticks(length)
        ax1.set_xticklabels(barcodes)
        ax1.legend(labels=['%s, %1.1f %%' % (l, s) for l, s in zip(barcodes, sizes)],
                   loc="upper right", bbox_to_anchor=(1.1, 1.175))

    plt.savefig(output_file)
    plt.close()

    barcode_table = pd.DataFrame({"barcode arrangement": count_sorted/sum(count_sorted)*100,
                                  "read count": count_sorted})
    barcode_table.sort_index(inplace=True)
    pd.options.display.float_format = '{:.2f}%'.format

    table_html = pd.DataFrame.to_html(barcode_table)

    return main, output_file, table_html, desc


def barcode_length_boxplot(result_dict, main, my_dpi, result_directory, desc):
    """
    Plot boxplot of the 1D pass and fail read length for each barcode indicated in the sample sheet
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'

    plt.figure(figsize=(12, 7), dpi=my_dpi)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    plt.subplot(gs[0])
    plt.subplots_adjust(bottom=0.015, top=1.0)

    ax = sns.boxplot(data=result_dict['basecaller.sequencing.summary.1d.extractor.barcode_selection_sequence_length_melted_dataframe'],
                     x='barcodes', y='length', hue='passes_filtering',
                     showfliers=False, palette={True: "yellowgreen", False: "orangered"},
                     hue_order=[True, False])

    handles, _ = ax.get_legend_handles_labels()
    plt.legend(bbox_to_anchor=(0.905, 0.98), loc=2, borderaxespad=0., labels=["Pass", "Fail"], handles=handles)
    plt.xlabel('Barcodes')
    plt.ylabel('Read length(bp)')
    plt.title(main)

    df = result_dict['basecaller.sequencing.summary.1d.extractor.barcode_selection_sequence_length_dataframe']
    all_read = df.describe().T
    read_pass = df.loc[df['passes_filtering'] == bool(True)].describe().T
    read_fail = df.loc[df['passes_filtering'] == bool(False)].describe().T
    concat = pd.concat([all_read, read_pass, read_fail], keys=['1D', '1D pass', '1D fail'])
    dataframe = concat.T

    dataframe.loc['count'] = dataframe.loc['count'].astype(int).astype(str)
    dataframe.iloc[1:] = dataframe.iloc[1:].applymap('{:.2f}'.format)
    table_html = pd.DataFrame.to_html(dataframe)

    plt.savefig(output_file)
    plt.close()

    return main, output_file, table_html, desc


def barcoded_phred_score_frequency(result_dict, main, my_dpi, result_directory, desc):
    """
    Plot boxplot of the 1D pass and fail read qscore for each barcode indicated in the sample sheet
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'

    plt.figure(figsize=(12, 7), dpi=my_dpi)
    plt.subplots_adjust(bottom=0.015, top=1.0)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    plt.subplot(gs[0])

    ax = sns.boxplot(data=result_dict['basecaller.sequencing.summary.1d.extractor.barcode_selection_sequence_phred_melted_dataframe'],
                     x='barcodes', y='qscore', hue='passes_filtering', showfliers=False,
                     palette={True: "yellowgreen", False: "orangered"}, hue_order=[True, False])
    handles, _ = ax.get_legend_handles_labels()
    plt.legend(bbox_to_anchor=(0.905, 0.98), loc=2, borderaxespad=0., labels=["Pass", "Fail"], handles=handles)
    plt.xlabel('Barcodes')
    plt.ylabel('Mean Phred score')
    plt.title(main)

    df = result_dict['basecaller.sequencing.summary.1d.extractor.barcode_selection_sequence_phred_dataframe']
    all_read = df.describe().T
    read_pass = df.loc[df['passes_filtering'] == bool(True)].describe().T
    read_fail = df.loc[df['passes_filtering'] == bool(False)].describe().T
    concat = pd.concat([all_read, read_pass, read_fail], keys=['1D', '1D pass', '1D fail'])
    dataframe = concat.T
    dataframe.loc['count'] = dataframe.loc['count'].astype(int).astype(str)
    dataframe.iloc[1:] = dataframe.iloc[1:].applymap('{:.2f}'.format)
    table_html = pd.DataFrame.to_html(dataframe)

    plt.savefig(output_file)
    plt.close()

    return main, output_file, table_html, desc

#
#  1Dsquare plots
#


def dsqr_read_count_histogram(result_dict, main, my_dpi, result_directory, desc):
    """
    Plots the count histograms of count  of the different
    types of reads: 1D, 1D square, 1D square pass, 1D square fail.
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'

    # Histogram completed with the number of basecalling reads (is.barcode == True)
    if 'basecaller.sequencing.summary.1dsqr.extractor.read.pass.barcoded.count' in result_dict or result_dict.keys().__contains__('basecaller.sequencing.summary.1dsqr.extractor.read.pass.barcoded.count'):

        plt.figure(figsize=(13, 7), dpi=my_dpi)

        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
        ax = plt.subplot(gs[0])
        plt.mec = 'black'
        plt.mfc = 'white'
        plt.subplots_adjust(bottom=0.015, top=1.0)

        read_type = [result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
                     result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.count'],
                     result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.count'],
                     result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.barcoded.count'],
                     result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.fail.count'],
                     result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.fail.barcoded.count']]

        label = ("1D", "1Dsquare", "1Dsquare pass", "1Dsquare pass barcoded", "1Dsquare fail","1Dsquare fail barcoded")
        nd = np.arange(len(read_type))
        bars = ax.bar(nd, read_type, align='center', color=["salmon", "orange", "yellowgreen", "lightgoldenrodyellow", "orangered","darksalmon"],
                      edgecolor="black", linewidth=1)

        array = \
            np.array([[result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
                       result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.count'],
                       result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.count'],
                       result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.barcoded.count'],
                       result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.count"],
                       result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.fail.barcoded.count']],
                      [result_dict['basecaller.sequencing.summary.1d.extractor.read.count.frequency'],
                       result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.count.frequency"],
                       result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.frequency"],
                       result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.barcoded.frequency"],
                       result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.frequency"],
                       result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.barcoded.frequency"]]])

        dataframe = pd.DataFrame(array, index=['count', 'frequency'],
                                 columns=["FastQ_entries", "1D square", "1D square pass","1D square pass barcoded", "1D square fail" , "1D square fail barcoded"])

    else:
        plt.figure(figsize=(12, 7), dpi=my_dpi)

        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
        ax = plt.subplot(gs[0])
        plt.mec = 'black'
        plt.mfc = 'white'
        plt.subplots_adjust(bottom=0.015, top=1.0)

        read_type = [result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
                     result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.count'],
                     result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.count'],
                     result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.fail.count']]
        label = ("1D", "1Dsquare", "1Dsquare pass", "1Dsquare fail")
        nd = np.arange(len(read_type))
        bars = ax.bar(nd, read_type, align='center', color=["salmon", "orange", "yellowgreen", "orangered"],
                      edgecolor="black", linewidth=1)

        array = \
            np.array([[result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
                       result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.count'],
                       result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.count'],
                       result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.count"]],
                      [result_dict['basecaller.sequencing.summary.1d.extractor.read.count.frequency'],
                       result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.count.frequency"],
                       result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.frequency"],
                       result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.frequency"]]])

        dataframe = pd.DataFrame(array, index=['count', 'frequency'],
                                 columns=["FastQ_entries", "1D", "1D pass", "1D fail"])

    plt.xticks(nd, label)
    plt.xlabel("Read type")
    plt.ylabel("Counts")
    plt.title(main)

    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2., 1 * height, '%d' % int(height), ha='center', va='bottom')

    plt.savefig(output_file)
    plt.close()

    dataframe.iloc[0] = dataframe.iloc[0].astype(int).astype(str)
    dataframe.iloc[1:] = dataframe.iloc[1:].applymap('{:.2f}'.format)
    table_html = pd.DataFrame.to_html(dataframe)

    return main, output_file, table_html, desc


def dsqr_read_length_multihistogram(result_dict, main, my_dpi, result_directory, desc):
    """
    Plots an histogram of the read length for the different types of read:
    1D, 1D square, 1D square pass, 1D square fail
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'

    read_1d = result_dict["basecaller.sequencing.summary.1d.extractor.sequence.length"]
    read_1dsqr = result_dict['basecaller.sequencing.summary.1dsqr.extractor.sequence.length']
    read_pass_1dsqr = result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.length']
    read_fail_1dsqr = result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.fail.length']

    minimum, maximum = \
        min(result_dict["basecaller.sequencing.summary.1d.extractor.sequence.length"]), \
        max(result_dict["basecaller.sequencing.summary.1d.extractor.sequence.length"])
    read_type = ['1D', '1Dsquare', '1Dsquare pass', '1Dsquare fail']

    plt.figure(figsize=(12, 7), dpi=my_dpi)
    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])
    plt.subplots_adjust(bottom=0.015, top=1)

    n, bins, patches = ax.hist([read_1d, read_1dsqr, read_pass_1dsqr, read_fail_1dsqr],
                               color=["salmon", "goldenrod", "yellowgreen", "orangered"],
                               edgecolor='black', label=read_type,
                               bins=2 ** np.linspace(_safe_log(minimum), _safe_log(maximum), 30))
    plt.legend()

    ax.set_xscale('log', basex=2)
    ax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.0f}'))
    ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.0f}'))
    ax.set_xticks(bins)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_rotation('vertical')
    ax.set_xlabel('Read length(bp)')
    ax.set_ylabel('Read number')
    ax.set_title(main)

    dataframe = \
        pd.DataFrame({"1D": read_1d, '1Dsquare': read_1dsqr,
                      "1Dsquare pass": read_pass_1dsqr,
                      "1Dsquare fail": read_fail_1dsqr})

    dataframe = dataframe[["1D", "1Dsquare", "1Dsquare pass", "1Dsquare fail"]]
    table_html = pd.DataFrame.to_html(_make_desribe_dataframe(dataframe))

    plt.savefig(output_file)
    plt.close()

    return main, output_file, table_html, desc


def dsqr_read_quality_multiboxplot(result_dict, main, my_dpi, result_directory, desc):
    """
    Plots a boxplot of reads quality per read type (1D square, 1D square pass, 1D square fail)
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[2, 1])
    plt.subplots_adjust(bottom=0.015, top=1.0)
    mean_qscore_1d = result_dict["basecaller.sequencing.summary.1d.extractor.mean.qscore"]

    mean_qscore_1dsqr = result_dict['basecaller.sequencing.summary.1dsqr.extractor.mean.qscore']
    read_pass_1dsqr = result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.qscore']
    read_fail_1dsqr = result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.fail.qscore']

    my_pal = {"1D": "salmon",
              "1Dsquare": "goldenrod",
              "1Dsquare pass": "yellowgreen",
              "1Dsquare fail": "orangered"}
    order = ["1D", "1Dsquare", "1Dsquare pass", "1Dsquare fail"]
    dataframe = \
        pd.DataFrame({"1D": mean_qscore_1d,
                      "1Dsquare": mean_qscore_1dsqr,
                      "1Dsquare pass": read_pass_1dsqr,
                      "1Dsquare fail": read_fail_1dsqr})

    ax = plt.subplot(gs[0])
    sns.boxplot(data=dataframe, ax=ax, palette=my_pal, order=order, linewidth=1)
    plt.ylabel('Mean Phred score')
    ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.0f}'))

    ax2 = plt.subplot(gs[1])
    sns.violinplot(data=dataframe, ax=ax2, palette=my_pal, inner=None, cut=0, order=order, linewidth=1)
    plt.subplots_adjust(bottom=0.015, top=1.0)
    plt.ylabel('Mean Phred score')
    ax2.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.0f}'))
    plt.title(main)

    dataframe = dataframe[["1D", "1Dsquare", "1Dsquare pass", "1Dsquare fail"]]
    table_html = pd.DataFrame.to_html(_make_desribe_dataframe(dataframe))

    plt.savefig(output_file)
    plt.close()

    return main, output_file, table_html, desc


def dsqr_phred_score_frequency(result_dict, main, my_dpi, result_directory, desc):
    """
    Plot the distribution of the phred score
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    plt.subplots_adjust(bottom=0.015, top=1.0)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    plt.subplot(gs[0])

    sns.distplot(result_dict['basecaller.sequencing.summary.1dsqr.extractor.mean.qscore'], bins=15, color='goldenrod',
                 hist_kws=dict(edgecolor="k", linewidth=1), hist=True, label="1Dsquare")

    plt.legend()
    plt.xlabel("Mean Phred score")
    plt.ylabel("Frequency")
    plt.title(main)

    dataframe = pd.DataFrame({"1Dsquare": result_dict['basecaller.sequencing.summary.1dsqr.extractor.mean.qscore']})
    rd = dataframe.describe().drop('count').round(2).reset_index()

    plt.axvline(x=result_dict['basecaller.sequencing.summary.1dsqr.extractor.mean.qscore'].describe()['50%'], color='goldenrod')

    plt.savefig(output_file)
    plt.close()

    table_html = pd.DataFrame.to_html(_make_desribe_dataframe(rd))

    return main, output_file, table_html, desc


def dsqr_allphred_score_frequency(result_dict, main, my_dpi, result_directory, desc):
    """
    Plot the distribution of the phred score per read type (1D square , 1D square pass, 1D square fail)
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    plt.subplots_adjust(bottom=0.015, top=1.0)
    gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[2, 1])
    qscore_1d = result_dict["basecaller.sequencing.summary.1d.extractor.mean.qscore"]
    qscore_1dsqr = result_dict['basecaller.sequencing.summary.1dsqr.extractor.mean.qscore']
    mean_qscore_read_pass = result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.qscore']
    mean_qscore_read_fail = result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.fail.qscore']

    ax = plt.subplot(gs[0])
    sns.distplot(result_dict['basecaller.sequencing.summary.1dsqr.extractor.mean.qscore'], ax=ax, bins=15, color='goldenrod',
                 hist_kws=dict(edgecolor="k", linewidth=1), hist=True, label="1Dsquare")

    ax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.0f}'))
    ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.2f}'))
    plt.axvline(x=result_dict['basecaller.sequencing.summary.1dsqr.extractor.mean.qscore'].describe()['50%'], color='goldenrod')
    plt.legend()
    plt.xlabel("Mean Phred score")
    plt.ylabel("Frequency")
    plt.title(main)

    ax2 = plt.subplot(gs[1])

    sns.distplot(mean_qscore_read_pass, ax=ax2, bins=15, hist_kws=dict(edgecolor="k", linewidth=1), color='yellowgreen',
                 hist=True, label='1Dsquare pass')
    sns.distplot(mean_qscore_read_fail, ax=ax2, bins=15, hist_kws=dict(edgecolor="k", linewidth=1), color='orangered',
                 label='1Dsquare fail', hist=True)

    plt.axvline(x=mean_qscore_read_pass.describe()['50%'], color='yellowgreen')
    plt.axvline(x=mean_qscore_read_fail.describe()['50%'], color='orangered')
    ax2.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.0f}'))
    ax2.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.2f}'))

    plt.legend()
    plt.xlabel("Mean Phred score")
    plt.ylabel("Frequency")
    plt.title(main)

    dataframe = \
        pd.DataFrame({"1D": qscore_1d,
                      "1Dsquare": qscore_1dsqr,
                      "1Dsquare pass": mean_qscore_read_pass,
                      "1Dsquare fail": mean_qscore_read_fail})

    dataframe = _make_desribe_dataframe(dataframe).drop('count')

    table_html = pd.DataFrame.to_html(dataframe)

    plt.savefig(output_file)
    plt.close()

    return main, output_file, table_html, desc


def scatterplot_1dsqr(result_dict, main, my_dpi, result_directory, desc):
    """
    Plot the scatter plot representing the relation between the phred score and the sequence length converted in log
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'
    plt.figure(figsize=(12, 7), dpi=my_dpi)
    ax = plt.gca()

    length_1dsqr_read_pass = result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.length']
    length_1dsqr_read_fail = result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.fail.length']

    qscore_1dsqr_read_pass = result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.qscore']
    qscore_1dsqr_read_fail = result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.fail.qscore']

    read_pass = plt.scatter(x=length_1dsqr_read_pass, y=qscore_1dsqr_read_pass, color="yellowgreen")
    read_fail = plt.scatter(x=length_1dsqr_read_fail, y=qscore_1dsqr_read_fail, color="orangered")
    plt.legend((read_pass, read_fail), ("1Dsquare pass", "1Dsquare fail"))
    plt.xlim(np.min(result_dict['basecaller.sequencing.summary.1dsqr.extractor.sequence.length']
                    .loc[result_dict['basecaller.sequencing.summary.1dsqr.extractor.sequence.length'] > 0]),
             np.max(result_dict['basecaller.sequencing.summary.1dsqr.extractor.sequence.length']))
    plt.xscale('log')
    plt.xlabel("Sequence length")
    plt.ylabel("Mean Phred score")
    ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.0f}'))
    plt.title(main)
    plt.savefig(output_file)
    plt.close()

    table_html = None
    return main, output_file, table_html, desc

#
# For eache barcode 1Dsqr
#


def barcode_percentage_pie_chart_1dsqr_pass(result_dict, main, barcode_selection, my_dpi, result_directory, desc):
    """
    Plots a pie chart of 1D square read pass percentage per barcode of a run.
    Needs the sample sheet file describing the barcodes to run.
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'
    plt.figure(figsize=(800 / my_dpi, 800 / my_dpi), dpi=my_dpi)
    for element in barcode_selection:

        if all(result_dict['basecaller.sequencing.summary.1dsqr.extractor.barcode.arrangement'] != element):
            print("The barcode {} doesn't exist".format(element))
            return False

    count_sorted = result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.barcoded"]
    barcodes = count_sorted.index.values.tolist()

    cs = plt.get_cmap('Spectral')(np.arange(len(barcodes)) / len(barcodes))

    sizes = [(100 * chiffre) / sum(count_sorted) for chiffre in count_sorted.values]
    if len(barcode_selection) <= 10:
        fig1, ax1 = plt.subplots()
        ax1.pie(sizes, labels=None, startangle=90, colors=cs, wedgeprops={'linewidth': 1, 'edgecolor': 'k'})
        ax1.axis('equal')
        ax1.legend(labels=['%s, %1.1f %%' % (l, s) for l, s in zip(barcodes, sizes)],
                   loc="upper right", bbox_to_anchor=(1.1, 1.175), edgecolor="black")

    else:
        fig = plt.figure(figsize=(20, 10))
        ax1 = fig.add_subplot(111)
        length = np.arange(0, len(count_sorted))
        ax1.set_title(main)
        ax1.bar(length, count_sorted, color=cs)
        ax1.set_xticks(length)
        ax1.set_xticklabels(barcodes)
        plt.legend(labels=['%s, %1.1f %%' % (l, s) for l, s in zip(barcodes, sizes)],
                   loc="upper right", bbox_to_anchor=(1.1, 1.175))

    plt.savefig(output_file)
    plt.close()
    barcode_table = pd.DataFrame({"barcode arrangement": count_sorted/sum(count_sorted)*100,
                                  "read count": count_sorted})
    barcode_table.sort_index(inplace=True)
    pd.options.display.float_format = '{:.2f}%'.format
    table_html = \
        pd.DataFrame.to_html(barcode_table)

    return main, output_file, table_html, desc


def barcode_percentage_pie_chart_1dsqr_fail(result_dict, main, barcode_selection, my_dpi, result_directory, desc):
    """
    Plots a pie chart of 1D square read fail percentage per barcode of a run.
    Needs the sample sheet file describing the barcodes to run
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'
    plt.figure(figsize=(800 / my_dpi, 800 / my_dpi), dpi=my_dpi)
    for element in barcode_selection:

        if all(result_dict['basecaller.sequencing.summary.1dsqr.extractor.barcode.arrangement'] != element):
            print("The barcode {} doesn't exist".format(element))

    count_sorted = result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.barcoded"]
    barcodes = count_sorted.index.values.tolist()

    cs = plt.get_cmap('Spectral')(np.arange(len(barcodes)) / len(barcodes))

    sizes = [(100 * chiffre) / sum(count_sorted) for chiffre in count_sorted.values]
    if len(barcode_selection) <= 10:
        fig1, ax1 = plt.subplots()
        ax1.pie(sizes, labels=None, startangle=90, colors=cs, wedgeprops={'linewidth': 1, 'edgecolor': 'k'})
        ax1.axis('equal')
        ax1.legend(labels=['%s, %1.1f %%' % (l, s) for l, s in zip(barcodes, sizes)],
                   loc="upper right", bbox_to_anchor=(1.1, 1.175), edgecolor="black")

    else:
        fig = plt.figure(figsize=(20, 10))
        ax1 = fig.add_subplot(111)
        length = np.arange(0, len(count_sorted))
        ax1.set_title(main)
        ax1.bar(length, count_sorted, color=cs)
        ax1.set_xticks(length)
        ax1.set_xticklabels(barcodes)
        plt.legend(labels=['%s, %1.1f %%' % (l, s) for l, s in zip(barcodes, sizes)],
                   loc="upper right", bbox_to_anchor=(1.1, 1.175))

    plt.savefig(output_file)
    plt.close()
    barcode_table = pd.DataFrame({"barcode arrangement": count_sorted/sum(count_sorted)*100,
                                  "read count": count_sorted})
    barcode_table.sort_index(inplace=True)
    pd.options.display.float_format = '{:.2f}%'.format

    table_html = \
        pd.DataFrame.to_html(barcode_table)

    return main, output_file, table_html, desc


def barcode_length_boxplot_1dsqr(result_dict, main, my_dpi, result_directory, desc):
    """
    Plot boxplot of the 1D square pass and fail read length for each barcode indicated in the sample sheet
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'

    plt.figure(figsize=(12, 7), dpi=my_dpi)
    plt.subplots_adjust(bottom=0.015, top=1.0)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    plt.subplot(gs[0])

    ax = \
        sns.boxplot(
            data=result_dict['basecaller.sequencing.summary.1dsqr.extractor.barcode_selection_sequence_length_melted_dataframe'],
            x='barcodes', y='length', hue='passes_filtering', showfliers=False,
            palette={True: "yellowgreen", False: "orangered"}, hue_order=[True, False])
    handles, _ = ax.get_legend_handles_labels()
    plt.legend(bbox_to_anchor=(0.905, 0.98), loc=2, borderaxespad=0., labels=["Pass", "Fail"], handles=handles)
    plt.xlabel('Barcodes')
    plt.ylabel('Read length(bp)')
    plt.title(main)

    df = result_dict['basecaller.sequencing.summary.1dsqr.extractor.barcode_selection_sequence_length_dataframe']
    all_read = df.describe().T
    read_pass = df.loc[df['passes_filtering'] == bool(True)].describe().T
    read_fail = df.loc[df['passes_filtering'] == bool(False)].describe().T
    concat = pd.concat([all_read, read_pass, read_fail], keys=['1Dsqr', '1Dsqr pass', '1Dsqr fail'])
    dataframe = concat.T

    dataframe.loc['count'] = dataframe.loc['count'].astype(int).astype(str)
    dataframe.iloc[1:] = dataframe.iloc[1:].applymap('{:.2f}'.format)
    table_html = pd.DataFrame.to_html(dataframe)

    plt.savefig(output_file)
    plt.close()
    return main, output_file, table_html, desc


def barcoded_phred_score_frequency_1dsqr(result_dict, main, my_dpi, result_directory, desc):
    """
    Plot boxplot of the 1D square pass and fail read qscore for each barcode indicated in the sample sheet
    """
    output_file = result_directory + '/' + '_'.join(main.split()) + '.png'

    plt.figure(figsize=(12, 7), dpi=my_dpi)
    plt.subplots_adjust(bottom=0.015, top=1.0)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    plt.subplot(gs[0])

    ax = \
        sns.boxplot(
            data=result_dict['basecaller.sequencing.summary.1dsqr.extractor.barcode_selection_sequence_phred_melted_dataframe'],
            x='barcodes', y='qscore', hue='passes_filtering', showfliers=False,
            palette={True: "yellowgreen", False: "orangered"}, hue_order=[True, False])
    handles, _ = ax.get_legend_handles_labels()
    plt.legend(bbox_to_anchor=(0.905, 0.98), loc=2, borderaxespad=0., labels=["Pass", "Fail"], handles=handles)
    plt.xlabel('Barcodes')
    plt.ylabel('Mean Phred score')
    plt.title(main)

    df = result_dict['basecaller.sequencing.summary.1dsqr.extractor.barcode_selection_sequence_phred_dataframe']
    all_read = df.describe().T
    read_pass = df.loc[df['passes_filtering'] == bool(True)].describe().T
    read_fail = df.loc[df['passes_filtering'] == bool(False)].describe().T
    concat = pd.concat([all_read, read_pass, read_fail], keys=['1Dsquare', '1Dsquare pass', '1Dsquare fail'])
    dataframe = concat.T

    dataframe.loc['count'] = dataframe.loc['count'].astype(int).astype(str)
    dataframe.iloc[1:] = dataframe.iloc[1:].applymap('{:.2f}'.format)
    table_html = pd.DataFrame.to_html(dataframe)

    plt.savefig(output_file)
    plt.close()

    return main, output_file, table_html, desc
