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
# Maintainer: Karine Dias
# Since version 0.1

# Class for generating Plotly and MPL graphs and statistics tables in HTML format, they use the result_dict or dataframe_dict_1dsqr dictionnaries.

import numpy as np
import pandas as pd
import plotly.graph_objs as go

from toulligqc.plotly_graph_common import _barcode_boxplot_graph
from toulligqc.plotly_graph_common import _create_and_save_div
from toulligqc.plotly_graph_common import _dataFrame_to_html
from toulligqc.plotly_graph_common import _format_float
from toulligqc.plotly_graph_common import _format_int
from toulligqc.plotly_graph_common import _over_time_graph
from toulligqc.plotly_graph_common import _phred_score_density
from toulligqc.plotly_graph_common import _pie_chart_graph
from toulligqc.plotly_graph_common import _quality_multiboxplot
from toulligqc.plotly_graph_common import _read_length_distribution
from toulligqc.plotly_graph_common import _scatterplot
from toulligqc.plotly_graph_common import _title
from toulligqc.plotly_graph_common import _transparent_colors
from toulligqc.plotly_graph_common import _xaxis
from toulligqc.plotly_graph_common import _yaxis
from toulligqc.plotly_graph_common import default_graph_layout
from toulligqc.plotly_graph_common import line_width
from toulligqc.plotly_graph_common import plotly_background_color
from toulligqc.plotly_graph_common import toulligqc_colors


#
#  1D² plots
#

def dsqr_read_count_histogram(result_dict, result_directory):
    """
    Plots the histogram of 1D² count of the different types of reads:
    1D² read return by Guppy
    1D² pass read return by Guppy (Qscore >= 7)
    1D² fail read return by Guppy (Qscore < 7)
    """

    graph_name = "1D² Read count histogram"

    # Histogram with barcoded read counts
    if 'basecaller.sequencing.summary.1dsqr.extractor.read.pass.barcoded.count' in result_dict:

        data = {
            'All reads': result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
            '1D² reads': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.count"],
            '1D² pass reads': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.count"],
            '1D² fail reads': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.count"],
            '1D² pass barcoded reads': result_dict[
                "basecaller.sequencing.summary.1dsqr.extractor.read.pass.barcoded.count"],
            '1D² fail barcoded reads': result_dict[
                "basecaller.sequencing.summary.1dsqr.extractor.read.fail.barcoded.count"]
        }

        colors = [toulligqc_colors["all"], toulligqc_colors["all_1d2"], toulligqc_colors["pass"],
                  toulligqc_colors["fail"], toulligqc_colors["barcode_pass"], toulligqc_colors["barcode_fail"]]

        trace = go.Bar(x=[*data], y=list(data.values()),
                       hovertemplate='<b>%{x}</b><br>%{y:,}<extra></extra>',
                       marker_color=_transparent_colors(colors, plotly_background_color, .5),
                       marker_line_color=colors,
                       marker_line_width=line_width)

        # Array of data for HTML table with barcode reads
        array = np.array(
            # count
            [[result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.count"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.count"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.count"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.barcoded.count"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.barcoded.count"]],
             # frequencies
             [result_dict["basecaller.sequencing.summary.1d.extractor.read.count.frequency"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.count.frequency"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.frequency"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.frequency"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.barcoded.frequency"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.barcoded.frequency"]]])

        dataframe = pd.DataFrame(array, index=['count', 'percent'],
                                 columns=["All reads", "1D² reads", "1D² pass reads", "1D² fail reads",
                                          "1D² pass barcoded",
                                          "1D² fail barcoded"])

    # Histogram without barcodes
    else:

        data = {
            'All reads': result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
            '1D² reads': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.count"],
            '1D² pass reads': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.count"],
            '1D² fail reads': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.count"]
        }

        colors = [toulligqc_colors["all"], toulligqc_colors["all_1d2"], toulligqc_colors["pass"],
                  toulligqc_colors["fail"]]

        trace = go.Bar(x=[*data], y=list(data.values()),
                       hovertemplate='<b>%{x}</b><br>%{y:,}<extra></extra>',
                       marker_color=_transparent_colors(colors, plotly_background_color, .5),
                       marker_line_color=colors,
                       marker_line_width=line_width)

        # Array of data for HTML table without barcode reads
        array = np.array([[result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
                           result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.count"],
                           result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.count"],
                           result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.count"]],
                          # frequencies
                          [result_dict["basecaller.sequencing.summary.1d.extractor.read.count.frequency"],
                           result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.count.frequency"],
                           result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.frequency"],
                           result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.frequency"]]])

        # Create dataframe with array data
        dataframe = pd.DataFrame(array, index=['count', 'percent'],
                                 columns=["All reads", "1D² reads", "1D² pass reads", "1D² fail reads"])

    layout = go.Layout(
        hovermode="x",
        **_title(graph_name),
        **default_graph_layout,
        **_xaxis('1D² Read type', dict(fixedrange=True, categoryorder="trace")),
        **_yaxis('Read count'))

    fig = go.Figure(data=trace, layout=layout)

    # HTML table
    dataframe.iloc[0] = dataframe.iloc[0].astype(int).apply(lambda x: _format_int(x))
    dataframe.iloc[1:] = dataframe.iloc[1:].applymap(_format_float)
    table_html = _dataFrame_to_html(dataframe)
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def dsqr_read_length_scatterplot(dataframe_dict_1dsqr, result_directory):
    graph_name = "1D² Distribution of read lengths"

    return _read_length_distribution(graph_name=graph_name,
                                     all_reads=dataframe_dict_1dsqr['all.reads.sequence.length'],
                                     pass_reads=dataframe_dict_1dsqr['pass.reads.sequence.length'],
                                     fail_reads=dataframe_dict_1dsqr['fail.reads.sequence.length'],
                                     all_color=toulligqc_colors['all'],
                                     pass_color=toulligqc_colors['pass'],
                                     fail_color=toulligqc_colors['fail'],
                                     xaxis_title='1D² Read length (bp)',
                                     result_directory=result_directory)


def dsqr_read_quality_multiboxplot(result_dict, dataframe_dict_1dsqr, result_directory):
    """
    Boxplot of PHRED score between read pass and read fail
    Violin plot of PHRED score between read pass and read fail
    """

    graph_name = "1D² PHRED score distribution"

    df = pd.DataFrame(
        {"1D²": dataframe_dict_1dsqr["all.reads.mean.qscore"],
         "1D² pass": dataframe_dict_1dsqr['pass.reads.mean.qscore'],
         "1D² fail": dataframe_dict_1dsqr['fail.reads.mean.qscore']
         })

    return _quality_multiboxplot(graph_name, result_directory, df, onedsquare=True)


def dsqr_allphred_score_frequency(result_dict, dataframe_dict_1dsqr, result_directory):
    """
    Plot the distribution of the phred score per read type (1D² , 1D² pass, 1D² fail)
    """

    graph_name = "1D² PHRED score density distribution"

    dataframe = \
        pd.DataFrame({"1D²": dataframe_dict_1dsqr["all.reads.mean.qscore"],
                      "1D² pass": dataframe_dict_1dsqr['pass.reads.mean.qscore'],
                      "1D² fail": dataframe_dict_1dsqr['fail.reads.mean.qscore']})

    return _phred_score_density(graph_name=graph_name,
                                dataframe=dataframe,
                                prefix="1D²",
                                all_color=toulligqc_colors['all'],
                                pass_color=toulligqc_colors['pass'],
                                fail_color=toulligqc_colors['fail'],
                                result_directory=result_directory)


def scatterplot_1dsqr(dataframe_dict_1dsqr, result_directory):
    """
    Plot the scatter plot representing the relation between the phred score and the sequence length in log
    """

    graph_name = "Correlation between 1D² read length and PHRED score"

    return _scatterplot(graph_name, dataframe_dict_1dsqr, result_directory, onedsquare=True)


#
# For each barcode 1D²
#

def barcode_percentage_pie_chart_1dsqr_pass(dataframe_dict_1dsqr, barcode_selection, result_directory):
    """
    Plots a pie chart of 1D² read pass percentage per barcode of a run.
    """

    graph_name = "1D² read pass barcode distribution"

    count_sorted = dataframe_dict_1dsqr["read.pass.barcoded"]

    return _pie_chart_graph(graph_name=graph_name,
                            count_sorted=count_sorted,
                            color_palette=toulligqc_colors['pie_chart_palette'],
                            one_d_square=True,
                            result_directory=result_directory)


def barcode_percentage_pie_chart_1dsqr_fail(dataframe_dict_1dsqr, barcode_selection, result_directory):
    """
    Plots a pie chart of 1D² read fail percentage per barcode of a run.
    Needs the samplesheet file describing the barcodes to run
    """

    graph_name = "1D² read fail barcode distribution"

    count_sorted = dataframe_dict_1dsqr["read.fail.barcoded"]

    return _pie_chart_graph(graph_name=graph_name,
                            count_sorted=count_sorted,
                            color_palette=toulligqc_colors['pie_chart_palette'],
                            one_d_square=True,
                            result_directory=result_directory)


def barcode_length_boxplot_1dsqr(dataframe_dict_1dsqr, result_directory):
    """
    Boxplots all the 1D² pass and fail read length for each barcode indicated in the sample sheet
    """

    graph_name = "1D² read size distribution for barcodes"

    df = dataframe_dict_1dsqr['barcode_selection_sequence_length_dataframe']

    return _barcode_boxplot_graph(graph_name=graph_name,
                                  df=df,
                                  barcode_selection=df.columns.drop('passes_filtering'),
                                  pass_color=toulligqc_colors['pass'],
                                  fail_color=toulligqc_colors['fail'],
                                  yaxis_title='Sequence length (bp)',
                                  legend_title='1D² read type',
                                  result_directory=result_directory)


def barcoded_phred_score_frequency_1dsqr(dataframe_dict_1dsqr, result_directory):
    """
    Plot boxplot of the 1D pass and fail read qscore for each barcode indicated in the sample sheet
    """

    graph_name = "1D² PHRED score distribution for barcodes"

    df = dataframe_dict_1dsqr['barcode_selection_sequence_phred_dataframe']

    return _barcode_boxplot_graph(graph_name=graph_name,
                                  df=df,
                                  barcode_selection=df.columns.drop('passes_filtering'),
                                  pass_color=toulligqc_colors['pass'],
                                  fail_color=toulligqc_colors['fail'],
                                  yaxis_title='PHRED score',
                                  legend_title='1D² read type',
                                  result_directory=result_directory)


def sequence_length_over_time_dsqr(dataframe_dict_1dsqr, result_directory):
    graph_name = "1D² Read length over time"

    return _over_time_graph(data_series=dataframe_dict_1dsqr['all.reads.sequence.length'],
                            time_series=dataframe_dict_1dsqr['all.reads.start.time1'],
                            result_directory=result_directory,
                            graph_name=graph_name,
                            color=toulligqc_colors['sequence_length_over_time'],
                            yaxis_title='Read length (bp)')


def phred_score_over_time_dsqr(result_dict, dataframe_dict_1dsqr, result_directory):
    graph_name = "1D² PHRED score over time"

    pass_min_qscore = 7
    key = 'sequencing.telemetry.extractor.pass.threshold.qscore'
    if key in result_dict:
        pass_min_qscore = float(result_dict[key])

    return _over_time_graph(data_series=dataframe_dict_1dsqr['all.reads.mean.qscore'],
                            time_series=dataframe_dict_1dsqr['all.reads.start.time1'],
                            result_directory=result_directory,
                            graph_name=graph_name,
                            color=toulligqc_colors['phred_score_over_time'],
                            yaxis_title='PHRED quality score',
                            yaxis_starts_zero=True,
                            green_zone_starts_at=pass_min_qscore,
                            green_zone_color=toulligqc_colors['green_zone_color'])


def speed_over_time_dsqr(dataframe_dict_1dsqr, result_directory):
    graph_name = "1D² translocation speed"

    speed = pd.Series(dataframe_dict_1dsqr['all.reads.sequence.length'] / dataframe_dict_1dsqr['all.reads.duration'])

    return _over_time_graph(data_series=speed,
                            time_series=dataframe_dict_1dsqr['all.reads.start.time1'],
                            result_directory=result_directory,
                            graph_name=graph_name,
                            color=toulligqc_colors['speed_over_time'],
                            yaxis_title='Speed (bases per second)',
                            green_zone_starts_at=300,
                            green_zone_color=toulligqc_colors['green_zone_color'])
