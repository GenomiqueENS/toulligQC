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

# Class for generating Plotly and MPL graphs and statistics tables in HTML format, they use the result_dict or dataframe_dict dictionnaries.

import copy
import numpy as np
import pandas as pd
import plotly.graph_objs as go

from toulligqc.plotly_graph_common import _barcode_boxplot_graph
from toulligqc.plotly_graph_common import _create_and_save_div
from toulligqc.plotly_graph_common import _dataFrame_to_html
from toulligqc.plotly_graph_common import _format_float
from toulligqc.plotly_graph_common import _format_int
from toulligqc.plotly_graph_common import _legend
from toulligqc.plotly_graph_common import _over_time_graph
from toulligqc.plotly_graph_common import _phred_score_density
from toulligqc.plotly_graph_common import _pie_chart_graph
from toulligqc.plotly_graph_common import _quality_multiboxplot
from toulligqc.plotly_graph_common import _read_length_distribution
from toulligqc.plotly_graph_common import _scatterplot
from toulligqc.plotly_graph_common import _smooth_data
from toulligqc.plotly_graph_common import _title
from toulligqc.plotly_graph_common import _transparent_colors
from toulligqc.plotly_graph_common import _xaxis
from toulligqc.plotly_graph_common import _yaxis
from toulligqc.plotly_graph_common import default_graph_layout
from toulligqc.plotly_graph_common import interpolation_points
from toulligqc.plotly_graph_common import line_width
from toulligqc.plotly_graph_common import plotly_background_color
from toulligqc.plotly_graph_common import toulligqc_colors


#
#  1D plots
#


def read_count_histogram(result_dict, result_directory):
    """
    Plots the histogram of count of the different types of reads:
    1D read return by Guppy
    1D pass read return by Guppy (Qscore >= 7)
    1D fail read return by Guppy (Qscore < 7)
    """

    graph_name = 'Read count histogram'

    # Histogram with barcoded read counts
    if 'basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count' in result_dict:

        data = {
            'All reads': result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
            'Pass reads': result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
            'Fail reads': result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"],
            'Pass barcoded reads': result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count"],
            'Fail barcoded reads': result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.count"]
        }

        colors = [toulligqc_colors["all"], toulligqc_colors["pass"], toulligqc_colors["fail"],
                  toulligqc_colors["barcode_pass"], toulligqc_colors["barcode_fail"]]

        trace = go.Bar(x=[*data], y=list(data.values()),
                       hovertemplate='<b>%{x}</b><br>%{y:,}<extra></extra>',
                       marker_color=_transparent_colors(colors, plotly_background_color, .5),
                       marker_line_color=colors,
                       marker_line_width=line_width)

        # Array of data for HTML table with barcode reads
        array = np.array(
            # count
            [[result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.count"]],
             # frequencies
             [result_dict["basecaller.sequencing.summary.1d.extractor.read.count.frequency"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.frequency"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.frequency"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.frequency"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.frequency"]]])

        dataframe = pd.DataFrame(array, index=['count', 'percent'],
                                 columns=["All reads", "Pass reads", "Fail reads", "Pass barcoded reads",
                                          "Fail barcoded reads"])

    # Histogram without barcodes
    else:

        data = {
            'All reads': result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
            'Pass reads': result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
            'Fail reads': result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"]
        }

        colors = [toulligqc_colors['all'], toulligqc_colors['pass'], toulligqc_colors['fail']]

        trace = go.Bar(x=[*data], y=list(data.values()),
                       hovertemplate='<b>%{x}</b><br>%{y:,}<extra></extra>',
                       marker_color=_transparent_colors(colors, plotly_background_color, .5),
                       marker_line_color=colors,
                       marker_line_width=line_width)

        # Array of data for HTML table without barcode reads
        array = np.array([[result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
                           result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
                           result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"]],
                          # frequencies
                          [result_dict["basecaller.sequencing.summary.1d.extractor.read.count.frequency"],
                           result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.frequency"],
                           result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.frequency"]]])

        # Create dataframe with array data
        dataframe = pd.DataFrame(array, index=['count', 'percent'],
                                 columns=["All reads", "Pass reads", "Fail reads"])

    layout = go.Layout(
        **_title(graph_name),
        **default_graph_layout,
        hovermode="x",
        **_xaxis('Read type', dict(fixedrange=True, categoryorder="total descending")),
        **_yaxis('Read count'))

    fig = go.Figure(data=trace, layout=layout)

    # HTML table
    dataframe.iloc[0] = dataframe.iloc[0].astype(int).apply(lambda x: _format_int(x))
    dataframe.iloc[1:] = dataframe.iloc[1:].applymap(_format_float)
    table_html = _dataFrame_to_html(dataframe)

    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def read_length_scatterplot(dataframe_dict, result_directory):
    graph_name = "Distribution of read lengths"

    return _read_length_distribution(graph_name=graph_name,
                                     all_reads=dataframe_dict['all.reads.sequence.length'],
                                     pass_reads=dataframe_dict['pass.reads.sequence.length'],
                                     fail_reads=dataframe_dict['fail.reads.sequence.length'],
                                     all_color=toulligqc_colors['all'],
                                     pass_color=toulligqc_colors['pass'],
                                     fail_color=toulligqc_colors['fail'],
                                     xaxis_title='Read length (bp)',
                                     result_directory=result_directory)


def yield_plot(df, result_directory, oneDsquare=False):
    """
    Plots the different reads (1D, 1D pass, 1D fail) produced along the run against the time(in hour)
    """

    graph_name = "Yield plot through time"

    if oneDsquare:
        start_time_column = 'start_time1'
    else:
        start_time_column = 'start_time'

    new_df = df.filter(['sequence_length', start_time_column, 'passes_filtering']).sort_values(by=start_time_column)
    new_df['start_time'] = new_df[start_time_column] / 3600

    all_reads_length_df = new_df.filter(['sequence_length', start_time_column])
    pass_reads_length_df = new_df[new_df['passes_filtering'] == True].filter(['sequence_length', start_time_column])
    fail_reads_length_df = new_df[new_df['passes_filtering'] == False].filter(['sequence_length', start_time_column])

    data = [(all_reads_length_df, 'All reads', toulligqc_colors['all']),
            (pass_reads_length_df, 'Pass reads', toulligqc_colors['pass']),
            (fail_reads_length_df, 'Fail reads', toulligqc_colors['fail'])]

    npoints, sigma = interpolation_points(new_df['start_time'], 'yield_plot')
    coef = max(all_reads_length_df[start_time_column]) / npoints

    fig = go.Figure()

    # Figures for cumulative base yield plot
    first = True
    for reads in [True, False]:
        smooth_data_dict = {}
        for d in data:

            if d[1] not in smooth_data_dict:
                if reads:
                    count_x, count_y, cum_count_y = _smooth_data(npoints=npoints, sigma=sigma,
                                                                 data=d[0][start_time_column])
                else:
                    count_x, count_y, cum_count_y = _smooth_data(npoints=npoints, sigma=sigma,
                                                                 data=d[0][start_time_column],
                                                                 weights=d[0]['sequence_length'])

                smooth_data_dict[d[1]] = (count_x, count_y, cum_count_y)

            count_x, count_y, cum_count_y = smooth_data_dict[d[1]]

            fig.add_trace(go.Scatter(x=count_x,
                                     y=cum_count_y,
                                     name=d[1],
                                     fill='tozeroy',
                                     marker_color=d[2],
                                     visible=first
                                     ))

        count_x, count_y, cum_count_y = smooth_data_dict[data[0][1]]
        for p in [50, 75, 90, 99]:
            y = cum_count_y
            ymax = max(y)
            index = (np.abs(y - ymax * p / 100)).argmin()
            x0 = count_x[index]
            fig.add_trace(go.Scatter(
                mode="lines+text",
                name=d[1],
                x=[x0, x0],
                y=[0, ymax],
                line=dict(color="gray", width=1, dash="dot"),
                text=["", str(p) + "% all reads"],
                textposition="top center",
                hoverinfo="skip",
                showlegend=False,
                visible=first
            ))
        first = False

        for d in data:
            count_x, count_y, cum_count_y = smooth_data_dict[d[1]]

            fig.add_trace(go.Scatter(x=count_x,
                                     y=count_y / coef,
                                     name=d[1],
                                     marker_color=d[2],
                                     fill='tozeroy',
                                     visible=False
                                     ))

    fig.update_layout(
        **_title(graph_name),
        **default_graph_layout,
        **_legend(args=dict(y=0.75)),
        hovermode='x',
        **_xaxis('Time (hours)', dict(rangemode="tozero")),
        **_yaxis('Read count', dict(fixedrange=False, rangemode="tozero")),
    )

    # Add buttons

    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="down",
                buttons=list([
                    dict(
                        args=[{'visible': [True, True, True,
                                           True, True, True, True,
                                           False, False, False,
                                           False, False, False,
                                           False, False, False, False,
                                           False, False, False]},
                              {'yaxis': {'title': '<b>Read count</b>', 'rangemode': "tozero"}}],
                        label="Cumulative reads",
                        method="update"
                    ),
                    dict(
                        args=[{'visible': [False, False, False,
                                           False, False, False, False,
                                           True, True, True,
                                           False, False, False,
                                           False, False, False, False,
                                           False, False, False
                                           ]},
                              {'yaxis': {'title': '<b>Read count per hour</b>', 'rangemode': "tozero"}}],
                        label="Yield reads",
                        method="update"
                    ),
                    dict(
                        args=[{'visible': [False, False, False,
                                           False, False, False, False,
                                           False, False, False,
                                           True, True, True,
                                           True, True, True, True,
                                           False, False, False]},
                              {'yaxis': {'title': '<b>Base count</b>', 'rangemode': "tozero"}}],
                        label="Cumulative bases",
                        method="update"
                    ),
                    dict(
                        args=[{'visible': [False, False, False,
                                           False, False, False, False,
                                           False, False, False,
                                           False, False, False,
                                           False, False, False, False,
                                           True, True, True]},
                              {'yaxis': {'title': '<b>Base count per hour</b>', 'rangemode': "tozero"}}],
                        label="Yield bases",
                        method="update"
                    )
                ]),
                pad={"r": 20, "t": 20, "l": 20, "b": 20},
                showactive=True,
                x=1.0,
                xanchor="left",
                y=1.25,
                yanchor="top"
            ),
        ]
    )
    table_html = None
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def read_quality_multiboxplot(dataframe_dict, result_directory):
    """
    Boxplot of PHRED score between read pass and read fail
    Violin plot of PHRED score between read pass and read fail
    """

    graph_name = "PHRED score distribution"

    df = pd.DataFrame(
        {"1D": dataframe_dict['all.reads.mean.qscore'],
         "1D pass": dataframe_dict['pass.reads.mean.qscore'],
         "1D fail": dataframe_dict['fail.reads.mean.qscore']
         })

    return _quality_multiboxplot(graph_name, result_directory, df, onedsquare=False)


def allphred_score_frequency(dataframe_dict, result_directory):
    """
    Plot the distribution of the phred score per read type (1D , 1D pass, 1D fail)
    """

    graph_name = "PHRED score density distribution"

    dataframe = \
        pd.DataFrame({"1D": dataframe_dict['all.reads.mean.qscore'],
                      "1D pass": dataframe_dict['pass.reads.mean.qscore'],
                      "1D fail": dataframe_dict['fail.reads.mean.qscore']})

    return _phred_score_density(graph_name=graph_name,
                                dataframe=dataframe,
                                prefix="1D",
                                all_color=toulligqc_colors['all'],
                                pass_color=toulligqc_colors['pass'],
                                fail_color=toulligqc_colors['fail'],
                                result_directory=result_directory)


def all_scatterplot(dataframe_dict, result_directory):
    """
    Plot the scatter plot representing the relation between the phred score and the sequence length in log
    """

    graph_name = "Correlation between read length and PHRED score"

    return _scatterplot(graph_name, dataframe_dict, result_directory)


def _compute_channel_map(df):

    max_channel = df['channel'].max()
    if max_channel > 512:
        # PromethION geometry
        # Array is simple blocks of 25*10 channels
        channel_array = np.hstack([
            np.arange(x, x + 250).reshape(25, 10)
            for x in range(1, 2752, 250)])
        ca = pd.DataFrame(channel_array).reset_index().melt('index')
        ca.columns = ['row', 'column', 'channel']

    elif max_channel <= 130:
        # Flongle geometry
        # channels are in a simple grid except two upper- and
        # lower-most channels on right-hand column are missing
        channel_array = np.concatenate([
            np.arange(1, 13), np.array([0]),
            np.arange(13, 25), np.array([0]),
            np.arange(25, 115), np.array([0]),
            np.arange(115, 127), np.array([0])]).reshape(10, 13)
        ca = pd.DataFrame(channel_array).reset_index().melt('index')
        ca.columns = ['row', 'column', 'channel']

    else:
        # MinION geometry
        # The array is composed of blocks of 64 channels, the low
        # halve is to the right (high to left), working inwards
        # in a 4 x 8 block
        def channel_block(i):
            m = np.zeros((4, 16), dtype=int)
            m[4::-1, :7:-1] = np.arange(i, i + 32).reshape(4, 8)
            m[4::-1, 0:8] = np.arange(i + 32, i + 64).reshape(4, 8)
            return m

        # the blocks ascend from high to low, except the
        # lowest block is at the top
        first_channel = list(range(65, 450, 64)) + [1]
        channel_array = np.vstack([channel_block(x) for x in first_channel])
        ca = pd.DataFrame(channel_array).reset_index().melt('index')
        ca.columns = ['row', 'column', 'channel']

    return ca


def _compute_channel_count(df, channel_map):

    # count pass reads per channel...
    #counts = df[df['passes_filtering']] \
    counts = df \
        .groupby('channel') \
        .size().to_frame('reads') \
        .reset_index()

    # ...and merge with the channel map
    counts = counts.merge(channel_map, on='channel', how='outer').fillna(0)

    max_row = counts['row'].max()
    max_col = counts['column'].max()
    data = counts.to_dict('split')['data']
    z = []
    ids = []
    for i in range(max_row + 1):
        z.append([np.nan] * (max_col + 1))
        ids.append([0] * (max_col + 1))

    max_value = 0
    for v in data:
        if v[0] == 0:
            z[v[2]][v[3]] = np.nan
            ids[v[2]][v[3]] = 'no channel'
        else:
            z[v[2]][v[3]] = v[1]
            ids[v[2]][v[3]] = str(v[0])
            max_value = max(max_value, v[1])

    return max_row, max_col, int(max_value), counts, z, ids


def plot_performance(df, result_directory):
    """
    Plots the channels occupancy by the reads
    @:param pore_measure: reads number per pore
    """

    graph_name = "Channel occupancy of the flowcell"

    # Compute geometry of the flowcell
    channel_map = _compute_channel_map(df)

    max_row, max_col, max_value, counts, z_pass, ids = _compute_channel_count(df[df['passes_filtering']], channel_map)
    max_row, max_col, max_value, counts, z_fail, ids = _compute_channel_count(df[~df['passes_filtering']], channel_map)
    max_row, max_col, max_value, counts, z_all, ids = _compute_channel_count(df, channel_map)

    # Compute fail ratio
    z_ratio = copy.deepcopy(z_all)
    for i in range(len(z_ratio)):
        for j in range(len(z_ratio[i])):
            if z_ratio[i][j] > 0:
                z_ratio[i][j] = z_fail[i][j] / z_ratio[i][j] * 100.0

    fig = go.Figure()
    fig.add_trace(go.Heatmap(x=list(range(1, max_col + 1)),
                             y=list(range(1, max_row + 1)),
                             z=z_all,
                             zmax=max_value,
                             zmin=0,
                             name='All reads',
                             colorbar=dict(title='Reads', len=0.6, yanchor='middle'),
                             hovertemplate='<b>Channel ID:</b> %{text}<br>'
                                           '<b>Row:</b> %{y}<br>'
                                           '<b>Column:</b> %{x}<br>'
                                           '<b>Reads:</b> %{z}<br>',
                             text=ids,
                             hoverongaps=False,
                             visible=True))

    fig.add_trace(go.Heatmap(x=list(range(1, max_col + 1)),
                             y=list(range(1, max_row + 1)),
                             z=z_pass,
                             zmax=max_value,
                             zmin=0,
                             name='Pass reads',
                             colorbar=dict(title='Reads', len=0.6, yanchor='middle'),
                             hovertemplate='<b>Channel ID:</b> %{text}<br>'
                                           '<b>Row:</b> %{y}<br>'
                                           '<b>Column:</b> %{x}<br>'
                                           '<b>Reads:</b> %{z}<br>',
                             text=ids,
                             hoverongaps=False,
                             visible=False))

    fig.add_trace(go.Heatmap(x=list(range(1, max_col + 1)),
                             y=list(range(1, max_row + 1)),
                             z=z_fail,
                             zmax=max_value,
                             zmin=0,
                             name='Fail reads',
                             colorbar=dict(title='Reads', len=0.6, yanchor='middle'),
                             hovertemplate='<b>Channel ID:</b> %{text}<br>'
                                           '<b>Row:</b> %{y}<br>'
                                           '<b>Column:</b> %{x}<br>'
                                           '<b>Reads:</b> %{z}<br>',
                             text=ids,
                             hoverongaps=False,
                             visible=False))

    fig.add_trace(go.Heatmap(x=list(range(1, max_col + 1)),
                             y=list(range(1, max_row + 1)),
                             z=z_ratio,
                             name='Percent of<br>fail reads',
                             colorscale='reds',
                             colorbar=dict(title='%', len=0.6, yanchor='middle'),
                             hovertemplate='<b>Channel ID:</b> %{text}<br>'
                                           '<b>Row:</b> %{y}<br>'
                                           '<b>Column:</b> %{x}<br>'
                                           '<b>Fail reads:</b> %{z:.1f}%<br>',
                             text=ids,
                             hoverongaps=False,
                             visible=False))

    graph_layout = dict(default_graph_layout)
    graph_layout['plot_bgcolor'] = 'rgba(0,0,0,0)'

    fig.update_layout(
        **_title(graph_name),
        **_legend(),
        **graph_layout,
        hovermode='x',
        **_xaxis('Columns', dict(fixedrange=True)),
        **_yaxis('Rows', dict(fixedrange=True))
    )

    # Add buttons
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="down",
                buttons=list([
                    dict(
                        args=[{'visible': [True, False, False, False]}, {'hovermode': 'x'}],
                        label="All reads",
                        method="update"
                    ),
                    dict(
                        args=[{'visible': [False, True, False, False]}, {'hovermode': 'x'}],
                        label="Pass reads",
                        method="update"
                    ),
                    dict(
                        args=[{'visible': [False, False, True, False]}, {'hovermode': 'x'}],
                        label="Fail reads",
                        method="update"
                    ),
                    dict(
                        args=[{'visible': [False, False, False, True]}, {'hovermode': 'x'}],
                        label="Fail percent",
                        method="update"
                    )
                ]),
                pad={"r": 20, "t": 20, "l": 20, "b": 20},
                showactive=True,
                x=1.0,
                xanchor="left",
                y=1.25,
                yanchor="top"
            ),
        ]
    )

    table_html = None
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


#
# For each barcode 1D
#


def barcode_percentage_pie_chart_pass(dataframe_dict, barcode_selection, result_directory):
    """
    Plots a pie chart of 1D read pass percentage per barcode of a run.
    """

    graph_name = "Pass barcoded reads distribution"

    count_sorted = dataframe_dict["read.pass.barcoded"]

    return _pie_chart_graph(graph_name=graph_name,
                            count_sorted=count_sorted,
                            color_palette=toulligqc_colors['pie_chart_palette'],
                            one_d_square=False,
                            result_directory=result_directory)


def barcode_percentage_pie_chart_fail(dataframe_dict, barcode_selection, result_directory):
    """
    Plots a pie chart of 1D read fail percentage per barcode of a run.
    Needs the samplesheet file describing the barcodes to run
    """

    graph_name = "Fail barcoded reads distribution"

    count_sorted = dataframe_dict["read.fail.barcoded"]

    return _pie_chart_graph(graph_name=graph_name,
                            count_sorted=count_sorted,
                            color_palette=toulligqc_colors['pie_chart_palette'],
                            one_d_square=False,
                            result_directory=result_directory)


def barcode_length_boxplot(datafame_dict, result_directory):
    """
    Boxplots all the 1D pass and fail read length for each barcode indicated in the sample sheet
    """

    graph_name = "Read size distribution for barcodes"

    df = datafame_dict['barcode_selection_sequence_length_dataframe']

    return _barcode_boxplot_graph(graph_name=graph_name,
                                  df=df,
                                  barcode_selection=df.columns.drop('passes_filtering'),
                                  pass_color=toulligqc_colors['pass'],
                                  fail_color=toulligqc_colors['fail'],
                                  yaxis_title="Sequence length (bp)",
                                  legend_title="Read type",
                                  result_directory=result_directory)


def barcoded_phred_score_frequency(dataframe_dict, result_directory):
    """
    Plot boxplot of the 1D pass and fail read qscore for each barcode indicated in the sample sheet
    """

    graph_name = "PHRED score distribution for barcodes"

    df = dataframe_dict['barcode_selection_sequence_phred_dataframe']

    return _barcode_boxplot_graph(graph_name=graph_name,
                                  df=df,
                                  barcode_selection=df.columns.drop('passes_filtering'),
                                  pass_color=toulligqc_colors['pass'],
                                  fail_color=toulligqc_colors['fail'],
                                  yaxis_title="PHRED score",
                                  legend_title="Read type",
                                  result_directory=result_directory)


def sequence_length_over_time(dataframe_dict, result_directory):
    graph_name = "Read length over time"

    return _over_time_graph(data_series=dataframe_dict['all.reads.sequence.length'],
                            time_series=dataframe_dict['all.reads.start.time'],
                            result_directory=result_directory,
                            graph_name=graph_name,
                            color=toulligqc_colors['sequence_length_over_time'],
                            yaxis_title='Read length (bp)')


def phred_score_over_time(dataframe_dict, result_dict, result_directory):
    graph_name = "PHRED score over time"

    pass_min_qscore = 7
    key = 'sequencing.telemetry.extractor.pass.threshold.qscore'
    if key in result_dict:
        pass_min_qscore = float(result_dict[key])

    qscore_series = dataframe_dict["all.reads.mean.qscore"]
    time_series = dataframe_dict['all.reads.start.time']

    return _over_time_graph(data_series=qscore_series,
                            time_series=time_series,
                            result_directory=result_directory,
                            graph_name=graph_name,
                            color=toulligqc_colors['phred_score_over_time'],
                            yaxis_title='PHRED quality score',
                            min_max=True,
                            yaxis_starts_zero=True,
                            green_zone_starts_at=pass_min_qscore,
                            green_zone_color=toulligqc_colors['green_zone_color'])


def speed_over_time(dataframe_dict, result_directory):
    graph_name = "Translocation speed"

    sequence_length_series = dataframe_dict['all.reads.sequence.length']
    duration_series = dataframe_dict['all.reads.duration']

    speed_series = pd.Series(sequence_length_series / duration_series)

    return _over_time_graph(data_series=speed_series,
                            time_series=dataframe_dict['all.reads.start.time'],
                            result_directory=result_directory,
                            graph_name=graph_name,
                            color=toulligqc_colors['speed_over_time'],
                            yaxis_title='Speed (bases per second)',
                            green_zone_starts_at=300,
                            green_zone_color=toulligqc_colors['green_zone_color'])
