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

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import norm
import plotly.graph_objs as go

# Import common function
from toulligqc.plotly_graph_common import _smooth_data
from toulligqc.plotly_graph_common import _precompute_boxplot_values
from toulligqc.plotly_graph_common import _dataFrame_to_html
from toulligqc.plotly_graph_common import _interpolate
from toulligqc.plotly_graph_common import _make_describe_dataframe
from toulligqc.plotly_graph_common import _transparent_colors
from toulligqc.plotly_graph_common import _create_and_save_div

# Import common constants
from toulligqc.plotly_graph_common import figure_image_width
from toulligqc.plotly_graph_common import figure_image_height
from toulligqc.plotly_graph_common import int_format_str
from toulligqc.plotly_graph_common import float_format_str
from toulligqc.plotly_graph_common import percent_format_str
from toulligqc.plotly_graph_common import toulligqc_colors
from toulligqc.plotly_graph_common import plotly_background_color
from toulligqc.plotly_graph_common import line_width
from toulligqc.plotly_graph_common import interpolation_threshold
from toulligqc.plotly_graph_common import legend_font_size
from toulligqc.plotly_graph_common import axis_font_size
from toulligqc.plotly_graph_common import on_chart_font_size
from toulligqc.plotly_graph_common import title_size

#
#  1D plots
#


def read_count_histogram(result_dict, dataframe_dict, main, my_dpi, result_directory, desc):
    """
    Plots the histogram of count of the different types of reads:
    1D read return by Guppy
    1D pass read return by Guppy (Qscore >= 7)
    1D fail read return by Guppy (Qscore < 7)
    """
    # Histogram with barcoded read counts
    if 'read.pass.barcoded.count' in dataframe_dict:

        data = {
            'Read Count': result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
            'Read Pass Count': result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
            'Read Fail Count': result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"],
            'Read Pass Barcoded Count': dataframe_dict["read.pass.barcoded.count"],
            'Read Fail Barcoded Count': dataframe_dict["read.fail.barcoded.count"]
        }

        colors = [toulligqc_colors["all"], toulligqc_colors["pass"], toulligqc_colors["fail"],
                  toulligqc_colors["barcode_pass"], toulligqc_colors["barcode_fail"]]

        trace = go.Bar(x=[*data], y=list(data.values()),
                                hovertext=["<b>Total number of reads</b>",
                                           "<b>Pass reads</b>",
                                           "<b>Fail reads</b>",
                                           "<b>Barcoded pass reads</b>",
                                           "<b>Barcoded fail reads</b>"],
                                name="Barcoded graph",
                                marker_color=_transparent_colors(colors, plotly_background_color, .5),
                                marker_line_color=colors,
                                marker_line_width=line_width)

        # Array of data for HTML table with barcode reads
        array = np.array(
            #count
            [[result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"],
              dataframe_dict["read.pass.barcoded.count"],
              dataframe_dict["read.fail.barcoded.count"]],
             #frequencies
             [result_dict["basecaller.sequencing.summary.1d.extractor.read.count.frequency"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.frequency"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.frequency"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.frequency"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.frequency"]]])

        dataframe = pd.DataFrame(array, index=['count', 'frequency'],
                                 columns=["Read count", "1D pass", "1D fail", "1D pass barcoded", "1D fail barcoded"])

    # Histogram without barcodes
    else:

        data = {
            'Read Count': result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
            'Read Pass Count': result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
            'Read Fail Count': result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"]
        }

        colors = [toulligqc_colors['all'], toulligqc_colors['pass'], toulligqc_colors['fail']]

        trace = go.Bar(x=[*data], y=list(data.values()),
                       hovertext=["<b>Total number of reads</b>",
                                  "<b>Pass reads</b>",
                                  "<b>Fail reads</b>"],
                       name="Barcoded graph",
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
        dataframe = pd.DataFrame(array, index=['count', 'frequency'],
                                 columns=["Read count", "1D pass", "1D fail"])

    layout = go.Layout(
        hovermode="x",
        title={
            'text': "<b>Read count histogram</b>",
            'y': 0.95,
            'x': 0,
            'xanchor': 'left',
            'yanchor': 'top',
            'font': dict(
                size=title_size,
                color="black")
        },
        xaxis=dict(title="<b>Read type</b>",
                   titlefont=dict(
                       size=axis_font_size,
                       color="black"
                   ),
                   categoryorder="total descending"
                   ),
        yaxis=dict(title="<b>Counts</b>",
                   titlefont=dict(
                       size=axis_font_size,
                       color="black"
                   )),
        width=figure_image_width,
        height=figure_image_height)

    fig = go.Figure(data=trace, layout=layout)

    # HTML table
    dataframe.iloc[0] = dataframe.iloc[0].astype(int).apply(lambda x: int_format_str.format(x))
    dataframe.iloc[1:] = dataframe.iloc[1:].applymap(float_format_str.format)
    table_html = _dataFrame_to_html(dataframe)

    div, output_file = _create_and_save_div(fig, result_directory, main)
    return main, output_file, table_html, desc, div


def read_length_scatterplot(result_dict, sequence_length_df, main, my_dpi, result_directory, desc):

    all_read = sequence_length_df.loc[sequence_length_df >= 10].dropna().values
    read_pass = result_dict['basecaller.sequencing.summary.1d.extractor.read.pass.length'].loc[result_dict['basecaller.sequencing.summary.1d.extractor.read.pass.length'] >= 10]
    read_fail = result_dict['basecaller.sequencing.summary.1d.extractor.read.fail.length'].loc[result_dict['basecaller.sequencing.summary.1d.extractor.read.fail.length'] >= 10]

    count_x1, count_y1 = _smooth_data(10000, 5, all_read)
    count_x2, count_y2 = _smooth_data(10000, 5, read_pass)
    count_x3, count_y3 = _smooth_data(10000, 5, read_fail)

    # Find 50 percentile for zoomed range on x axis
    max_x_range = max(np.percentile(count_x1, 50), np.percentile(count_x2, 50), np.percentile(count_x3, 50))

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=count_x1,
                             y=count_y1,
                             name='All reads',
                             hoverinfo='x+y',
                             fill='tozeroy',
                             marker_color=toulligqc_colors['all']
                             ))
    fig.add_trace(go.Scatter(x=count_x2,
                               y=count_y2,
                               name='Pass reads',
                               hoverinfo='x+y',
                               fill='tozeroy',
                               marker_color=toulligqc_colors['pass']
                               ))
    fig.add_trace(go.Scatter(x=count_x3,
                               y=count_y3,
                               name='Fail reads',
                               hoverinfo='x+y',
                               fill='tozeroy',
                               marker_color=toulligqc_colors['fail']
                               ))

    fig.update_layout(
        title={
            'text': "<b>Distribution of read lengths</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        size=title_size,
                        color="black")},
        xaxis=dict(
            title="<b>Read length (bp)</b>",
            titlefont_size=axis_font_size,
            range=[0, max_x_range]
        ),
        yaxis=dict(
            title='<b>Density</b>',
            titlefont_size=axis_font_size,
            tickfont_size=axis_font_size
        ),
        legend=dict(
            x=1.02,
            y=0.95,
            title_text="<b>Legend</b>",
            title=dict(font=dict(size=legend_font_size)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=legend_font_size)
        ),
        hovermode='x',
        height=figure_image_height, width=figure_image_width
    )

    table_html = None
    div, output_file = _create_and_save_div(fig, result_directory, main)
    return main, output_file, table_html, desc, div


def yield_plot(result_dict, main, my_dpi, result_directory, desc):
    """
    Plots the different reads (1D, 1D pass, 1D fail) produced along the run against the time(in hour)
    """
    all_read = result_dict['basecaller.sequencing.summary.1d.extractor.start.time.sorted']
    read_pass = result_dict['basecaller.sequencing.summary.1d.extractor.read.pass.sorted']
    read_fail = result_dict['basecaller.sequencing.summary.1d.extractor.read.fail.sorted']

    count_x1, count_y1 = _smooth_data(10000, 5, all_read)
    count_x2, count_y2 = _smooth_data(10000, 5, read_pass)
    count_x3, count_y3 = _smooth_data(10000, 5, read_fail)

    fig = go.Figure()

    fig.add_trace(go.Scatter(x=count_x1,
                             y=count_y1,
                               name='All reads',
                               marker_color=toulligqc_colors['all'],
                               fill='tozeroy'
                               ))

    fig.add_trace(go.Scatter(x=count_x2,
                             y=count_y2,
                               name='Pass reads',
                               marker_color=toulligqc_colors['pass'],
                               fill='tozeroy'
                               ))

    fig.add_trace(go.Scatter(x=count_x3,
                             y=count_y3,
                               name='Fail reads',
                               marker_color=toulligqc_colors['fail'],
                               fill='tozeroy'
                               ))
    # Figures for cumulative yield plot
    fig.add_trace(go.Scatter(x=count_x1,
                             y=np.cumsum(count_y1),
                             name='All reads',
                             hoverinfo='x+y',
                             fill='tozeroy',
                             marker_color=toulligqc_colors['all'],
                             visible=False
                             ))
    fig.add_trace(go.Scatter(x=count_x2,
                               y=np.cumsum(count_y2),
                               name='Pass reads',
                               hoverinfo='x+y',
                               fill='tozeroy',
                               marker_color=toulligqc_colors['pass'],
                               visible=False
                               ))
    fig.add_trace(go.Scatter(x=count_x3,
                               y=np.cumsum(count_y3),
                               name='Fail reads',
                               hoverinfo='x+y',
                               fill='tozeroy',
                               marker_color=toulligqc_colors['fail'],
                               visible=False
                               ))

    fig.update_layout(
        title={
            'text': "<b>Yield plot through experiment time</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        size=title_size,
                        color="black")},
        xaxis=dict(
            title="<b>Time (hours)</b>",
            titlefont_size=axis_font_size
        ),
        yaxis=dict(
            title='<b>Density</b>',
            titlefont_size=axis_font_size,
            tickfont_size=axis_font_size,
        ),
        legend=dict(
            x=1.02,
            y=0.95,
            title_text="<b>Legend</b>",
            title=dict(font=dict(size=legend_font_size)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=legend_font_size)
        ),
        hovermode='x',
        height=figure_image_height, width=figure_image_width
    )

    # Add buttons
    fig.update_layout(
        updatemenus=[
            dict(
                type = "buttons",
                direction = "left",
                buttons=list([
                    dict(
                        args=[{'visible': [True, True, True, False, False, False]}],
                        label="Yield plot",
                        method="update"
                    ),
                    dict(
                        args=[{'visible': [False, False, False, True, True, True]}],
                        label="Cumulative yield plot",
                        method="update"
                    )
                ]),
                pad={"r": 20, "t": 20, "l":20, "b":20},
                showactive=True,
                x=1.0,
                xanchor="left",
                y=1.25,
                yanchor="top"
            ),
        ]
    )

    table_html = None
    div, output_file = _create_and_save_div(fig, result_directory, main)
    return main, output_file, table_html, desc, div


def read_quality_multiboxplot(result_dict, main, my_dpi, result_directory, desc):
    """
    Boxplot of PHRED score between read pass and read fail
    Violin plot of PHRED score between read pass and read fail
    """
    df = pd.DataFrame(
        {"1D": result_dict["basecaller.sequencing.summary.1d.extractor.mean.qscore"],
         "1D pass": result_dict['basecaller.sequencing.summary.1d.extractor.read.pass.qscore'],
         "1D fail": result_dict['basecaller.sequencing.summary.1d.extractor.read.fail.qscore']
         })

    # If more than 10.000 reads, interpolate data
    if len(df["1D"]) > interpolation_threshold:
        dataframe = pd.DataFrame({
        "1D" : _interpolate(df["1D"], 1000),
        "1D pass" : _interpolate(df["1D pass"], 1000),
        "1D fail" : _interpolate(df["1D fail"], 1000)
    })
    else:
        dataframe = df
    names = {"1D": "All reads",
             "1D pass": "Read pass",
             "1D fail": "Read fail"}

    colors = {"1D": toulligqc_colors['all'],
              "1D pass": toulligqc_colors['pass'],
              "1D fail": toulligqc_colors['fail']}

    # Max yaxis value for displaying same scale between plots
    max_yaxis = (dataframe.max(skipna=True, numeric_only=True).values.max() + 2.0)
    min_yaxis = (dataframe.min(skipna=True, numeric_only=True).values.min() - 2.0)

    fig = go.Figure()

    for column in dataframe.columns:
        d=_precompute_boxplot_values(dataframe[column])
        fig.add_trace(go.Box(
            q1=[d['q1']], median=[d['median']], q3=[d['q3']], lowerfence=[d['lowerfence']], upperfence=[d['upperfence']],
            name=names[column],
            x0=names[column],
            marker=dict(
                opacity=0.3,
                color=colors[column]

            ),
            boxmean=False,
            showlegend=True
        ))

        fig.add_trace(go.Violin(y=dataframe[column],
                            name=names[column],
                            meanline_visible=True,
                      marker=dict(color=colors[column]),
                      visible = False))

    fig.update_layout(
        title={
            'text': "<b>PHRED score distribution of all read types</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        size=title_size,
                        color="black")},
        xaxis=dict(
            title="<b>Read type</b>",
            titlefont_size=axis_font_size
        ),
        yaxis=dict(
            title='<b>PHRED score</b>',
            titlefont_size=axis_font_size,
            tickfont_size=axis_font_size,
            range=[min_yaxis, max_yaxis]
        ),
        legend=dict(
            x=1.02,
            y=0.95,
            title_text="<b>Legend</b>",
            title=dict(font=dict(size=legend_font_size)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=legend_font_size)
        ),
        hovermode='x',
        height=figure_image_height, width=figure_image_width
    )

    # Add buttons
    fig.update_layout(
        updatemenus=[
            dict(
                type = "buttons",
                direction = "left",
                buttons=list([
                    dict(
                        args=[{'visible': [True, False]}],
                        label="Boxplot",
                        method="update"
                    ),
                    dict(
                        args=[{'visible': [False, True]}],
                        label="Violin plot",
                        method="update"
                    )
                ]),
                pad={"r": 20, "t": 20, "l":20, "b":20},
                showactive=True,
                x=1.0,
                xanchor="left",
                y=1.25,
                yanchor="top"
            ),
        ]
    )

    df = df[["1D", "1D pass", "1D fail"]]
    table_html = _dataFrame_to_html(_make_describe_dataframe(df))

    div, output_file = _create_and_save_div(fig, result_directory, main)
    return main, output_file, table_html, desc, div



def allphred_score_frequency(result_dict, main, my_dpi, result_directory, desc):
    """
    Plot the distribution of the phred score per read type (1D , 1D pass, 1D fail)
    """
    dataframe = \
        pd.DataFrame({"1D": result_dict["basecaller.sequencing.summary.1d.extractor.mean.qscore"],
                      "1D pass": result_dict['basecaller.sequencing.summary.1d.extractor.read.pass.qscore'],
                      "1D fail": result_dict['basecaller.sequencing.summary.1d.extractor.read.fail.qscore']})

    # If more than 10.000 reads, interpolate data
    if len(dataframe["1D"]) > interpolation_threshold:
        phred_score_pass = _interpolate(dataframe["1D pass"], npoints=5000)
        phred_score_fail = _interpolate(dataframe["1D fail"], npoints=5000)
    else:
        phred_score_pass = dataframe["1D pass"]
        phred_score_fail = dataframe["1D fail"]

    arr_1D_pass = np.array(pd.Series(phred_score_pass).dropna())
    x = np.linspace(0, max(arr_1D_pass), 200)
    mu, std = norm.fit(arr_1D_pass)
    pdf_1D_pass = norm.pdf(x, mu, std)

    arr_1D_fail = np.array(pd.Series(phred_score_fail).dropna())
    x2 = np.linspace(0, max(arr_1D_fail), 200)
    mu2, std2 = norm.fit(arr_1D_fail)
    pdf_1D_fail = norm.pdf(x2, mu2, std2)

    fig = go.Figure()
    fig.add_trace(go.Histogram(x=phred_score_pass, name="Read pass", marker_color=toulligqc_colors['pass'], histnorm='probability density'))
    fig.add_trace(go.Histogram(x=phred_score_fail, name="Read fail", marker_color=toulligqc_colors['fail'], histnorm='probability density'))
    fig.add_trace(go.Scatter(x=x, y=pdf_1D_pass, mode="lines", name='Density curve of read pass', line=dict(color=toulligqc_colors['pass'], width=3, shape="spline", smoothing=0.5)))
    fig.add_trace(go.Scatter(x=x2, y=pdf_1D_fail, mode="lines", name='Density curve of read fail', line=dict(color=toulligqc_colors['fail'], width=3, shape="spline", smoothing=0.5)))

    fig.update_layout(
        title={
            'text': "<b>PHRED Score Density Distribution</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        size=title_size,
                        color="black")},
        xaxis=dict(
            title="<b>PHRED score</b>",
            titlefont_size=axis_font_size
        ),
        yaxis=dict(
            title='<b>Density probability</b>',
            titlefont_size=axis_font_size,
            tickfont_size=axis_font_size,
        ),
        legend=dict(
            x=1.02,
            y=0.95,
            title_text="<b>Legend</b>",
            title=dict(font=dict(size=legend_font_size)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=legend_font_size)
        ),
        barmode='group',
        hovermode='x',
        height=figure_image_height, width=figure_image_width
    )

    dataframe = _make_describe_dataframe(dataframe).drop('count')
    table_html = _dataFrame_to_html(dataframe)

    div, output_file = _create_and_save_div(fig, result_directory, main)
    return main, output_file, table_html, desc, div


def all_scatterplot(result_dict, main, my_dpi, result_directory, desc):
    """
    Plot the scatter plot representing the relation between the phred score and the sequence length in log
    """
    read_pass_length = result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.length"]
    read_pass_qscore = result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.qscore"]
    read_fail_length = result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.length"]
    read_fail_qscore = result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.qscore"]

    # If more than 10.000 reads, interpolate data
    if len(read_pass_length) > interpolation_threshold:
        pass_data = _interpolate(read_pass_length, 4000, y=read_pass_qscore, interp_type="nearest")
        fail_data = _interpolate(read_fail_length, 4000, y=read_fail_qscore, interp_type="nearest")
    else:
        pass_data = [read_pass_length, read_pass_qscore]
        fail_data = [read_fail_length, read_fail_qscore]
    fig = go.Figure()

    fig.add_trace(go.Scatter(x=pass_data[0],
                             y=pass_data[1],
                             name="Pass reads",
                             marker_color=toulligqc_colors['pass'],
                             mode="markers"
                             ))

    fig.add_trace(go.Scatter(x=fail_data[0],
                             y=fail_data[1],
                             name="Fail reads",
                             marker_color=toulligqc_colors['fail'],
                             mode="markers"
                             ))

    fig.update_layout(
        title={
            'text': "<b>Correlation between read length and PHRED score</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        size=title_size,
                        color="black")},
        xaxis=dict(
            title="<b>Sequence length (bp)</b>",
            titlefont_size=axis_font_size
        ),
        yaxis=dict(
            title='<b>PHRED score</b>',
            titlefont_size=axis_font_size,
            tickfont_size=axis_font_size,
        ),
        legend=dict(
            x=1.02,
            y=.5,
            title_text="<b>Read Type</b>",
            title=dict(font=dict(size=legend_font_size)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=legend_font_size)
        ),
        height=figure_image_height, width=figure_image_width
    )
    # Trim x axis to avoid negative values
    if max(read_pass_length) >= max(read_fail_length):
        max_val = max(read_pass_length)
    max_val = max(read_fail_length)

    fig.update_xaxes(range=[0, max_val])

    table_html = None
    div, output_file = _create_and_save_div(fig, result_directory, main)
    return main, output_file, table_html, desc, div


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
    plt.figure(figsize=(figure_image_width / my_dpi, figure_image_height / my_dpi), dpi=my_dpi)
    sns.heatmap(d, fmt="", linewidths=.5, cmap="YlGnBu", annot_kws={"size": 7},
                cbar_kws={'label': 'Read number per pore channel', "orientation": "horizontal"})

    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

    table_html = None
    return main, output_file, table_html, desc

#
# For each barcode 1D
#


def barcode_percentage_pie_chart_pass(result_dict, dataframe_dict, main, barcode_selection, my_dpi, result_directory, desc):
    """
    Plots a pie chart of 1D read pass percentage per barcode of a run.
    """
    for element in barcode_selection:

        if all(dataframe_dict['barcode.arrangement'] != element):
            print("The barcode {} doesn't exist".format(element))
            return False

    count_sorted = dataframe_dict["read.pass.barcoded"]
    labels = count_sorted.index.values.tolist()

    fig = go.Figure(data=[go.Pie(labels=labels,
                                 values=count_sorted)])
    if len(labels) <= 12:
        palette = ["f3a683", "f7d794", "778beb", "e77f67", "cf6a87", "786fa6", "f8a5c2", "63cdda", "ea8685", "596275", "#b8e994", "#78e08f"]
        fig.update_traces(hoverinfo='label+percent', textinfo='percent', textfont_size=on_chart_font_size,
                  marker=dict(colors=palette, line=dict(width=line_width)))
    else:
        fig.update_traces(hoverinfo='label+percent', textinfo='percent', textfont_size=on_chart_font_size,
                  marker=dict(line=dict(width=line_width)))
    fig.update_traces(textposition='inside')
    fig.update_layout(uniformtext_minsize=12, uniformtext_mode='hide')
    fig.update_layout(
        title={
            'text': "<b>Read Pass Barcode Distribution</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        size=title_size,
                        color="black")},
        legend=dict(
            x=1.02,
            y=.5,
            title_text="<b>Barcodes</b>",
            title=dict(font=dict(size=legend_font_size)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=legend_font_size)
        ),
        height=figure_image_height, width=figure_image_width
    )

    barcode_table = pd.DataFrame({"barcode arrangement": count_sorted/sum(count_sorted)*100,
                                 "read count": count_sorted})
    barcode_table.sort_index(inplace=True)
    pd.options.display.float_format = percent_format_str.format
    barcode_table["read count"] = barcode_table["read count"].astype(int).apply(lambda x: int_format_str.format(x))
    table_html = _dataFrame_to_html(barcode_table)

    div, output_file = _create_and_save_div(fig, result_directory, main)
    return main, output_file, table_html, desc, div


def barcode_percentage_pie_chart_fail(result_dict, dataframe_dict, main, barcode_selection, my_dpi, result_directory, desc):
    """
    Plots a pie chart of 1D read fail percentage per barcode of a run.
    Needs the samplesheet file describing the barcodes to run
    """
    for element in barcode_selection:

        if all(dataframe_dict['barcode.arrangement'] != element):
            print("The barcode {} doesn't exist".format(element))
            return False

    count_sorted = dataframe_dict["read.fail.barcoded"]
    labels = count_sorted.index.values.tolist()

    fig = go.Figure(data=[go.Pie(labels=labels,
                                 values=count_sorted)])
    if len(labels) <= 12:
        palette = ["f3a683", "f7d794", "778beb", "e77f67", "cf6a87", "786fa6", "f8a5c2", "63cdda", "ea8685", "596275"]
        fig.update_traces(hoverinfo='label+percent', textinfo='percent', textfont_size=14,
                  marker=dict(colors=palette, line=dict(width=line_width)))
    else:
        fig.update_traces(hoverinfo='label+percent', textinfo='percent', textfont_size=14,
                  marker=dict(line=dict(width=line_width)))
    fig.update_traces(textposition='inside')
    fig.update_layout(uniformtext_minsize=12, uniformtext_mode='hide')
    fig.update_layout(
        title={
            'text': "<b>Read Fail Barcode Distribution</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        size=title_size,
                        color="black")},
        legend=dict(
            x=1.02,
            y=.5,
            title_text="<b>Barcodes</b>",
            title=dict(font=dict(size=legend_font_size)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=legend_font_size)
        ),
        height=figure_image_height, width=figure_image_width
    )

    barcode_table = pd.DataFrame({"barcode arrangement": count_sorted/sum(count_sorted)*100,
                                  "read count": count_sorted})
    barcode_table.sort_index(inplace=True)
    pd.options.display.float_format = percent_format_str.format
    barcode_table["read count"] = barcode_table["read count"].astype(int).apply(lambda x: int_format_str.format(x))
    table_html = _dataFrame_to_html(barcode_table)

    div, output_file = _create_and_save_div(fig, result_directory, main)
    return main, output_file, table_html, desc, div


def barcode_length_boxplot(result_dict, datafame_dict, main, my_dpi, result_directory, desc):
    """
    Boxplots all the 1D pass and fail read length for each barcode indicated in the sample sheet
    """
    df = datafame_dict['barcode_selection_sequence_length_dataframe']

    # Sort reads by read type and drop read type column
    read_pass_length = df.loc[df['passes_filtering']
                              == bool(True)].drop(columns='passes_filtering')
    read_fail_length = df.loc[df['passes_filtering']
                              == bool(False)].drop(columns='passes_filtering')

    # Remove negative values from all columns of the dataframes
    for col in read_pass_length.columns:
        read_pass_length[col] = read_pass_length[col].loc[read_pass_length[col] > 0]
    for col in read_fail_length.columns:
        read_fail_length[col] = read_fail_length[col].loc[read_fail_length[col] > 0]

    fig = go.Figure()

    for col in read_pass_length.columns:
        d = _precompute_boxplot_values(read_pass_length[col])
        fig.add_trace(go.Box(
            q1=[d['q1']], median=[d['median']], q3=[d['q3']], lowerfence=[d['lowerfence']], upperfence=[d['upperfence']],
            name=col,
            x0=col,
            marker_color=toulligqc_colors['pass'],
            legendgroup="pass",
            offsetgroup="pass"
        ))

    for col in read_fail_length.columns:
        d = _precompute_boxplot_values(read_fail_length[col])
        fig.add_trace(go.Box(
            q1=[d['q1']], median=[d['median']], q3=[d['q3']], lowerfence=[d['lowerfence']], upperfence=[d['upperfence']],
            name=col,
            x0=col,
            marker_color=toulligqc_colors['fail'],
            legendgroup="fail",
            offsetgroup="fail"
        ))

    fig.update_layout(
        title={
            'text': "<b>Read size distribution for each barcode</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        size=title_size,
                        color="black")},
        xaxis=dict(
            title="<b>Barcodes</b>",
            titlefont_size=axis_font_size
        ),
        yaxis=dict(
            title='<b>Sequence length (bp)</b>',
            titlefont_size=axis_font_size,
            tickfont_size=axis_font_size,
        ),
        legend=dict(
            x=1.02,
            y=.5,
            title_text="<b>Read Type</b>",
            title=dict(font=dict(size=legend_font_size)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=legend_font_size)
        ),
        boxmode='group',
        boxgap=0.4,
        boxgroupgap=0,
        height=figure_image_height, width=figure_image_width
    )

    all_read = df.describe().T
    read_pass = df.loc[df['passes_filtering'] == bool(True)].describe().T
    read_fail = df.loc[df['passes_filtering'] == bool(False)].describe().T
    concat = pd.concat([all_read, read_pass, read_fail],
                       keys=['1D', '1D pass', '1D fail'])
    dataframe = concat.T

    dataframe.loc['count'] = dataframe.loc['count'].astype(int).apply(lambda x: int_format_str.format(x))
    dataframe.iloc[1:] = dataframe.iloc[1:].applymap(float_format_str.format)
    table_html = _dataFrame_to_html(dataframe)

    table_html = None
    div, output_file = _create_and_save_div(fig, result_directory, main)
    return main, output_file, table_html, desc, div


def barcoded_phred_score_frequency(barcode_selection, dataframe_dict, main, my_dpi, result_directory, desc):
    """
    Plot boxplot of the 1D pass and fail read qscore for each barcode indicated in the sample sheet
    """
    df = dataframe_dict['barcode_selection_sequence_phred_melted_dataframe']
    barcode_list = barcode_selection

    # Sort reads by read type and drop read type column
    read_pass_qscore = df.loc[df['passes_filtering'] == bool(True)].drop(columns='passes_filtering')
    read_fail_qscore = df.loc[df['passes_filtering'] == bool(False)].drop(columns='passes_filtering')

    fig = go.Figure()

    for barcode in barcode_list:
        final_df = read_pass_qscore.loc[read_pass_qscore['barcodes'] == barcode].dropna()
        d = _precompute_boxplot_values(final_df['qscore'])
        fig.add_trace(go.Box(
                             q1=[d['q1']], median=[d['median']], q3=[d['q3']], lowerfence=[d['lowerfence']], upperfence=[d['upperfence']],
                             name=barcode,
                             x0=barcode,
                             marker_color=toulligqc_colors['pass'],
                             legendgroup="pass",
                             offsetgroup="pass"
                             ))

    for barcode in barcode_list:
        final_df = read_fail_qscore.loc[read_fail_qscore['barcodes'] == barcode].dropna()
        d = _precompute_boxplot_values(final_df['qscore'])
        fig.add_trace(go.Box(
                             q1=[d['q1']], median=[d['median']], q3=[d['q3']], lowerfence=[d['lowerfence']], upperfence=[d['upperfence']],
                             name=barcode,
                             x0=barcode,
                             marker_color=toulligqc_colors['fail'],
                             legendgroup="fail",
                             offsetgroup="fail"
                             ))

    fig.update_layout(
        title={
            'text': "<b>PHRED score distribution for each barcode</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        size=title_size,
                        color="black")},
        xaxis=dict(
            title="<b>Barcodes</b>",
            titlefont_size=axis_font_size
        ),
        yaxis=dict(
            title='<b>PHRED score</b>',
            titlefont_size=axis_font_size,
            tickfont_size=axis_font_size,
        ),
        legend=dict(
            x=1.02,
            y=.5,
            title_text="<b>Read Type</b>",
            title=dict(font=dict(size=legend_font_size)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=legend_font_size)
        ),
        boxmode='group',
        boxgap=0.4,
        boxgroupgap=0,
        height=figure_image_height, width=figure_image_width
    )

    all_read = df.describe().T
    read_pass = df.loc[df['passes_filtering'] == bool(True)].describe().T
    read_fail = df.loc[df['passes_filtering'] == bool(False)].describe().T
    concat = pd.concat([all_read, read_pass, read_fail], keys=['1D', '1D pass', '1D fail'])
    dataframe = concat.T
    dataframe.loc['count'] = dataframe.loc['count'].astype(int).apply(lambda x: int_format_str.format(x))
    dataframe.iloc[1:] = dataframe.iloc[1:].applymap(float_format_str.format)
    table_html = _dataFrame_to_html(dataframe)

    table_html = None
    div, output_file = _create_and_save_div(fig, result_directory, main)
    return main, output_file, table_html, desc, div


def sequence_length_over_time(time_df, dataframe_dict, main, my_dpi, result_directory, desc):

        time = [t/3600 for t in time_df.dropna()]
        time = np.array(sorted(time))

        length = dataframe_dict.get('sequence.length')

         # If more than 10.000 reads, interpolate data
        if len(length) > interpolation_threshold:
            df_time, df_length = _interpolate(time, 200, length, "linear")
        else:
            df_time = time
            df_length = length

        fig = go.Figure()

        fig.add_trace(go.Scatter(
        x=df_time,
        y=df_length,
        fill='tozeroy',
        mode='lines',
        line=dict(color=toulligqc_colors['sequence_length_over_time'],
                  width=line_width,
                  shape="spline")))

        fig.update_layout(
                title={
                'text': "<b>Read length over experiment time</b>",
                'y':0.95,
                'x':0,
                'xanchor': 'left',
                'yanchor': 'top',
                'font' : dict(
                size=title_size,
                color="black")},
            xaxis=dict(
                title="<b>Experiment time (hours)</b>",
                titlefont_size=axis_font_size
                ),
            yaxis=dict(
                title='<b>Read length (bp)</b>',
                titlefont_size=axis_font_size,
                tickfont_size=axis_font_size,
            ),
            hovermode='x',
            height=figure_image_height,
            width=figure_image_width
        )

        table_html = None
        div, output_file = _create_and_save_div(fig, result_directory, main)
        return main, output_file, table_html, desc, div


def phred_score_over_time(qscore_df, time_df, main, my_dpi, result_directory, desc):

        # Time data
        time = [t/3600 for t in time_df.dropna()]
        time = np.array(sorted(time))

        # Qscore data
        qscore = qscore_df.dropna()

        #If more than 10.000 reads, interpolate data
        if len(qscore) > interpolation_threshold:
            df_time, df_qscore = _interpolate(time, 100, qscore, "nearest")
        else:
            df_time = time
            df_qscore = qscore

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=df_time,
            y=df_qscore,
            fill='tozeroy',
            mode='lines',
            line=dict(color=toulligqc_colors['phred_score_over_time'],
                      width=line_width,
                      shape="spline")))

        fig.update_layout(
                title={
                'text': "<b>PHRED score over experiment time</b>",
                'y':0.95,
                'x':0,
                'xanchor': 'left',
                'yanchor': 'top',
                'font' : dict(
                size=title_size,
                color="black")},
            xaxis=dict(
                title="<b>Experiment time (hours)</b>",
                titlefont_size=axis_font_size
                ),
            yaxis=dict(
                title='<b>PHRED quality score</b>',
                titlefont_size=axis_font_size,
                tickfont_size=axis_font_size,
            ),
            height=figure_image_height, width=figure_image_width
        )

        table_html = None
        div, output_file = _create_and_save_div(fig, result_directory, main)
        return main, output_file, table_html, desc, div


def speed_over_time(duration_df, sequence_length_df, time_df, main, my_dpi, result_directory, desc):

        speed = pd.Series(sequence_length_df / duration_df)

        time = [t/3600 for t in time_df]
        time = np.array(sorted(time))

        # If more than 10.000 reads, interpolate data
        if len(time) > interpolation_threshold:
            time_df, speed_df = _interpolate(time, 200, speed, "linear")
        else:
            time_df = time
            speed_df = speed

        fig = go.Figure()

        fig.add_trace(go.Scatter(
        x=time_df,
        y=speed_df,
        fill='tozeroy',
        mode='lines',
        line=dict(color=toulligqc_colors['speed_over_time'],
                  width=line_width,
                  shape="spline")))

        fig.update_layout(
                title={
                'text': "<b>Speed over experiment time</b>",
                'y':0.95,
                'x':0,
                'xanchor': 'left',
                'yanchor': 'top',
                'font' : dict(
                size=title_size,
                color="black")},
            xaxis=dict(
                title="<b>Experiment time (hours)</b>",
                titlefont_size=axis_font_size
                ),
            yaxis=dict(
                title='<b>Speed (bases per second)</b>',
                titlefont_size=axis_font_size,
                tickfont_size=axis_font_size,
            ),
            hovermode='x',
            height=figure_image_height,
            width=figure_image_width
        )
        fig.update_yaxes(type="log")

        table_html = None
        div, output_file = _create_and_save_div(fig, result_directory, main)
        return main, output_file, table_html, desc, div


def nseq_over_time(time_df, main, my_dpi, result_directory, desc):

        time = [t/3600 for t in time_df]
        time = pd.Series(time)

        # create custom xaxis points to reduce graph size
        time_points = np.linspace(min(time), max(time), 50)
        n_seq = time.groupby(pd.cut(time, time_points, right=True)).count()

        fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=time_points,
            y=list(n_seq.values),
            mode='lines',
            fill="tozeroy",
            line=dict(color=toulligqc_colors['nseq_over_time'],
                      width=line_width,
                      shape="spline")))

        fig.update_layout(
                title={
                'text': "<b>Number of sequences through experiment time</b>",
                'y':0.95,
                'x':0,
                'xanchor': 'left',
                'yanchor': 'top',
                'font' : dict(
                size=title_size,
                color="black")},
            xaxis=dict(
                title="<b>Experiment time (hours)</b>",
                titlefont_size=axis_font_size
                ),
            yaxis=dict(
                title='<b>Number of sequences</b>',
                titlefont_size=axis_font_size,
                tickfont_size=axis_font_size,
            ),
            hovermode='x',
            height=figure_image_height,
            width=figure_image_width
        )

        table_html = None
        div, output_file = _create_and_save_div(fig, result_directory, main)
        return main, output_file, table_html, desc, div
