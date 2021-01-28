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
from scipy.stats import norm

from toulligqc.plotly_graph_common import _create_and_save_div
from toulligqc.plotly_graph_common import _dataFrame_to_html
from toulligqc.plotly_graph_common import _interpolate
from toulligqc.plotly_graph_common import _make_describe_dataframe
from toulligqc.plotly_graph_common import _precompute_boxplot_values
# Import common function
from toulligqc.plotly_graph_common import _smooth_data
from toulligqc.plotly_graph_common import _transparent_colors
from toulligqc.plotly_graph_common import axis_font_size
from toulligqc.plotly_graph_common import figure_image_height
# Import common constants
from toulligqc.plotly_graph_common import figure_image_width
from toulligqc.plotly_graph_common import float_format_str
from toulligqc.plotly_graph_common import graph_font
from toulligqc.plotly_graph_common import int_format_str
from toulligqc.plotly_graph_common import interpolation_threshold
from toulligqc.plotly_graph_common import legend_font_size
from toulligqc.plotly_graph_common import line_width
from toulligqc.plotly_graph_common import on_chart_font_size
from toulligqc.plotly_graph_common import plotly_background_color
from toulligqc.plotly_graph_common import title_size
from toulligqc.plotly_graph_common import toulligqc_colors


#
#  1D² plots
#

def dsqr_read_count_histogram(result_dict, dataframe_dict_1dsqr, result_directory):
    """
    Plots the histogram of 1D² count of the different types of reads:
    1D² read return by Guppy
    1D² pass read return by Guppy (Qscore >= 7)
    1D² fail read return by Guppy (Qscore < 7)
    """

    graph_name = "1D² Read count histogram"

    # Histogram with barcoded read counts
    if 'read.pass.barcoded.count' in dataframe_dict_1dsqr:

        data = {
            '1D read count': result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
            '1D² read count': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.count"],
            '1D² read pass count': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.count"],
            '1D² read fail count': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.count"],
            '1D² read pass barcoded count': dataframe_dict_1dsqr["read.pass.barcoded.count"],
            '1D² read fail barcoded count': dataframe_dict_1dsqr["read.fail.barcoded.count"]
        }

        colors = [toulligqc_colors["all"], toulligqc_colors["all_1d2"], toulligqc_colors["pass"],
                  toulligqc_colors["fail"], toulligqc_colors["barcode_pass"], toulligqc_colors["barcode_fail"]]

        trace = go.Bar(x=[*data], y=list(data.values()),
                       hovertext=["<b>Total number of reads</b>",
                                  "<b>Number of 1D² reads</b>",
                                  "<b>1D² pass reads</b>",
                                  "<b>1D² fail reads</b>",
                                  "<b>1D² Barcoded pass reads</b>",
                                  "<b>1D² Barcoded fail reads</b>"],
                       name="Barcoded graph",
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
              dataframe_dict_1dsqr["read.pass.barcoded.count"],
              dataframe_dict_1dsqr["read.fail.barcoded.count"]],
             # frequencies
             [result_dict["basecaller.sequencing.summary.1d.extractor.read.count.frequency"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.count.frequency"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.frequency"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.frequency"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.barcoded.frequency"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.barcoded.frequency"]]])

        dataframe = pd.DataFrame(array, index=['count', 'frequency'],
                                 columns=["Read count", "1D² Read count", "1D² pass", "1D² fail", "1D² pass barcoded",
                                          "1D² fail barcoded"])

    # Histogram without barcodes
    else:

        data = {
            'Read Count': result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
            '1D² Read Count': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.count"],
            'Read Pass Count': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.count"],
            'Read Fail Count': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.count"]
        }

        colors = [toulligqc_colors["all"], toulligqc_colors["all_1d2"], toulligqc_colors["pass"],
                  toulligqc_colors["fail"]]

        trace = go.Bar(x=[*data], y=list(data.values()),
                       hovertext=["<b>Total number of reads</b>",
                                  "<b>Number of 1D² reads</b>",
                                  "<b>1D² pass reads</b>",
                                  "<b>1D² fail reads</b>"],
                       name="Barcoded graph",
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
        dataframe = pd.DataFrame(array, index=['count', 'frequency'],
                                 columns=["Read count", "1D² Read count", "1D² pass", "1D² fail"])

    layout = go.Layout(
        hovermode="x",
        title={
            'text': "<b>" + graph_name + "</b>",
            'y': 0.95,
            'x': 0,
            'xanchor': 'left',
            'yanchor': 'top',
            'font': dict(
                size=title_size,
                color="black")
        },
        xaxis=dict(title="<b>1D² Read type</b>",
                   titlefont=dict(
                       size=axis_font_size,
                       color="black"
                   ),
                   categoryorder="trace"
                   ),
        yaxis=dict(title="<b>Counts</b>",
                   titlefont=dict(
                       size=axis_font_size,
                       color="black"
                   )),
        font=dict(family=graph_font),
        width=figure_image_width,
        height=figure_image_height)

    fig = go.Figure(data=trace, layout=layout)

    # HTML table
    dataframe.iloc[0] = dataframe.iloc[0].astype(int).apply(lambda x: int_format_str.format(x))
    dataframe.iloc[1:] = dataframe.iloc[1:].applymap(float_format_str.format)
    table_html = _dataFrame_to_html(dataframe)
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def dsqr_read_length_scatterplot(result_dict, sequence_length_1dsqr, result_directory):
    graph_name = "1D² Distribution of read lengths"

    all_read = sequence_length_1dsqr.loc[sequence_length_1dsqr >= 10].dropna().values
    read_pass = result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.length'].loc[
        result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.length'] >= 10]
    read_fail = result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.fail.length'].loc[
        result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.fail.length'] >= 10]

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
            'text': "<b>" + graph_name + "</b>",
            'y': 0.95,
            'x': 0,
            'xanchor': 'left',
            'yanchor': 'top',
            'font': dict(
                size=title_size,
                color="black")},
        xaxis=dict(
            title="<b>1D² Read length (bp)</b>",
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
            font=dict(size=15)
        ),
        hovermode='x',
        font=dict(family=graph_font),
        height=figure_image_height,
        width=figure_image_width
    )

    table_html = None
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def dsqr_read_quality_multiboxplot(result_dict, dataframe_dict_1dsqr, result_directory):
    """
    Boxplot of PHRED score between read pass and read fail
    Violin plot of PHRED score between read pass and read fail
    """

    graph_name = "1D² PHRED score distribution of all read types"

    df = pd.DataFrame(
        {"1D²": dataframe_dict_1dsqr["mean.qscore"],
         "1D² pass": result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.qscore'],
         "1D² fail": result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.fail.qscore']
         })

    # If more than 10.000 reads, interpolate data
    if len(df["1D²"]) > interpolation_threshold:
        dataframe = pd.DataFrame({
            "1D²": _interpolate(df["1D²"], 1000),
            "1D² pass": _interpolate(df["1D² pass"], 1000),
            "1D² fail": _interpolate(df["1D² fail"], 1000)
        })
    else:
        dataframe = df
    names = {"1D²": "All reads",
             "1D² pass": "Read pass",
             "1D² fail": "Read fail"}

    colors = {"1D²": toulligqc_colors['all'],
              "1D² pass": toulligqc_colors['pass'],
              "1D² fail": toulligqc_colors['fail']}

    # Max yaxis value for displaying same scale between plots
    max_yaxis = (dataframe.max(skipna=True, numeric_only=True).values.max() + 2.0)
    min_yaxis = (dataframe.min(skipna=True, numeric_only=True).values.min() - 2.0)

    fig = go.Figure()

    for column in dataframe.columns:
        d = _precompute_boxplot_values(dataframe[column])
        fig.add_trace(go.Box(
            q1=[d['q1']], median=[d['median']], q3=[d['q3']], lowerfence=[d['lowerfence']],
            upperfence=[d['upperfence']],
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
                                visible=False))

    fig.update_layout(
        title={
            'text': "<b>" + graph_name + "</b>",
            'y': 0.95,
            'x': 0,
            'xanchor': 'left',
            'yanchor': 'top',
            'font': dict(
                size=title_size,
                color="black")},
        xaxis=dict(
            title="<b>1D² Read type</b>",
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
            font=dict(size=15)
        ),
        hovermode='x',
        font=dict(family=graph_font),
        height=figure_image_height,
        width=figure_image_width
    )

    # Add buttons
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="left",
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
                pad={"r": 20, "t": 20, "l": 20, "b": 20},
                showactive=True,
                x=1.0,
                xanchor="left",
                y=1.25,
                yanchor="top"
            ),
        ]
    )

    df = df[["1D²", "1D² pass", "1D² fail"]]
    table_html = _dataFrame_to_html(_make_describe_dataframe(df))
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def dsqr_allphred_score_frequency(result_dict, dataframe_dict_1dsqr, result_directory):
    """
    Plot the distribution of the phred score per read type (1D² , 1D² pass, 1D² fail)
    """

    graph_name = "1D² PHRED score density distribution"

    dataframe = \
        pd.DataFrame({"1D²": dataframe_dict_1dsqr["mean.qscore"],
                      "1D² pass": result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.qscore'],
                      "1D² fail": result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.fail.qscore']})

    # If more than 10.000 reads, interpolate data
    if len(dataframe["1D²"]) > interpolation_threshold:
        phred_score_pass = _interpolate(dataframe["1D² pass"], npoints=5000)
        phred_score_fail = _interpolate(dataframe["1D² fail"], npoints=5000)
    else:
        phred_score_pass = dataframe["1D² pass"]
        phred_score_fail = dataframe["1D² fail"]

    arr_1Dsqr_pass = np.array(pd.Series(phred_score_pass).dropna())
    x = np.linspace(0, max(arr_1Dsqr_pass), 200)
    mu, std = norm.fit(arr_1Dsqr_pass)
    pdf_1Dsqr_pass = norm.pdf(x, mu, std)

    arr_1Dsqr_fail = np.array(pd.Series(phred_score_fail).dropna())
    x2 = np.linspace(0, max(arr_1Dsqr_fail), 200)
    mu2, std2 = norm.fit(arr_1Dsqr_fail)
    pdf_1Dsqr_fail = norm.pdf(x2, mu2, std2)

    fig = go.Figure()
    fig.add_trace(go.Histogram(x=phred_score_pass, name="Read pass", marker_color=toulligqc_colors['pass'],
                               histnorm='probability density'))
    fig.add_trace(go.Histogram(x=phred_score_fail, name="Read fail", marker_color=toulligqc_colors['fail'],
                               histnorm='probability density'))
    fig.add_trace(go.Scatter(x=x, y=pdf_1Dsqr_pass, mode="lines", name='Density curve of read pass',
                             line=dict(color=toulligqc_colors['pass'], width=3, shape="spline", smoothing=0.5)))
    fig.add_trace(go.Scatter(x=x2, y=pdf_1Dsqr_fail, mode="lines", name='Density curve of read fail',
                             line=dict(color=toulligqc_colors['fail'], width=3, shape="spline", smoothing=0.5)))

    fig.update_layout(
        title={
            'text': "<b>" + graph_name + "</b>",
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
            font=dict(size=15)
        ),
        barmode='group',
        hovermode='x',
        font=dict(family=graph_font),
        height=figure_image_height,
        width=figure_image_width
    )

    dataframe = _make_describe_dataframe(dataframe).drop('count')

    table_html = _dataFrame_to_html(dataframe)
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def scatterplot_1dsqr(result_dict, result_directory):
    """
    Plot the scatter plot representing the relation between the phred score and the sequence length in log
    """

    graph_name = "Correlation between 1D² read length and PHRED score"

    read_pass_length = result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.length"]
    read_pass_qscore = result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.qscore"]
    read_fail_length = result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.length"]
    read_fail_qscore = result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.qscore"]

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
                             name='Fail reads',
                             marker_color=toulligqc_colors['fail'],
                             mode="markers"
                             ))

    fig.update_layout(
        title={
            'text': "<b>" + graph_name + "</b>",
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
            title_text="<b>1D² Read type</b>",
            title=dict(font=dict(size=legend_font_size)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=15)
        ),
        font=dict(family=graph_font),
        height=figure_image_height,
        width=figure_image_width
    )
    # Trim x axis to avoid negative values
    if max(read_pass_length) >= max(read_fail_length):
        max_val = max(read_pass_length)
    max_val = max(read_fail_length)

    fig.update_xaxes(range=[0, max_val])

    table_html = None
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


#
# For each barcode 1D²
#

def barcode_percentage_pie_chart_1dsqr_pass(dataframe_dict_1dsqr, barcode_selection, result_directory):
    """
    Plots a pie chart of 1D² read pass percentage per barcode of a run.
    """

    graph_name = "1D² read pass barcode distribution"

    for element in barcode_selection:

        if all(dataframe_dict_1dsqr['barcode.arrangement'] != element):
            print("The barcode {} doesn't exist".format(element))
            return False

    count_sorted = dataframe_dict_1dsqr["read.pass.barcoded"]
    labels = count_sorted.index.values.tolist()

    fig = go.Figure(data=[go.Pie(labels=labels,
                                 values=count_sorted)])
    if len(labels) <= 12:
        palette = ["f3a683", "f7d794", "778beb", "e77f67", "cf6a87", "786fa6", "f8a5c2", "63cdda", "ea8685", "596275",
                   "#b8e994", "#78e08f"]
        fig.update_traces(hoverinfo='label+percent', textinfo='percent', textfont_size=on_chart_font_size,
                          marker=dict(colors=palette, line=dict(width=line_width, color='#808080')))
    else:
        fig.update_traces(hoverinfo='label+percent', textinfo='percent', textfont_size=on_chart_font_size,
                          marker=dict(line=dict(width=line_width, color='#808080')))
    fig.update_traces(textposition='inside')
    fig.update_layout(uniformtext_minsize=12, uniformtext_mode='hide')
    fig.update_layout(
        title={
            'text': "<b>" + graph_name + "</b>",
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
            font=dict(size=15)
        ),
        font=dict(family=graph_font),
        height=figure_image_height,
        width=figure_image_width
    )

    barcode_table = pd.DataFrame({"Barcode arrangement (%)": count_sorted / sum(count_sorted) * 100,
                                  "1D² read count": count_sorted})
    barcode_table.sort_index(inplace=True)
    pd.options.display.float_format = float_format_str.format
    barcode_table["1D² read count"] = barcode_table["1D² read count"].astype(int).apply(
        lambda x: int_format_str.format(x))
    table_html = _dataFrame_to_html(barcode_table)
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def barcode_percentage_pie_chart_1dsqr_fail(dataframe_dict_1dsqr, barcode_selection, result_directory):
    """
    Plots a pie chart of 1D² read fail percentage per barcode of a run.
    Needs the samplesheet file describing the barcodes to run
    """

    graph_name = "1D² read pass barcode distribution"

    for element in barcode_selection:

        if all(dataframe_dict_1dsqr['barcode.arrangement'] != element):
            print("The barcode {} doesn't exist".format(element))
            return False

    count_sorted = dataframe_dict_1dsqr["read.fail.barcoded"]
    labels = count_sorted.index.values.tolist()

    fig = go.Figure(data=[go.Pie(labels=labels,
                                 values=count_sorted)])
    if len(labels) <= 12:
        palette = ["f3a683", "f7d794", "778beb", "e77f67", "cf6a87", "786fa6", "f8a5c2", "63cdda", "ea8685", "596275"]
        fig.update_traces(hoverinfo='label+percent', textinfo='percent', textfont_size=on_chart_font_size,
                          marker=dict(colors=palette, line=dict(width=line_width, color='#808080')))
    else:
        fig.update_traces(hoverinfo='label+percent', textinfo='percent', textfont_size=on_chart_font_size,
                          marker=dict(line=dict(width=line_width, color='#808080')))
    fig.update_traces(textposition='inside')
    fig.update_layout(uniformtext_minsize=12, uniformtext_mode='hide')
    fig.update_layout(
        title={
            'text': "<b>" + graph_name + "</b>",
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
            font=dict(size=15)
        ),
        font=dict(family=graph_font),
        height=figure_image_height,
        width=figure_image_width
    )

    barcode_table = pd.DataFrame({"Barcode arrangement (%)": count_sorted / sum(count_sorted) * 100,
                                  "1D² read count": count_sorted})
    barcode_table.sort_index(inplace=True)
    pd.options.display.float_format = float_format_str.format
    barcode_table["1D² read count"] = barcode_table["1D² read count"].astype(int).apply(
        lambda x: int_format_str.format(x))
    table_html = _dataFrame_to_html(barcode_table)
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def barcode_length_boxplot_1dsqr(dataframe_dict_1dsqr, result_directory):
    """
    Boxplots all the 1D² pass and fail read length for each barcode indicated in the sample sheet
    """

    graph_name = "1D² read size distribution for each barcode"

    df = dataframe_dict_1dsqr['barcode_selection_sequence_length_dataframe']

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

    first = True
    for col in read_pass_length.columns:
        d = _precompute_boxplot_values(read_pass_length[col])
        fig.add_trace(go.Box(
            q1=[d['q1']], median=[d['median']], q3=[d['q3']], lowerfence=[d['lowerfence']],
            upperfence=[d['upperfence']],
            name="Read pass",
            x0=col,
            marker_color=toulligqc_colors['pass'],
            offsetgroup="pass",
            showlegend=first
        ))
        if first:
            first=False

    first = True
    for col in read_fail_length.columns:
        d = _precompute_boxplot_values(read_fail_length[col])
        fig.add_trace(go.Box(
            q1=[d['q1']], median=[d['median']], q3=[d['q3']], lowerfence=[d['lowerfence']],
            upperfence=[d['upperfence']],
            name="Read fail",
            x0=col,
            marker_color=toulligqc_colors['fail'],
            offsetgroup="fail",
            showlegend=first
        ))
        if first:
            first = False

    fig.update_layout(
        title={
            'text': "<b>" + graph_name + "</b>",
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
            title_text="<b>1D² read type</b>",
            title=dict(font=dict(size=legend_font_size)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=15)
        ),
        boxmode='group',
        boxgap=0.4,
        boxgroupgap=0,
        font=dict(family=graph_font),
        height=figure_image_height,
        width=figure_image_width
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
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def barcoded_phred_score_frequency_1dsqr(barcode_selection, dataframe_dict_1dsqr, result_directory):
    """
    Plot boxplot of the 1D pass and fail read qscore for each barcode indicated in the sample sheet
    """

    graph_name = "1D² PHRED score distribution for each barcode"

    df = dataframe_dict_1dsqr['barcode_selection_sequence_phred_melted_dataframe']
    barcode_list = barcode_selection

    # Sort reads by read type and drop read type column
    read_pass_qscore = df.loc[df['passes_filtering'] == bool(True)].drop(columns='passes_filtering')
    read_fail_qscore = df.loc[df['passes_filtering'] == bool(False)].drop(columns='passes_filtering')

    fig = go.Figure()

    first = True
    for barcode in barcode_list:
        final_df = read_pass_qscore.loc[read_pass_qscore['barcodes'] == barcode].dropna()
        d = _precompute_boxplot_values(final_df['qscore'])
        fig.add_trace(go.Box(
            q1=[d['q1']], median=[d['median']], q3=[d['q3']], lowerfence=[d['lowerfence']],
            upperfence=[d['upperfence']],
            name="Read pass",
            x0=barcode,
            marker_color=toulligqc_colors['pass'],
            offsetgroup="pass",
            showlegend=first
        ))
        if first:
            first=False

    first = True
    for barcode in barcode_list:
        final_df = read_fail_qscore.loc[read_fail_qscore['barcodes'] == barcode].dropna()
        d = _precompute_boxplot_values(final_df['qscore'])
        fig.add_trace(go.Box(
            q1=[d['q1']], median=[d['median']], q3=[d['q3']], lowerfence=[d['lowerfence']],
            upperfence=[d['upperfence']],
            name="Read fail",
            x0=barcode,
            marker_color=toulligqc_colors['fail'],
            offsetgroup="fail",
            showlegend=first
        ))
        if first:
            first=False

    fig.update_layout(
        title={
            'text': "<b>" + graph_name + "</b>",
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
            title_text="<b>1D² read type</b>",
            title=dict(font=dict(size=legend_font_size)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=15)
        ),
        boxmode='group',
        boxgap=0.4,
        boxgroupgap=0,
        font=dict(family=graph_font),
        height=figure_image_height,
        width=figure_image_width
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
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def sequence_length_over_time_dsqr(time_df, sequence_length_df, result_directory):
    graph_name = "1D² Read length over experiment time"

    time = [t / 3600 for t in time_df.dropna()]
    time = np.array(sorted(time))
    length = sequence_length_df

    # If more than 10.000 reads, interpolate data
    if len(length) > interpolation_threshold:
        df_time, df_length = _interpolate(time, 200, length, "nearest")
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
            'text': "<b>" + graph_name + "</b>",
            'y': 0.95,
            'x': 0,
            'xanchor': 'left',
            'yanchor': 'top',
            'font': dict(
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
        font=dict(family=graph_font),
        height=figure_image_height,
        width=figure_image_width
    )

    table_html = None
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def phred_score_over_time_dsqr(qscore_df, time_df, result_directory):
    graph_name = "1D² PHRED score over experiment time"

    # Time data
    time = [t / 3600 for t in time_df.dropna()]
    time = np.array(sorted(time))

    # Qscore data
    qscore = qscore_df.dropna()

    # If more than 10.000 reads, interpolate data
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
            'text': "<b>" + graph_name + "</b>",
            'y': 0.95,
            'x': 0,
            'xanchor': 'left',
            'yanchor': 'top',
            'font': dict(
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
        font=dict(family=graph_font),
        height=figure_image_height,
        width=figure_image_width
    )

    table_html = None
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def speed_over_time_dsqr(duration_df, sequence_length_df, time_df, result_directory):
    graph_name = "1D² Speed over experiment time"

    speed = pd.Series(sequence_length_df / duration_df)

    time = [t / 3600 for t in time_df]
    time = np.array(sorted(time))

    # If more than 10.000 reads, interpolate data
    if len(time) > interpolation_threshold:
        time_df, speed_df = _interpolate(time, 200, speed, "nearest")
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
                  shape="spline"))
    )

    fig.update_layout(
        title={
            'text': "<b>" + graph_name + "</b>",
            'y': 0.95,
            'x': 0,
            'xanchor': 'left',
            'yanchor': 'top',
            'font': dict(
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
        font=dict(family=graph_font),
        height=figure_image_height,
        width=figure_image_width
    )
    fig.update_yaxes(type="log")

    table_html = None
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def nseq_over_time_dsqr(time_df, result_directory):
    graph_name = "Number of sequences through experiment time"

    time = [t / 3600 for t in time_df]
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
                  shape="spline")
    ))

    fig.update_layout(
        title={
            'text': "<b>" + graph_name + "</b>",
            'y': 0.95,
            'x': 0,
            'xanchor': 'left',
            'yanchor': 'top',
            'font': dict(
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
        font=dict(family=graph_font),
        height=figure_image_height,
        width=figure_image_width
    )

    table_html = None
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div
