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

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib.pyplot import table

from plotly.subplots import make_subplots
from scipy.interpolate import interp1d
from scipy.stats import norm
from sklearn.utils import resample
import plotly.graph_objs as go
import plotly.offline as py
import plotly.colors as colors
from collections import defaultdict
from scipy.ndimage.filters import gaussian_filter1d

def _make_desribe_dataframe(value):
    """
    Creation of a statistics table printed with the graph in report.html
    :param value: information measured (series)
    """

    desc = value.describe()
    desc.loc['count'] = desc.loc['count'].astype(int).astype(str)
    desc.iloc[1:] = desc.iloc[1:].applymap(lambda x: '%.2f' % x)
    desc.rename({'50%': 'median'}, axis='index', inplace=True)

    return desc

#
#  1D² plots
#


def dsqr_read_count_histogram(result_dict, dataframe_dict_1dsqr, main, my_dpi, result_directory, desc):
    """
    Plots the histogram of 1D² count of the different types of reads:
    1D² read return by Guppy
    1D² pass read return by Guppy (Qscore >= 7)
    1D² fail read return by Guppy (Qscore < 7)
    """
    output_file = result_directory + '/' + '_'.join(main.split())

    # Histogram with barcoded read counts
    if 'read.pass.barcoded.count' in dataframe_dict_1dsqr:

        data = {
            '1D² Read Count': result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
            '1D² Read Count': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.count"],
            '1D² Read Pass Count': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.count"],
            '1D² Read Pass Barcoded Count': dataframe_dict_1dsqr["read.pass.barcoded.count"],
            '1D² Read Fail Count': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.count"],
            '1D² Read Fail Barcoded Count': dataframe_dict_1dsqr["read.fail.barcoded.count"]
        }

        colors = ["#54A8FC", "#cb4335 ", "salmon", '#ffa931', "#50c878", "SlateBlue"]

        trace = go.Bar(x=[*data], y=list(data.values()),
                                hovertext=["<b>Total number of reads</b>",
                                           "<b>Number of 1D² reads</b>",
                                           "<b>1D² Reads of qscore > 7</b>",
                                           "<b>1D² Barcoded reads with qscore > 7</b>",
                                           "<b>1D² Reads of qscore < 7</b>",
                                           "<b>1D² Barcoded reads with qscore < 7</b>"],
                                name="Barcoded graph",
                                marker_color=colors,
                                marker_line_color="black",
                                marker_line_width=1.5, opacity=0.9)

        # Array of data for HTML table with barcode reads
        array = np.array(
            #count
            [[result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.count"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.count"],
              dataframe_dict_1dsqr["read.pass.barcoded.count"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.count"],
              dataframe_dict_1dsqr["read.fail.barcoded.count"]],
             #frequencies
             [result_dict["basecaller.sequencing.summary.1d.extractor.read.count.frequency"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.count.frequency"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.frequency"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.barcoded.frequency"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.frequency"],
              result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.barcoded.frequency"]]])

        dataframe = pd.DataFrame(array, index=['count', 'frequency'],
                                 columns=["Read count", "1D² Read count","1D² pass", "1D² pass barcoded", "1D² fail", "1D² fail barcoded"])

    # Histogram without barcodes
    else:

        data = {
            'Read Count': result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
            '1D² Read Count': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.count"],
            'Read Pass Count': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.count"],
            'Read Fail Count': result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.count"]
        }

        colors = ["#54A8FC", "#cb4335", "salmon", "#50c878"]

        trace = go.Bar(x=[*data], y=list(data.values()),
                       hovertext=["<b>Total number of reads</b>",
                                  "<b>Reads of qscore > 7</b>",
                                  "<b>Barcoded reads with qscore > 7</b>",
                                  "<b>Reads of qscore < 7</b>",
                                  "<b>Barcoded reads with qscore < 7</b>"],
                       name="Barcoded graph",
                       marker_color=colors,
                       marker_line_color="black",
                       marker_line_width=1.5, opacity=0.9)

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
            'text': "<b>1D² Read count histogram</b>",
            'y': 0.95,
            'x': 0,
            'xanchor': 'left',
            'yanchor': 'top',
            'font': dict(
                family="Calibri, sans",
                size=25,
                color="black")
        },
        xaxis=dict(title="<b>1D² Read type</b>",
                   linecolor="black",
                   titlefont=dict(
                       family="Calibri, sans",
                       size=18,
                       color="black"
                   ),
                   categoryorder="trace"
                   ),
        yaxis=dict(title="<b>Counts</b>",
                   linecolor="black",
                   titlefont=dict(
                       family="Calibri, sans",
                       size=18,
                       color="black"
                   )),
        width=800,
        height=600,
        plot_bgcolor="white",
        paper_bgcolor="white"
    )

    fig = go.Figure(data=trace, layout=layout)
    div = py.plot(fig,
                  include_plotlyjs=False,
                  output_type='div',
                  auto_open=False,
                  show_link=False)
    py.plot(fig, filename=output_file, output_type="file", include_plotlyjs="directory", auto_open=False)

    # HTML table
    dataframe.iloc[0] = dataframe.iloc[0].astype(int).astype(str)
    dataframe.iloc[1:] = dataframe.iloc[1:].applymap('{:.2f}'.format)
    table_html = pd.DataFrame.to_html(dataframe)

    return main, output_file, table_html, desc, div

def dsqr_read_length_multihistogram(result_dict, sequence_length_1dsqr, main, my_dpi, result_directory, desc):

    output_file = result_directory + '/' + '_'.join(main.split())
    
    all_read = sequence_length_1dsqr.loc[sequence_length_1dsqr >= 10].values
    read_pass = result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.length'].loc[result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.length'] >= 10]
    read_fail = result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.fail.length'].loc[result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.fail.length'] >= 10]
  
    fig = go.Figure()
    
    fig.add_trace(go.Histogram(x=all_read,
                               name='All reads',
                               nbinsx=500,
                               marker_color='#fca311'  # yellow
                               ))

    fig.add_trace(go.Histogram(x=read_pass,
                               name='Pass reads',
                            nbinsx=500,
                               marker_color='#51a96d'  # green
                               ))

    fig.add_trace(go.Histogram(x=read_fail,
                            nbinsx=500,
                               name='Fail reads',
                               marker_color='#d90429'  # red
                               ))
    
    fig.update_layout(
        title={
            'text': "<b>Distribution of 1D² read length</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        family="Calibri, sans",
                        size=26,
                        color="black")},
        xaxis=dict(
            title="<b>Read length (bp)</b>",
            titlefont_size=16,
            range=[0, 5000]
        ),
        yaxis=dict(
            title='<b>Number of sequences</b>',
            titlefont_size=16,
            tickfont_size=14,
        ),
        legend=dict(
            x=1.02,
            y=0.95,
            title_text="<b>Legend</b>",
            title=dict(font=dict(size=16)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=15)
        ),
        hovermode='x',
        height=800, width=1400
    )

    div = py.plot(fig,
                  include_plotlyjs=False,
                  output_type='div',
                  auto_open=False,
                  show_link=False)
    py.plot(fig, filename=output_file, output_type="file", include_plotlyjs="directory", auto_open=False)

    table_html = None

    return main, output_file, table_html, desc, div

def dsqr_read_quality_multiboxplot(result_dict, dataframe_dict_1dsqr, main, my_dpi, result_directory, desc):
    """
    Boxplot of PHRED score between read pass and read fail
    Violin plot of PHRED score between read pass and read fail
    """
    output_file = result_directory + '/' + '_'.join(main.split())

    df = pd.DataFrame(
        {"1D²": dataframe_dict_1dsqr["mean.qscore"],
         "1D² pass": result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.qscore'],
         "1D² fail": result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.fail.qscore']
         })

    # If more than 10.000 reads, interpolate data
    if len(df["1D²"]) > 10000:
        dataframe = pd.DataFrame({
        "1D²" : _interpolate(df["1D²"], 1000),
        "1D² pass" : _interpolate(df["1D² pass"], 1000),
        "1D² fail" : _interpolate(df["1D² fail"], 1000)
    })
    else:
        dataframe = df
    names = {"1D²": "All reads",
             "1D² pass": "Read pass",
             "1D² fail": "Read fail"}

    colors = {"1D²": '#fca311',
              "1D² pass": '#51a96d',
              "1D² fail": '#d90429'}

    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=("<b>PHRED score boxplot</b>",
                                        "<b>PHRED score violin plot</b>"),
                        horizontal_spacing=0.15)

    for column in dataframe.columns:
        fig.append_trace(go.Box(
            y=dataframe[column],
            name=names[column],
            marker=dict(
                opacity=0.3,
                color=colors[column]

            ),
            boxmean=True,
            showlegend=False
        ), row=1, col=1)

        fig.add_trace(go.Violin(y=dataframe[column],
                            name=names[column],
                            meanline_visible=True,
                      marker=dict(color=colors[column])),
                      row=1, col=2)


    fig.update_layout(
        title={
            'text': "<b>1D² PHRED score distribution of all read types</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        family="Calibri, sans",
                        size=26,
                        color="black")},
        xaxis=dict(
            title="<b>1D² </b>",
            titlefont_size=16
        ),
        yaxis=dict(
            title='<b>PHRED score</b>',
            titlefont_size=16,
            tickfont_size=14,
        ),
        legend=dict(
            x=1.02,
            y=0.95,
            title_text="<b>Legend</b>",
            title=dict(font=dict(size=16)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=15)
        ),
        hovermode='x',
        height=800, width=1400
    )

    div = py.plot(fig,
                  include_plotlyjs=False,
                  output_type='div',
                  auto_open=False,
                  show_link=False)
    py.plot(fig, filename=output_file, output_type="file", include_plotlyjs="directory", auto_open=False)

    df = df[["1D²", "1D² pass", "1D² fail"]]
    table_html = pd.DataFrame.to_html(_make_desribe_dataframe(df))

    return main, output_file, table_html, desc, div

def dsqr_allphred_score_frequency(result_dict, dataframe_dict_1dsqr, main, my_dpi, result_directory, desc):
    """
    Plot the distribution of the phred score per read type (1D² , 1D² pass, 1D² fail)
    """
    output_file = result_directory + '/' + '_'.join(main.split())

    dataframe = \
        pd.DataFrame({"1D²": dataframe_dict_1dsqr["mean.qscore"],
                      "1D² pass": result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.pass.qscore'],
                      "1D² fail": result_dict['basecaller.sequencing.summary.1dsqr.extractor.read.fail.qscore']})

    # If more than 10.000 reads, interpolate data
    if len(dataframe["1D²"]) > 10000:
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
    fig.add_trace(go.Histogram(x=phred_score_pass, name="Read pass", marker_color="#4A69FF", histnorm='probability density'))
    fig.add_trace(go.Histogram(x=phred_score_fail, name="Read fail", marker_color="#CE3D1D", histnorm='probability density'))
    fig.add_trace(go.Scatter(x=x, y=pdf_1Dsqr_pass, mode="lines", name='Density curve of read pass', line=dict(color='#4A69FF', width=3, shape="spline", smoothing=0.5)))
    fig.add_trace(go.Scatter(x=x2, y=pdf_1Dsqr_fail, mode="lines", name='Density curve of read fail', line=dict(color='#CE3D1D', width=3, shape="spline", smoothing=0.5)))

    fig.update_layout(
        title={
            'text': "<b>1D² PHRED Score Density Distribution</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        family="Calibri, sans",
                        size=26,
                        color="black")},
        xaxis=dict(
            title="<b>PHRED score</b>",
            titlefont_size=16
        ),
        yaxis=dict(
            title='<b>Density probability</b>',
            titlefont_size=16,
            tickfont_size=14,
        ),
        legend=dict(
            x=1.02,
            y=0.95,
            title_text="<b>Legend</b>",
            title=dict(font=dict(size=16)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=15)
        ),
        barmode='group',
        hovermode='x',
        height=800, width=1400
    )

    div = py.plot(fig,
                  include_plotlyjs=False,
                  output_type='div',
                  auto_open=False,
                  show_link=False)
    py.plot(fig, filename=output_file, output_type="file", include_plotlyjs="directory", auto_open=False)

    dataframe = _make_desribe_dataframe(dataframe).drop('count')

    table_html = pd.DataFrame.to_html(dataframe)

    return main, output_file, table_html, desc, div

def scatterplot_1dsqr(result_dict, main, my_dpi, result_directory, desc):
    """
    Plot the scatter plot representing the relation between the phred score and the sequence length in log
    """
    output_file = result_directory + '/' + '_'.join(main.split())

    read_pass_length = result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.length"]
    read_pass_qscore = result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.qscore"]
    read_fail_length = result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.length"]
    read_fail_qscore = result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.qscore"]

    # If more than 10.000 reads, interpolate data
    if len(read_pass_length) > 10000:
        pass_data = _interpolate(read_pass_length, 4000, y=read_pass_qscore, interp_type="nearest")
        fail_data = _interpolate(read_fail_length, 4000, y=read_fail_qscore, interp_type="nearest")
    else:
        pass_data = [read_pass_length, read_pass_qscore]
        fail_data = [read_fail_length, read_fail_qscore]
    fig = go.Figure()

    fig.add_trace(go.Scatter(x=pass_data[0],
                             y=pass_data[1],
                             name="Pass reads",
                             marker_color="#51a96d",
                             mode="markers"
                             ))

    fig.add_trace(go.Scatter(x=fail_data[0],
                             y=fail_data[1],
                             name='Fail reads',
                             marker_color='#d90429',
                             mode="markers"
                             ))

    fig.update_layout(
        title={
            'text': "<b>Correlation between 1D² read length and PHRED score</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        family="Calibri, sans",
                        size=26,
                        color="black")},
        xaxis=dict(
            title="<b>Sequence length (bp)</b>",
            titlefont_size=16
        ),
        yaxis=dict(
            title='<b>PHRED score</b>',
            titlefont_size=16,
            tickfont_size=14,
        ),
        legend=dict(
            x=1.02,
            y=.5,
            title_text="<b>1D² Read Type</b>",
            title=dict(font=dict(size=16)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=15)
        ),
        height=800, width=1400
    )
    # Trim x axis to avoid negative values
    if max(read_pass_length) >= max(read_fail_length):
        max_val = max(read_pass_length)
    max_val = max(read_fail_length)

    fig.update_xaxes(range=[0, max_val])

    div = py.plot(fig,
                  include_plotlyjs=False,
                  output_type='div',
                  auto_open=False,
                  show_link=False)
    py.plot(fig, filename=output_file, output_type="file", include_plotlyjs="directory", auto_open=False)

    table_html = None

    return main, output_file, table_html, desc, div

#
# For each barcode 1D²
#


def barcode_percentage_pie_chart_1dsqr_pass(result_dict, dataframe_dict_1dsqr, main, barcode_selection, my_dpi, result_directory, desc):
    """
    Plots a pie chart of 1D² read pass percentage per barcode of a run.
    """
    output_file = result_directory + '/' + '_'.join(main.split())

    for element in barcode_selection:

        if all(dataframe_dict_1dsqr['barcode.arrangement'] != element):
            print("The barcode {} doesn't exist".format(element))
            return False

    count_sorted = dataframe_dict_1dsqr["read.pass.barcoded"]
    labels = count_sorted.index.values.tolist()

    fig = go.Figure(data=[go.Pie(labels=labels,
                                 values=count_sorted)])
    if len(labels) <= 12:
        palette = ["f3a683", "f7d794", "778beb", "e77f67", "cf6a87", "786fa6", "f8a5c2", "63cdda", "ea8685", "596275", "#b8e994", "#78e08f"]
        fig.update_traces(hoverinfo='label+percent', textinfo='percent', textfont_size=14,
                  marker=dict(colors=palette, line=dict(color='#2a2a2a', width=.5)))
    else:
        fig.update_traces(hoverinfo='label+percent', textinfo='percent', textfont_size=14,
                  marker=dict(line=dict(color='#2a2a2a', width=.5)))
    fig.update_traces(textposition='inside')
    fig.update_layout(uniformtext_minsize=12, uniformtext_mode='hide')
    fig.update_layout(
        title={
            'text': "<b>1D² Read Pass Barcode Distribution</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        family="Calibri, sans",
                        size=26,
                        color="black")},
        legend=dict(
            x=1.02,
            y=.5,
            title_text="<b>Barcodes</b>",
            title=dict(font=dict(size=16)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=15)
        ),
        height=500, width=875
    )

    div = py.plot(fig,
                  include_plotlyjs=False,
                  output_type='div',
                  auto_open=False,
                  show_link=False)
    py.plot(fig, filename=output_file, output_type="file", include_plotlyjs="directory", auto_open=False)

    barcode_table = pd.DataFrame({"barcode arrangement": count_sorted/sum(count_sorted)*100,
                                 "1D² read count": count_sorted})
    barcode_table.sort_index(inplace=True)
    pd.options.display.float_format = '{:.2f}%'.format
    table_html = pd.DataFrame.to_html(barcode_table)

    return main, output_file, table_html, desc, div


def barcode_percentage_pie_chart_1dsqr_fail(result_dict, dataframe_dict_1dsqr, main, barcode_selection, my_dpi, result_directory, desc):
    """
    Plots a pie chart of 1D² read fail percentage per barcode of a run.
    Needs the samplesheet file describing the barcodes to run
    """
    output_file = result_directory + '/' + '_'.join(main.split())

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
        fig.update_traces(hoverinfo='label+percent', textinfo='percent', textfont_size=14,
                  marker=dict(colors=palette, line=dict(color='#2a2a2a', width=.5)))
    else:
        fig.update_traces(hoverinfo='label+percent', textinfo='percent', textfont_size=14,
                  marker=dict(line=dict(color='#2a2a2a', width=.5)))
    fig.update_traces(textposition='inside')
    fig.update_layout(uniformtext_minsize=12, uniformtext_mode='hide')
    fig.update_layout(
        title={
            'text': "<b>1D² Read Pass Barcode Distribution</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        family="Calibri, sans",
                        size=26,
                        color="black")},
        legend=dict(
            x=1.02,
            y=.5,
            title_text="<b>Barcodes</b>",
            title=dict(font=dict(size=16)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=15)
        ),
        height=500, width=875
    )

    div = py.plot(fig,
                  include_plotlyjs=False,
                  output_type='div',
                  auto_open=False,
                  show_link=False)
    py.plot(fig, filename=output_file, output_type="file", include_plotlyjs="directory", auto_open=False)

    barcode_table = pd.DataFrame({"barcode arrangement": count_sorted/sum(count_sorted)*100,
                                  "1D² read count": count_sorted})
    barcode_table.sort_index(inplace=True)
    pd.options.display.float_format = '{:.2f}%'.format

    table_html = pd.DataFrame.to_html(barcode_table)

    return main, output_file, table_html, desc, div

def barcode_length_boxplot_1dsqr(result_dict, dataframe_dict_1dsqr, main, my_dpi, result_directory, desc):
    """
    Boxplots all the 1D² pass and fail read length for each barcode indicated in the sample sheet
    """
    output_file = result_directory + '/' + '_'.join(main.split())

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

    for col in read_pass_length.columns:
        fig.add_trace(go.Box(
            y=read_pass_length[col],
            name=col,
            marker_color='#51a96d',
            legendgroup="pass",
            offsetgroup="pass"
        ))

    for col in read_fail_length.columns:
        fig.add_trace(go.Box(
            y=read_fail_length[col],
            name=col,
            marker_color='#d90429',
            legendgroup="fail",
            offsetgroup="fail"
        ))

    fig.update_layout(
        title={
            'text': "<b>1D² Read size distribution for each barcode</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        family="Calibri, sans",
                        size=26,
                        color="black")},
        xaxis=dict(
            title="<b>Barcodes</b>",
            titlefont_size=16
        ),
        yaxis=dict(
            title='<b>Sequence length (bp)</b>',
            titlefont_size=16,
            tickfont_size=14,
        ),
        legend=dict(
            x=1.02,
            y=.5,
            title_text="<b>1D² Read Type</b>",
            title=dict(font=dict(size=16)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=15)
        ),
        boxmode='group',
        boxgap=0.4,
        boxgroupgap=0,
        height=700, width=1400
    )

    div = py.plot(fig,
                  include_plotlyjs=False,
                  output_type='div',
                  auto_open=False,
                  show_link=False)
    py.plot(fig, filename=output_file, output_type="file", include_plotlyjs="directory", auto_open=False)

    all_read = df.describe().T
    read_pass = df.loc[df['passes_filtering'] == bool(True)].describe().T
    read_fail = df.loc[df['passes_filtering'] == bool(False)].describe().T
    concat = pd.concat([all_read, read_pass, read_fail],
                       keys=['1D', '1D pass', '1D fail'])
    dataframe = concat.T

    dataframe.loc['count'] = dataframe.loc['count'].astype(int).astype(str)
    dataframe.iloc[1:] = dataframe.iloc[1:].applymap('{:.2f}'.format)
    table_html = pd.DataFrame.to_html(dataframe)

    table_html = None

    return main, output_file, table_html, desc, div

def barcoded_phred_score_frequency_1dsqr(barcode_selection, dataframe_dict_1dsqr, main, my_dpi, result_directory, desc):
    """
    Plot boxplot of the 1D pass and fail read qscore for each barcode indicated in the sample sheet
    """
    output_file = result_directory + '/' + '_'.join(main.split())

    df = dataframe_dict_1dsqr['barcode_selection_sequence_phred_melted_dataframe']
    barcode_list = barcode_selection

    # Sort reads by read type and drop read type column
    read_pass_qscore = df.loc[df['passes_filtering'] == bool(True)].drop(columns='passes_filtering')
    read_fail_qscore = df.loc[df['passes_filtering'] == bool(False)].drop(columns='passes_filtering')

    fig = go.Figure()

    for barcode in barcode_list:
        final_df = read_pass_qscore.loc[read_pass_qscore['barcodes'] == barcode].dropna()
        fig.add_trace(go.Box(
                             y=final_df['qscore'],
                             name=barcode,
                             marker_color='#51a96d',
                             legendgroup="pass",
                             offsetgroup="pass"
                             ))

    for barcode in barcode_list:
        final_df = read_fail_qscore.loc[read_fail_qscore['barcodes'] == barcode].dropna()
        fig.add_trace(go.Box(
                             y=final_df['qscore'],
                             name=barcode,
                             marker_color='#d90429',
                             legendgroup="fail",
                             offsetgroup="fail"
                             ))

    fig.update_layout(
        title={
            'text': "<b>1D² PHRED score distribution for each barcode</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        family="Calibri, sans",
                        size=26,
                        color="black")},
        xaxis=dict(
            title="<b>Barcodes</b>",
            titlefont_size=16
        ),
        yaxis=dict(
            title='<b>PHRED score</b>',
            titlefont_size=16,
            tickfont_size=14,
        ),
        legend=dict(
            x=1.02,
            y=.5,
            title_text="<b>1D² Read Type</b>",
            title=dict(font=dict(size=16)),
            bgcolor='white',
            bordercolor='white',
            font=dict(size=15)
        ),
        boxmode='group',
        boxgap=0.4,
        boxgroupgap=0,
        height=700, width=1400
    )

    div = py.plot(fig,
                  include_plotlyjs=False,
                  output_type='div',
                  auto_open=False,
                  show_link=False)
    py.plot(fig, filename=output_file, output_type="file", include_plotlyjs="directory", auto_open=False)

    all_read = df.describe().T
    read_pass = df.loc[df['passes_filtering'] == bool(True)].describe().T
    read_fail = df.loc[df['passes_filtering'] == bool(False)].describe().T
    concat = pd.concat([all_read, read_pass, read_fail], keys=['1D', '1D pass', '1D fail'])
    dataframe = concat.T
    dataframe.loc['count'] = dataframe.loc['count'].astype(int).astype(str)
    dataframe.iloc[1:] = dataframe.iloc[1:].applymap('{:.2f}'.format)
    table_html = pd.DataFrame.to_html(dataframe)

    return main, output_file, table_html, desc, div

def sequence_length_over_time_dsqr(time_df, sequence_length_df, main, my_dpi, result_directory, desc):

        output_file = result_directory + '/' + '_'.join(main.split())

        time = [t/3600 for t in time_df.dropna()]
        time = np.array(sorted(time))
        length = sequence_length_df

         # If more than 10.000 reads, interpolate data
        if len(length) > 10000:
            df_time, df_length = _interpolate(time, 100, length, "nearest")
        else:
            df_time = time
            df_length = length

        fig = go.Figure()

        fig.add_trace(go.Scatter(
        x=df_time,
        y=df_length,
        fill='tozeroy',
        fillcolor="#2f8769",
        mode='lines',
        name='interpolation curve',
        line=dict(color='#205b47', width=3, shape="spline", smoothing=0.5))
        )

        fig.update_layout(
                title={
                'text': "<b>1D² Read length over experiment time</b>",
                'y':0.95,
                'x':0,
                'xanchor': 'left',
                'yanchor': 'top',
                'font' : dict(
                family="Calibri, sans",
                size=26,
                color="black")},
            xaxis=dict(
                title="<b>Experiment time (hours)</b>",
                titlefont_size=16
                ),
            yaxis=dict(
                title='<b>Read length (bp)</b>',
                titlefont_size=16,
                tickfont_size=14,
            ),
            legend=dict(
                x=1.0,
                y=0.95,
                title_text="<b>Legend</b>",
                title=dict(font=dict(size=16)),
                bgcolor='rgba(255, 255, 255, 0)',
                bordercolor='rgba(255, 255, 255, 0)',
                font=dict(size=15)
            ),
            hovermode=False,
            height=800, width=1400
        )

        div = py.plot(fig,
                            include_plotlyjs=False,
                            output_type='div',
                            auto_open=False,
                            show_link=False)
        py.plot(fig, filename=output_file, output_type="file", include_plotlyjs="directory", auto_open=False)

        table_html = None

        return main, output_file, table_html, desc, div

def phred_score_over_time_dsqr(qscore_df, time_df, main, my_dpi, result_directory, desc):

        output_file = result_directory + '/' + '_'.join(main.split())

        # Time data
        time = [t/3600 for t in time_df.dropna()]
        time = np.array(sorted(time))

        # Qscore data
        qscore = qscore_df.dropna()

        #If more than 10.000 reads, interpolate data
        if len(qscore) > 10000:
            df_time, df_qscore = _interpolate(time, 50, qscore, "nearest")
        else:
            df_time = time
            df_score = qscore

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=df_time,
            y=df_qscore,
            fill='tozeroy',
            fillcolor="#adccf3",
            mode='lines',
            name='interpolation curve',
            line=dict(color='#7aaceb', width=3, shape="spline", smoothing=0.5),
            marker=dict(
                size=10,
                color="blue")))

        fig.update_layout(
                title={
                'text': "<b>1D² PHRED score over experiment time</b>",
                'y':0.95,
                'x':0,
                'xanchor': 'left',
                'yanchor': 'top',
                'font' : dict(
                family="Calibri, sans",
                size=26,
                color="black")},
            xaxis=dict(
                title="<b>Experiment time (hours)</b>",
                titlefont_size=16
                ),
            yaxis=dict(
                title='<b>PHRED quality score</b>',
                titlefont_size=16,
                tickfont_size=14,
            ),
            height=800, width=1400
        )

        div = py.plot(fig,
                            include_plotlyjs=False,
                            output_type='div',
                            auto_open=False,
                            show_link=False)
        py.plot(fig, filename=output_file, output_type="file", include_plotlyjs="directory", auto_open=False)
        table_html = None

        return main, output_file, table_html, desc, div

def speed_over_time_dsqr(duration_df, sequence_length_df, time_df, main, my_dpi, result_directory, desc):

        output_file = result_directory + '/' + '_'.join(main.split())

        speed = pd.Series(sequence_length_df / duration_df)

        time = [t/3600 for t in time_df]
        time = np.array(sorted(time))

        # If more than 10.000 reads, interpolate data
        if len(time) > 10000:
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
        line=dict(color='#AE3F7B', width=3, shape="linear"))
        )

        fig.update_layout(
                title={
                'text': "<b>1D² Speed over experiment time</b>",
                'y':0.95,
                'x':0,
                'xanchor': 'left',
                'yanchor': 'top',
                'font' : dict(
                family="Calibri, sans",
                size=26,
                color="black")},
            xaxis=dict(
                title="<b>Experiment time (hours)</b>",
                titlefont_size=16
                ),
            yaxis=dict(
                title='<b>Speed (bases per second)</b>',
                titlefont_size=16,
                tickfont_size=14,
            ),
            legend=dict(
                x=1.0,
                y=0.95,
                title_text="<b>Legend</b>",
                title=dict(font=dict(size=16)),
                bgcolor='rgba(255, 255, 255, 0)',
                bordercolor='rgba(255, 255, 255, 0)',
                font=dict(size=15)
            ),
            hovermode='x',
            height=800, width=1400
        )
        fig.update_yaxes(type="log")

        div = py.plot(fig,
                            include_plotlyjs=False,
                            output_type='div',
                            auto_open=False,
                            show_link=False)
        py.plot(fig, filename=output_file, output_type="file", include_plotlyjs="directory", auto_open=False)

        table_html = None

        return main, output_file, table_html, desc, div

def nseq_over_time_dsqr(time_df, main, my_dpi, result_directory, desc):

        output_file = result_directory + '/' + '_'.join(main.split())

        time = [t/3600 for t in time_df]
        time = pd.Series(time)

        # create custom xaxis points to reduce graph size
        time_points = np.linspace(min(time), max(time), 50)
        n_seq = time.groupby(pd.cut(time, time_points, right=True)).count()

        fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=time_points,
            y=list(n_seq.values), mode='lines',
            fill="tozeroy",
            fillcolor="#f4d2a7",
            line=dict(color='#edb773', width=3, shape="spline", smoothing=0.7)
        ))

        fig.update_layout(
                title={
                'text': "<b>Number of sequences through experiment time</b>",
                'y':0.95,
                'x':0,
                'xanchor': 'left',
                'yanchor': 'top',
                'font' : dict(
                family="Calibri, sans",
                size=26,
                color="black")},
            xaxis=dict(
                title="<b>Experiment time (hours)</b>",
                titlefont_size=16
                ),
            yaxis=dict(
                title='<b>Number of sequences</b>',
                titlefont_size=16,
                tickfont_size=14,
            ),
            legend=dict(
                x=1.0,
                y=0.95,
                title_text="<b>Legend</b>",
                title=dict(font=dict(size=16)),
                bgcolor='rgba(255, 255, 255, 0)',
                bordercolor='rgba(255, 255, 255, 0)',
                font=dict(size=15)
            ),
            hovermode='x',
            height=800, width=1400
        )

        div = py.plot(fig,
                            include_plotlyjs=False,
                            output_type='div',
                            auto_open=False,
                            show_link=False)
        py.plot(fig, filename=output_file, output_type="file", include_plotlyjs="directory", auto_open=False)

        table_html = None

        return main, output_file, table_html, desc, div

def _interpolate(x, npoints:int, y=None, interp_type=None, axis=-1):
    """
    Function returning an interpolated version of data passed as input
    :param x: array of data
    :param npoints: number of desired points after interpolation (int)
    :param y: second array in case of 2D data
    :param interp_type: string specifying the type of interpolation (i.e. linear, nearest, cubic, quadratic etc.)
    :param axis: number specifying the axis of y along which to interpolate. Default = -1
    """
    # In case of single array of data, use
    if y is None:
        return np.sort(resample(x, n_samples=npoints, random_state=1))
    else:
        f = interp1d(x, y, kind=interp_type, axis=axis)
        x_int = np.linspace(min(x), max(x), npoints)
        y_int = f(x_int)
        return pd.Series(x_int), pd.Series(y_int)
