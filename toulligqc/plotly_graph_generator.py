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

figure_image_width = 1000
figure_image_height = 600

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
#  1D plots
#


def read_count_histogram(result_dict, dataframe_dict, main, my_dpi, result_directory, desc):
    """
    Plots the histogram of count of the different types of reads:
    1D read return by Guppy
    1D pass read return by Guppy (Qscore >= 7)
    1D fail read return by Guppy (Qscore < 7)
    """

    output_file = result_directory + '/' + '_'.join(main.split())

    # Histogram with barcoded read counts
    if 'read.pass.barcoded.count' in dataframe_dict:

        data = {
            'Read Count': result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
            'Read Pass Count': result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
            'Read Pass Barcoded Count': dataframe_dict["read.pass.barcoded.count"],
            'Read Fail Count': result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"],
            'Read Fail Barcoded Count': dataframe_dict["read.fail.barcoded.count"]
        }
        colors = ["#54A8FC", "salmon", '#ffa931', "#50c878", "SlateBlue"]

        trace = go.Bar(x=[*data], y=list(data.values()),
                                hovertext=["<b>Total number of reads</b>",
                                           "<b>Reads of qscore > 7</b>",
                                           "<b>Barcoded reads with qscore > 7</b>",
                                           "<b>Reads of qscore < 7</b>",
                                           "<b>Barcoded reads with qscore < 7</b>"],
                                #hoverinfo="x",
                                name="Barcoded graph",
                                marker_color=colors,
                                marker_line_color="black",
                                marker_line_width=1.5, opacity=0.9)

        # Array of data for HTML table with barcode reads
        array = np.array(
            #count
            [[result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
              dataframe_dict["read.pass.barcoded.count"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"],
              dataframe_dict["read.fail.barcoded.count"]],
             #frequencies
             [result_dict["basecaller.sequencing.summary.1d.extractor.read.count.frequency"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.frequency"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.frequency"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.frequency"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.frequency"]]])

        dataframe = pd.DataFrame(array, index=['count', 'frequency'],
                                 columns=["Read count", "1D pass", "1D pass barcoded", "1D fail", "1D fail barcoded"])

    # Histogram without barcodes
    else:

        data = {
            'Read Count': result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
            'Read Pass Count': result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
            'Read Fail Count': result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"]
        }

        colors = ["#54A8FC", "salmon", "#50c878"]

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
                family="Calibri, sans",
                size=25,
                color="black")
        },
        xaxis=dict(title="<b>Read type</b>",
                   linecolor="black",
                   titlefont=dict(
                       family="Calibri",
                       size=18,
                       color="black"
                   ),
                   categoryorder="total descending"
                   ),
        yaxis=dict(title="<b>Counts</b>",
                   linecolor="black",
                   titlefont=dict(
                       family="Calibri",
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


def read_length_multihistogram(result_dict, sequence_length_df, main, my_dpi, result_directory, desc):

    output_file = result_directory + '/' + '_'.join(main.split())
    
    all_read = sequence_length_df.loc[sequence_length_df >= 10].values
    read_pass = result_dict['basecaller.sequencing.summary.1d.extractor.read.pass.length'].loc[result_dict['basecaller.sequencing.summary.1d.extractor.read.pass.length'] >= 10]
    read_fail = result_dict['basecaller.sequencing.summary.1d.extractor.read.fail.length'].loc[result_dict['basecaller.sequencing.summary.1d.extractor.read.fail.length'] >= 10]
  
    fig = go.Figure()
    
    fig.add_trace(go.Histogram(x=all_read,
                               name='All reads',
                               nbinsx=500,
                               marker_color='#fca311'  # yellow
                               ))

    fig.add_trace(go.Histogram(x=read_pass,
                               name='Pass reads',
                            nbinsx=20,
                               marker_color='#51a96d'  # green
                               ))

    fig.add_trace(go.Histogram(x=read_fail,
                            nbinsx=500,
                               name='Fail reads',
                               marker_color='#d90429'  # red
                               ))
    
    fig.update_layout(
        title={
            'text': "<b>Distribution of read length</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        family="Open Sans",
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


def yield_plot(result_dict, main, my_dpi, result_directory, desc):
    """
    Plots the different reads (1D, 1D pass, 1D fail) produced along the run against the time(in hour)
    """
    output_file = result_directory + '/' + '_'.join(main.split())
    
    all_read = result_dict['basecaller.sequencing.summary.1d.extractor.start.time.sorted']
    read_pass = result_dict['basecaller.sequencing.summary.1d.extractor.read.pass.sorted']
    read_fail = result_dict['basecaller.sequencing.summary.1d.extractor.read.fail.sorted']

    # If more than 10.000 reads, interpolate data
    if len(all_read) > 10000:
        all_read = _interpolate(x=all_read, npoints=200)
        read_pass = _interpolate(x=read_pass, npoints=200)
        read_fail = _interpolate(x=read_fail, npoints=200)
        
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(x=all_read,
                               name='All reads',
                               marker_color='#fca311',
                               mode="lines"
                               ))

    fig.add_trace(go.Scatter(x=read_pass,
                               name='Pass reads',
                               marker_color='#51a96d'
                               ))

    fig.add_trace(go.Scatter(x=read_fail,
                               name='Fail reads',
                               marker_color='#d90429'
                               ))
    
    fig.update_layout(
        title={
            'text': "<b>Yield plot through experiment time</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        family="Calibri, sans",
                        size=26,
                        color="black")},
        xaxis=dict(
            title="<b>Time (hours)</b>",
            titlefont_size=16
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


def read_quality_multiboxplot(result_dict, main, my_dpi, result_directory, desc):
    """
    Boxplot of PHRED score between read pass and read fail
    Violin plot of PHRED score between read pass and read fail
    """
    output_file = result_directory + '/' + '_'.join(main.split())

    df = pd.DataFrame(
        {"1D": result_dict["basecaller.sequencing.summary.1d.extractor.mean.qscore"],
         "1D pass": result_dict['basecaller.sequencing.summary.1d.extractor.read.pass.qscore'],
         "1D fail": result_dict['basecaller.sequencing.summary.1d.extractor.read.fail.qscore']
         })
    
    # If more than 10.000 reads, interpolate data
    if len(df["1D"]) > 10000:
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
    
    colors = {"1D": '#fca311',
              "1D pass": '#51a96d',
              "1D fail": '#d90429'}

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
            'text': "<b>PHRED score distribution of all read types</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'left',
                    'yanchor': 'top',
                    'font': dict(
                        family="Calibri, sans",
                        size=26,
                        color="black")},
        xaxis=dict(
            title="<b>Read type</b>",
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

    df = df[["1D", "1D pass", "1D fail"]]
    table_html = pd.DataFrame.to_html(_make_desribe_dataframe(df))

    return main, output_file, table_html, desc, div


def allphred_score_frequency(result_dict, main, my_dpi, result_directory, desc):
    """
    Plot the distribution of the phred score per read type (1D , 1D pass, 1D fail)
    """
    output_file = result_directory + '/' + '_'.join(main.split())

    dataframe = \
        pd.DataFrame({"1D": result_dict["basecaller.sequencing.summary.1d.extractor.mean.qscore"],
                      "1D pass": result_dict['basecaller.sequencing.summary.1d.extractor.read.pass.qscore'],
                      "1D fail": result_dict['basecaller.sequencing.summary.1d.extractor.read.fail.qscore']})
        
    # If more than 10.000 reads, interpolate data
    if len(dataframe["1D"]) > 10000:
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
    fig.add_trace(go.Histogram(x=phred_score_pass, name="Read pass", marker_color="#4A69FF", histnorm='probability density'))
    fig.add_trace(go.Histogram(x=phred_score_fail, name="Read fail", marker_color="#CE3D1D", histnorm='probability density'))
    fig.add_trace(go.Scatter(x=x, y=pdf_1D_pass, mode="lines", name='Density curve of read pass', line=dict(color='#4A69FF', width=3, shape="spline", smoothing=0.5)))
    fig.add_trace(go.Scatter(x=x2, y=pdf_1D_fail, mode="lines", name='Density curve of read fail', line=dict(color='#CE3D1D', width=3, shape="spline", smoothing=0.5)))

    fig.update_layout(
        title={
            'text': "<b>PHRED Score Density Distribution</b>",
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


def all_scatterplot(result_dict, main, my_dpi, result_directory, desc):
    """
    Plot the scatter plot representing the relation between the phred score and the sequence length in log
    """
    output_file = result_directory + '/' + '_'.join(main.split())

    read_pass_length = result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.length"]
    read_pass_qscore = result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.qscore"]
    read_fail_length = result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.length"]
    read_fail_qscore = result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.qscore"]
    
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
                             name="Fail reads",
                             marker_color="#d90429",
                             mode="markers"
                             ))

    fig.update_layout(
        title={
            'text': "<b>Correlation between read length and PHRED score</b>",
            'y': 0.95,
            'x': 0,
                    'xanchor': 'center',
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
            title_text="<b>Read Type</b>",
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


def channel_count_histogram(Guppy_log, main, my_dpi, result_directory, desc):
    """
    Plots an histogram of the channel count according to the channel number (not use anymore)
    """
    output_file = result_directory + '/' + '_'.join(main.split())
    plt.figure(figsize=(figure_image_width / my_dpi, figure_image_height / my_dpi), dpi=my_dpi)
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
    ax = plt.subplot(gs[0])
    ax.hist(Guppy_log['channel'], edgecolor='black',
            bins=range(min(Guppy_log['channel']), max(Guppy_log['channel']) + 64, 64))
    ax.set_xlabel("Channel number")
    ax.set_ylabel("Count")

    channel_count = Guppy_log['channel']
    total_number_reads_per_channel = pd.value_counts(channel_count)
    plt.subplot(gs[1])

    dataframe = table(ax, np.round(total_number_reads_per_channel
                                   .describe().drop(['mean', 'std', '50%', '75%', '25%']), 2), loc='center')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.axis('off')

    dataframe.set_fontsize(12)
    dataframe.scale(1, 1.2)

    plt.tight_layout()
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
    output_file = result_directory + '/' + '_'.join(main.split())

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
        fig.update_traces(hoverinfo='label+percent', textinfo='percent', textfont_size=14,
                  marker=dict(colors=palette, line=dict(color='#2a2a2a', width=.5)))
    else:
        fig.update_traces(hoverinfo='label+percent', textinfo='percent', textfont_size=14,
                  marker=dict(line=dict(color='#2a2a2a', width=.5)))
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
                                 "read count": count_sorted})
    barcode_table.sort_index(inplace=True)
    pd.options.display.float_format = '{:.2f}%'.format
    table_html = pd.DataFrame.to_html(barcode_table)

    return main, output_file, table_html, desc, div


def barcode_percentage_pie_chart_fail(result_dict, dataframe_dict, main, barcode_selection, my_dpi, result_directory, desc):
    """
    Plots a pie chart of 1D read fail percentage per barcode of a run.
    Needs the samplesheet file describing the barcodes to run
    """
    output_file = result_directory + '/' + '_'.join(main.split())

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
                  marker=dict(colors=palette, line=dict(color='#2a2a2a', width=.5)))
    else:
        fig.update_traces(hoverinfo='label+percent', textinfo='percent', textfont_size=14,
                  marker=dict(line=dict(color='#2a2a2a', width=.5)))
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
                                  "read count": count_sorted})
    barcode_table.sort_index(inplace=True)
    pd.options.display.float_format = '{:.2f}%'.format

    table_html = pd.DataFrame.to_html(barcode_table)

    return main, output_file, table_html, desc, div


def barcode_length_boxplot(result_dict, datafame_dict, main, my_dpi, result_directory, desc):
    """
    Boxplots all the 1D pass and fail read length for each barcode indicated in the sample sheet
    """
    output_file = result_directory + '/' + '_'.join(main.split())
    
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
        fig.add_trace(go.Box(
            y=read_pass_length[col],
            name='Barcode ' + col,
            marker_color='#51a96d',
            legendgroup="pass",
            offsetgroup="pass"
        ))

    for col in read_fail_length.columns:
        fig.add_trace(go.Box(
            y=read_fail_length[col],
            name='Barcode ' + col,
            marker_color='#d90429',
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
            title_text="<b>Read Type</b>",
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


def barcoded_phred_score_frequency(barcode_selection, dataframe_dict, main, my_dpi, result_directory, desc):
    """
    Plot boxplot of the 1D pass and fail read qscore for each barcode indicated in the sample sheet
    """
    output_file = result_directory + '/' + '_'.join(main.split())

    df = dataframe_dict['barcode_selection_sequence_phred_melted_dataframe']
    barcode_list = barcode_selection
    
    # Sort reads by read type and drop read type column 
    read_pass_qscore = df.loc[df['passes_filtering'] == bool(True)].drop(columns='passes_filtering') #Df
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
            'text': "<b>PHRED score distribution for each barcode</b>",
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
            title_text="<b>Read Type</b>",
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


def sequence_length_over_time(time_df, dataframe_dict, main, my_dpi, result_directory, desc):
        
        output_file = result_directory + '/' + '_'.join(main.split())
        
        time = [t/3600 for t in time_df.dropna()]
        time = np.array(sorted(time))

        length = dataframe_dict.get('sequence.length')

         # If more than 10.000 reads, interpolate data
        if len(length) > 10000:
            df_time, df_length = _interpolate(time, 200, length, "linear")
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
                'text': "<b>Read length over experiment time</b>",
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


def phred_score_over_time(qscore_df, time_df, main, my_dpi, result_directory, desc):
        
        output_file = result_directory + '/' + '_'.join(main.split())

        # Time data
        time = [t/3600 for t in time_df.dropna()]
        time = np.array(sorted(time))
        
        # Qscore data
        qscore = qscore_df.dropna()

        #If more than 10.000 reads, interpolate data
        if len(qscore) > 10000:
            df_time, df_qscore = _interpolate(time, 100, qscore, "nearest")
        else:
            df_time = time
            df_qscore = qscore

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=df_time,
            y=df_qscore,fill='tozeroy',
            fillcolor="#adccf3",
            mode='lines',
            name='interpolation curve',
            line=dict(color='#7aaceb', width=3, shape="spline", smoothing=0.5),
            marker=dict(
                size=10,
                color="blue")))
        
        fig.update_layout(    
                title={
                'text': "<b>PHRED score over experiment time</b>",
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


def length_over_time_slider(time_df, dataframe_dict, main, my_dpi, result_directory, desc):
        
        output_file = result_directory + '/' + '_'.join(main.split())

        time = [t/3600 for t in time_df.dropna()]
        time = np.array(sorted(time))
        
        fig = go.Figure()
        
        length = dataframe_dict.get('sequence.length')
        f = interp1d(time, length, kind="nearest")
        # Add traces, one for each slider step
        for step in range(5, 505, 5):
            x_int = np.linspace(time[0],time[-1], step)
            y_int = f(x_int)
        
            fig.add_trace(go.Scatter(
                visible=False,
                x=x_int,
                y=y_int,
                mode='lines',
                line=dict(color='#2D85E2', width=2.5, shape="spline", smoothing=0.5))
            )
        
        # Make 60th trace visible
        fig.data[60].visible = True
        
        fig.update_layout(    
                title={
                'text': "<b>Interpolated read length over experiment time</b>",
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
            height=800, width=1400
        )
        
        # Create and add slider
        npoints = []
        for i in range(len(fig.data)):
            step = dict(
                method="update",
                args=[{"visible": [False] * len(fig.data)}],  # layout attribute
            )
            step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
            npoints.append(step)

        sliders = [dict(
            active=50,
            currentvalue={"prefix": "Number of points: "},
            pad={"t": 100},
            steps=npoints
        )]

        fig.update_layout(
            sliders=sliders
        )
        
        # Edit slider labels
        fig['layout']['sliders'][0]['currentvalue']['prefix']='Number of values : '
        for i, points in enumerate(range(5, 505, 5), start = 0):
            fig['layout']['sliders'][0]['steps'][i]['label']=points
        
        div = py.plot(fig,
                            include_plotlyjs=False,
                            output_type='div',
                            auto_open=False,
                            show_link=False)
        py.plot(fig, filename=output_file, output_type="file", include_plotlyjs="directory", auto_open=False)

        table_html = None

        return main, output_file, table_html, desc, div

    
def speed_over_time(duration_df, sequence_length_df, time_df, main, my_dpi, result_directory, desc):
        
        output_file = result_directory + '/' + '_'.join(main.split())
    
        speed = pd.Series(sequence_length_df / duration_df)
        
        time = [t/3600 for t in time_df]
        time = np.array(sorted(time))

        # If more than 10.000 reads, interpolate data
        if len(time) > 10000:
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
        line=dict(color='#AE3F7B', width=3, shape="linear"))
        )

        fig.update_layout(    
                title={
                'text': "<b>Speed over experiment time</b>",
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
    

def nseq_over_time(time_df, main, my_dpi, result_directory, desc):
        
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

