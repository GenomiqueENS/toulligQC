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
# First author: Lionel Ferrato-Berberian, Karine Dias, Laurent Jourdren
# Maintainer: Karine Dias
# Since version 2.0

# This module contains common methods for plotly modules.

import numpy as np
import pandas as pd
import plotly.offline as py
from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter1d
from sklearn.utils import resample

figure_image_width = 1000
figure_image_height = 562
int_format_str = '{:,d}'
float_format_str = '{:.2f}'
percent_format_str = '{:.2f}%'
line_width = 2
interpolation_threshold = 10000

toulligqc_colors = {'all': '#fca311',  # Yellow
                    'all_1d2': '#fca311',  # Yellow
                    'pass': '#51a96d',  # Green
                    'fail': '#d90429',  # Red
                    'barcode_pass': '#51a96d',  # Green
                    'barcode_fail': '#d90429',  # Red
                    'sequence_length_over_time': '#205b47',
                    'phred_score_over_time': '#7aaceb',
                    'speed_over_time': '#AE3F7B',
                    'nseq_over_time': '#edb773',
                    }

plotly_background_color = '#e5ecf6'
legend_font_size = 16
axis_font_size = 14
on_chart_font_size = 15
title_size = 24
graph_font = 'Open sans, Helvetica, Arial, sans-serif'
image_dpi = 100


def _make_describe_dataframe(value):
    """
    Creation of a statistics table printed with the graph in report.html
    :param value: information measured (series)
    """

    desc = value.describe()
    desc.loc['count'] = desc.loc['count'].astype(int).apply(lambda x: int_format_str.format(x))
    desc.iloc[1:] = desc.iloc[1:].applymap(lambda x: float_format_str.format(x))
    desc.rename({'50%': 'median'}, axis='index', inplace=True)

    return desc


def _interpolate(x, npoints: int, y=None, interp_type=None, axis=-1):
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


def _smooth_data(npoints: int, sigma: int, data):
    """
    Function for smmothing data with numpy histogram function
    Returns a tuple of smooth data (ndarray)
    :param data: must be array-like data
    :param npoints: number of desired points for smoothing
    :param sigma: sigma value of the gaussian filter
    """
    bins = np.linspace(np.nanmin(data), np.nanmax(data), num=npoints)
    count_y, count_x = np.histogram(a=data, bins=bins, density=True)
    # Removes the first value of count_x1
    count_x = count_x[1:]
    count_y = gaussian_filter1d(count_y * len(data), sigma=sigma)
    return count_x, count_y


def _precompute_boxplot_values(y):
    """
    Precompute values for boxplot to avoid data storage in boxplot.
    """

    q1 = y.quantile(.25)
    q3 = y.quantile(.75)
    iqr = q3 - q1
    upper_fence = q3 + (1.5 * iqr)
    lower_fence = q1 - (1.5 * iqr)
    import math
    notchspan = 1.57 * iqr / math.sqrt(len(y))

    return dict(min=min(y),
                lowerfence=max(lower_fence, float(min(y))),
                q1=q1,
                median=y.quantile(.5),
                q3=q3,
                upperfence=min(upper_fence, float(max(y))),
                max=max(y),
                notchspan=notchspan)


def _dataFrame_to_html(df):
    return pd.DataFrame.to_html(df, border="")


def _transparent_colors(colors, background_color, a):
    result = []

    br = int(background_color[1:3], 16)
    bg = int(background_color[3:5], 16)
    bb = int(background_color[5:7], 16)

    for c in colors:
        r = int(c[1:3], 16)
        g = int(c[3:5], 16)
        b = int(c[5:7], 16)
        new_c = '#' + \
                _transparent_component(r, br, a) + \
                _transparent_component(g, bg, a) + \
                _transparent_component(b, bb, a)
        result.append(new_c)

    return result


def _transparent_component(c, b, a):
    v = (1 - a) * c + a * b
    r = hex(int(v))[2:]

    if len(r) == 1:
        return '0' + r
    return r


def _create_and_save_div(fig, result_directory, main):
    output_file = result_directory + '/' + '_'.join(main.split())

    div = py.plot(fig,
                  include_plotlyjs=False,
                  output_type='div',
                  auto_open=False,
                  show_link=False)
    py.plot(fig,
            filename=output_file,
            output_type="file",
            include_plotlyjs="directory",
            auto_open=False)

    return div, output_file
