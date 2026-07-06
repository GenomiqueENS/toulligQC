#
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
# For more information on the ToulligQC project and its aims,
# visit the home page at:
#
#      https://github.com/GenomiqueENS/toulligQC
#

# This module contains common methods for plotly modules.

import pkgutil
from collections import defaultdict

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as py
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from sklearn.utils import resample

figure_image_width = 1000
figure_image_height = 562
percent_format_str = "{:.2f}%"
line_width = 2

toulligqc_colors = {
    "all": "#fca311",  # Yellow
    "all_1d2": "#fca311",  # Yellow
    "pass": "#51a96d",  # Green
    "fail": "#d90429",  # Red
    "barcode_pass": "#79bf90",  # Green
    "barcode_fail": "#fb1941",  # Red
    "sequence_length_over_time": "#205b47",
    "phred_score_over_time": "#7aaceb",
    "speed_over_time": "#AE3F7B",
    "nseq_over_time": "#edb773",
    "pie_chart_palette": [
        "#f3a683",
        "#f7d794",
        "#778beb",
        "#e77f67",
        "#cf6a87",
        "#786fa6",
        "#f8a5c2",
        "#63cdda",
        "#ea8685",
        "#596275",
    ],
    "green_zone_color": "rgba(0,100,0,.1)",
}

plotly_background_color = "#e5ecf6"
legend_font_size = 16
axis_title_font_size = 14
axis_font_size = 12
on_chart_font_size = 15
title_size = 24
graph_font = "Helvetica, Arial, sans-serif"
image_dpi = 100
default_graph_layout = dict(
    font=dict(family=graph_font), height=figure_image_height, width=figure_image_width
)

interpolation_point_count_dict = {
    "read_length_distribution": (None, 10000, 3),
    "yield_plot": (None, 10000, 3),
    "phred_score_density": (None, 1000, 3),
    "over_time_graph": (None, 1000, 3),
    "scatterplot": (10000, 4000, 3),
    "phred_violin": (10000, 4000, 3),
}

help_url = "https://htmlpreview.github.io/?https://github.com/GenomiqueENS/toulligQC/master/docs/help.html"


def help_html_link(title: str, javascript: bool = True) -> str:

    if javascript:
        return '<a onclick="window.open(\'{}#{}\');" style="cursor: pointer;" id="help_link">ⓘ</a>'.format(
            help_url, title.strip().replace(" ", "_").lower()
        )

    return '<a href="{}#{}" target="_blank" id="help_link">ⓘ</a>'.format(
        help_url, title.strip().replace(" ", "_").lower()
    )


def _format_int(i: int) -> str:
    return f"{i:,d}"


def _format_float(f) -> str | int:
    try:
        s = str(f)
        i = int(s.split(".")[0])
        f = float("0." + s.split(".")[1])
    except ValueError:
        return 0

    return f"{i:,d}" + f"{f:.2f}"[1:]


def _format_percent(f: float) -> str:
    return percent_format_str.format(f)


def _title(title: str) -> dict:
    return dict(
        title=dict(
            text=f"<b>{title}</b> <b>{help_html_link(title, False)}<b>",
            y=0.95,
            x=0,
            xanchor="left",
            yanchor="top",
            font=dict(size=title_size, color="black"),
        )
    )


def _legend(legend_title: str = "Legend", args=None) -> dict:
    legend_dict = dict(
        x=1.02,
        y=0.95,
        title_text="<b>" + legend_title + "</b>",
        title=dict(font=dict(size=legend_font_size)),
        bgcolor="white",
        bordercolor="white",
        font=dict(size=legend_font_size),
    )

    if args is not None:
        legend_dict.update(dict(**args))

    return dict(legend=legend_dict)


def _xaxis(title: str, args=None) -> dict:
    axis_dict = dict(
        title="<b>" + title + "</b>",
        titlefont_size=axis_title_font_size,
        tickfont_size=axis_font_size,
    )

    if args is not None:
        axis_dict.update(dict(**args))

    return dict(xaxis=axis_dict)


def _yaxis(title: str, args=None) -> dict:
    axis_dict = dict(
        title="<b>" + title + "</b>",
        titlefont_size=axis_title_font_size,
        tickfont_size=axis_font_size,
        fixedrange=True,
    )

    if args is not None:
        axis_dict.update(dict(**args))

    return dict(yaxis=axis_dict)


def _make_describe_dataframe(value) -> pd.DataFrame:
    """Create a statistics table printed with the graph in report.html.

    Args:
        value: Information measured (a pd.Series or pd.DataFrame).

    Returns:
        A pd.DataFrame with the formatted describe statistics.
    """

    desc = value.describe()
    desc = desc.astype(object)
    desc.loc["count"] = [_format_int(int(v)) for v in desc.loc["count"]]
    desc.iloc[1:] = desc.iloc[1:].map(lambda x: _format_float(x))
    desc.rename({"50%": "median"}, axis="index", inplace=True)

    return desc


def _interpolate(x, npoints: int, y=None, interp_type=None, axis=-1):
    """Return an interpolated version of the input data.

    Args:
        x: Array of data.
        npoints: Number of desired points after interpolation.
        y: Second array in case of 2D data.
        interp_type: Type of interpolation (e.g. linear, nearest, cubic,
            quadratic).
        axis: Axis of y along which to interpolate. Defaults to -1.

    Returns:
        The interpolated data: a single sorted array when ``y`` is None, or a
        tuple ``(x_int, y_int)`` of pd.Series otherwise.
    """
    # In case of single array of data, use
    if y is None:
        return np.sort(resample(x, n_samples=npoints, random_state=1))

    else:
        f = interp1d(x, y, kind=interp_type, axis=axis)
        x_int = np.linspace(min(x), max(x), npoints)
        y_int = f(x_int)
        return pd.Series(x_int), pd.Series(y_int)


def _smooth_data(
    npoints: int,
    sigma: int,
    data,
    min_arg=None,
    max_arg=None,
    weights=None,
    density: bool = False,
) -> tuple:
    """Smooth data with the numpy histogram function.

    Args:
        npoints: Number of desired points for smoothing.
        sigma: Sigma value of the gaussian filter.
        data: Array-like data to smooth.
        min_arg: Lower bound of the histogram range; computed from ``data`` if
            None.
        max_arg: Upper bound of the histogram range; computed from ``data`` if
            None.
        weights: Optional weights passed to the histogram.
        density: Whether to normalize the histogram to a density.

    Returns:
        A tuple of smoothed data (ndarray).
    """

    if min_arg is None:
        min_arg = 0 if len(data) == 0 else np.nanmin(data)

    if max_arg is None:
        max_arg = 0 if len(data) == 0 else np.nanmax(data)

    # Compute the bin
    bins = np.linspace(min_arg, max_arg, num=npoints)

    # Compute the histogram
    y, bin_edges = np.histogram(a=data, bins=bins, weights=weights, density=density)

    # Cumulative Y
    cum_y = np.cumsum(y)

    # Center histogram
    x = bin_edges[:-1] + np.diff(bin_edges) / 2

    if min_arg == 0:
        x = np.insert(x, 0, 0)
        y = np.insert(y, 0, 0)

    if density:
        y = gaussian_filter1d(y * len(data), sigma=sigma)
        cum_y = gaussian_filter1d(cum_y * len(data), sigma=sigma)
    else:
        y = gaussian_filter1d(y, sigma=sigma)
        cum_y = gaussian_filter1d(cum_y, sigma=sigma)

    return x, y, cum_y


def _precompute_boxplot_values(y) -> dict:
    """
    Precompute values for boxplot to avoid data storage in boxplot.
    https://github.com/plotly/plotly.js/blob/master/src/traces/box/calc.js
    """

    y = y.dropna()

    if len(y) == 0:
        return dict(
            min=0, lowerfence=0, q1=0, median=0, q3=0, upperfence=0, max=0, notchspan=0
        )

    q1 = y.quantile(0.25)
    q3 = y.quantile(0.75)
    iqr = q3 - q1
    upper_fence = q3 + (1.5 * iqr)
    lower_fence = q1 - (1.5 * iqr)
    import math

    notchspan = 1.57 * iqr / math.sqrt(len(y))

    return dict(
        min=min(y),
        lowerfence=max(lower_fence, float(min(y))),
        q1=q1,
        median=y.quantile(0.5),
        q3=q3,
        upperfence=min(upper_fence, float(max(y))),
        max=max(y),
        notchspan=notchspan,
    )


def _dataFrame_to_html(df: pd.DataFrame) -> str:
    return pd.DataFrame.to_html(df, border="")


def _transparent_colors(colors: list, background_color: str, a: float) -> list:
    result = []

    br = int(background_color[1:3], 16)
    bg = int(background_color[3:5], 16)
    bb = int(background_color[5:7], 16)

    for c in colors:
        r = int(c[1:3], 16)
        g = int(c[3:5], 16)
        b = int(c[5:7], 16)
        new_c = (
            "#"
            + _transparent_component(r, br, a)
            + _transparent_component(g, bg, a)
            + _transparent_component(b, bb, a)
        )
        result.append(new_c)

    return result


def _transparent_component(c: float, b: float, a: float) -> str:
    v = (1 - a) * c + a * b
    r = hex(int(v))[2:]

    if len(r) == 1:
        return "0" + r
    return r


def _copy_latest_minjs(result_directory: str, js_file: str) -> None:
    with open(result_directory + "/" + js_file, "w+") as f:
        plotly_min_js = pkgutil.get_data(
            __name__, "resources/plotly-latest.min.js"
        ).decode("utf8")
        f.write(plotly_min_js)


def _create_and_save_div(
    fig, result_directory: str | None, main: str
) -> tuple[str, str | None]:
    div = py.plot(
        fig, include_plotlyjs=False, output_type="div", auto_open=False, show_link=False
    )

    if result_directory is not None:
        output_file = result_directory + "/" + "_".join(main.split())
        js_file = "plotly.min.js"
        py.plot(
            fig,
            filename=output_file,
            output_type="file",
            include_plotlyjs=js_file,
            auto_open=False,
        )
        _copy_latest_minjs(result_directory, js_file)
    else:
        output_file = None

    return div, output_file


def _over_time_graph(
    data_series: pd.Series,
    time_series: pd.Series,
    result_directory: str | None,
    graph_name: str,
    color: str,
    yaxis_title: str,
    log: bool = False,
    sigma: int = 1,
    quartiles: bool = True,
    min_max: bool = False,
    yaxis_starts_zero: bool = False,
    green_zone_starts_at: float | None = None,
    green_zone_color: str = "rgba(0,100,0,.1)",
) -> tuple[str, str | None, str | None, str]:
    time_bins, sigma = interpolation_points(time_series, "over_time_graph")

    t = (time_series / 3600).values
    x = np.linspace(t.min(), t.max(), num=time_bins)
    t = np.digitize(t, bins=x, right=True)

    bin_dict = defaultdict(list)
    for bin_idx, val in zip(t, data_series):
        b = x[bin_idx]
        bin_dict[b].append(val)

    percentiles = (0, 25, 50, 75, 100)
    y = []
    for i in range(len(percentiles)):
        y.append([])

    for b in x:
        if b in bin_dict:
            for i, v in enumerate(percentiles):
                y[i].append(np.percentile(bin_dict[b], v))
        else:
            for i in range(len(percentiles)):
                y[i].append(np.nan)

    for i, v in enumerate(y):
        y[i] = gaussian_filter1d(v, sigma=sigma)

    fig = go.Figure()

    # define the green zone if required
    if green_zone_starts_at is not None:
        min_x = min(time_series) / 3600
        max_x = max(time_series) / 3600
        if min_max:
            max_y = max(y[4]) * 1.05
        else:
            max_y = max(y[3]) * 1.05
        fig.add_trace(
            go.Scatter(
                mode="lines",
                x=[min_x, max_x],
                y=[max_y, max_y],
                line=dict(width=0),
                hoverinfo="skip",
                showlegend=False,
            )
        )
        fig.add_trace(
            go.Scatter(
                mode="lines",
                name="Median target",
                x=[min_x, max_x],
                y=[green_zone_starts_at, green_zone_starts_at],
                fill="tonexty",
                fillcolor=green_zone_color,
                line=dict(width=0),
                hoverinfo="skip",
            )
        )

    if quartiles:
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y[1],
                name="25% quartile",
                mode="lines",
                fill="none",
                connectgaps=True,
                line=dict(color=color, width=line_width, shape="spline"),
            )
        )

        fig.add_trace(
            go.Scatter(
                x=x,
                y=y[3],
                name="75% quartile",
                mode="lines",
                fill="tonexty",
                connectgaps=True,
                line=dict(color=color, width=line_width, shape="spline"),
            )
        )

        fig.add_trace(
            go.Scatter(
                x=x,
                y=y[2],
                name="Median",
                mode="lines",
                line=dict(color="black", width=line_width, shape="spline"),
            )
        )
        if min_max:
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y[0],
                    name="Min",
                    mode="lines",
                    line=dict(color="black", width=int(line_width / 2), shape="spline"),
                )
            )
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y[4],
                    name="Max",
                    mode="lines",
                    line=dict(color="black", width=int(line_width / 2), shape="spline"),
                )
            )

    else:
        # No quartile lines, only median
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y[2],
                name="Median",
                mode="lines",
                fill="tozeroy",
                line=dict(color=color, width=line_width, shape="spline"),
            )
        )

    # set minimal value of y axis to 0 is required
    if yaxis_starts_zero:
        range_mode = "tozero"
    else:
        range_mode = "normal"

    fig.update_layout(
        **_title(graph_name),
        **_legend(),
        **default_graph_layout,
        hovermode="x",
        **_xaxis("Experiment time (hours)"),
        **_yaxis(yaxis_title, dict(rangemode=range_mode)),
    )

    if log:
        fig.update_yaxes(type="log")

    table_html = None
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def _barcode_boxplot_graph(
    graph_name: str,
    df: pd.DataFrame,
    barcode_selection,
    pass_color: str,
    fail_color: str,
    yaxis_title: str,
    legend_title: str,
    result_directory: str | None,
    barcode_alias: dict | None = None,
) -> tuple[str, str | None, str | None, str]:
    # Sort reads by read type and drop read type column
    pass_df = df.loc[df["passes_filtering"]].drop(columns="passes_filtering")
    fail_df = df.loc[~df["passes_filtering"]].drop(columns="passes_filtering")

    fig = go.Figure()

    for read_type in ("Pass", "Fail"):
        if read_type == "Pass":
            df = pass_df
            color = pass_color
        else:
            df = fail_df
            color = fail_color

        first = True
        for barcode in sorted(barcode_selection):
            df[barcode] = df[barcode].loc[df[barcode] > 0]
            d = _precompute_boxplot_values(df[barcode])
            fig.add_trace(
                go.Box(
                    q1=[d["q1"]],
                    median=[d["median"]],
                    q3=[d["q3"]],
                    lowerfence=[d["lowerfence"]],
                    upperfence=[d["upperfence"]],
                    name=read_type + " reads",
                    x0=barcode_alias.get(barcode, barcode)
                    if barcode_alias
                    else barcode,
                    marker_color=color,
                    offsetgroup=read_type.lower(),
                    showlegend=first,
                )
            )
            if first:
                first = False

    fig.update_layout(
        **_title(graph_name),
        **_legend(legend_title),
        **default_graph_layout,
        **_xaxis("Barcodes", dict(fixedrange=True)),
        **_yaxis(yaxis_title, dict(rangemode="tozero")),
        boxmode="group",
        boxgap=0.4,
        boxgroupgap=0,
    )

    # all_read = all_df.describe().T
    # read_pass = pass_df.describe().T
    # read_fail = fail_df.describe().T
    # concat = pd.concat([all_read, read_pass, read_fail],
    #                    keys=['1D', '1D pass', '1D fail'])
    # dataframe = concat.T
    #
    # dataframe.loc['count'] = dataframe.loc['count'].astype(int).apply(lambda x: int_format_str.format(x))
    # dataframe.iloc[1:] = dataframe.iloc[1:].applymap(float_format_str.format)
    # table_html = _dataFrame_to_html(dataframe)

    table_html = None
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def _pie_chart_graph(
    graph_name: str,
    count_sorted: list,
    color_palette: list,
    one_d_square: bool,
    result_directory: str | None,
    barcode_alias: dict | None = None,
) -> tuple[str, str | None, str | None, str]:
    read_count_sorted = count_sorted[0]
    base_count_sorted = count_sorted[1]
    labels = read_count_sorted.index.values.tolist()
    if barcode_alias:
        labels = [barcode_alias.get(label, label) for label in labels]

    fig = go.Figure()

    if len(labels) <= len(color_palette):
        pie_marker = dict(
            colors=color_palette, line=dict(width=line_width, color="#808080")
        )
        bar_colors = color_palette
    else:
        pie_marker = dict(line=dict(width=line_width, color="#808080"))
        bar_colors = color_palette[0]

    # reads Pie chart
    fig.add_trace(
        go.Pie(
            labels=labels,
            values=read_count_sorted,
            hoverinfo="label+percent",
            textinfo="percent",
            textfont_size=14,
            marker=pie_marker,
            textposition="inside",
            hovertemplate="<b>%{label}</b><br>%{percent:.1%} (%{value:,})<extra></extra>",
            visible=True,
        )
    )
    # Bases Pie chart
    fig.add_trace(
        go.Pie(
            labels=labels,
            values=base_count_sorted,
            hoverinfo="label+percent",
            textinfo="percent",
            textfont_size=14,
            marker=pie_marker,
            textposition="inside",
            hovertemplate="<b>%{label}</b><br>%{percent:.1%} (%{value:,})<extra></extra>",
            visible=False,
        )
    )
    # Reads Histogram
    fig.add_trace(
        go.Bar(
            x=labels,
            y=read_count_sorted,
            marker_color=bar_colors,
            marker_line_color="gray",
            marker_line_width=line_width,
            hovertemplate="<b>%{x}</b><br>%{y:,}<extra></extra>",
            visible=False,
        )
    )
    # Bases Histogram
    fig.add_trace(
        go.Bar(
            x=labels,
            y=base_count_sorted,
            marker_color=bar_colors,
            marker_line_color="gray",
            marker_line_width=line_width,
            hovertemplate="<b>%{x}</b><br>%{y:,}<extra></extra>",
            visible=False,
        )
    )

    # Layout
    fig.update_layout(
        **_title(graph_name),
        **default_graph_layout,
        **_legend("Barcodes", args=dict(y=0.75)),
        uniformtext_minsize=12,
        uniformtext_mode="hide",
        xaxis={"visible": False},
        yaxis={"visible": False},
        plot_bgcolor="white",
    )

    # Add buttons
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="down",
                buttons=list(
                    [
                        dict(
                            args=[
                                {"visible": [True, False, False, False]},
                                {
                                    "xaxis": {"visible": False},
                                    "yaxis": {"visible": False},
                                    "plot_bgcolor": "white",
                                },
                            ],
                            label="Reads Pie chart",
                            method="update",
                        ),
                        dict(
                            args=[
                                {"visible": [False, False, True, False]},
                                {
                                    **_xaxis("Barcodes", dict(visible=True)),
                                    **_yaxis("Read count", dict(visible=True)),
                                    "plot_bgcolor": plotly_background_color,
                                },
                            ],
                            label="Reads Histogram",
                            method="update",
                        ),
                        dict(
                            args=[
                                {"visible": [False, True, False, False]},
                                {
                                    "xaxis": {"visible": False},
                                    "yaxis": {"visible": False},
                                    "plot_bgcolor": "white",
                                },
                            ],
                            label="Bases Pie chart",
                            method="update",
                        ),
                        dict(
                            args=[
                                {"visible": [False, False, False, True]},
                                {
                                    **_xaxis("Barcodes", dict(visible=True)),
                                    **_yaxis("Base count", dict(visible=True)),
                                    "plot_bgcolor": plotly_background_color,
                                },
                            ],
                            label="Bases Histogram",
                            method="update",
                        ),
                    ]
                ),
                pad={"r": 20, "t": 20, "l": 40, "b": 20},
                showactive=True,
                x=1.0,
                xanchor="left",
                y=1.25,
                yanchor="top",
            )
        ]
    )

    if one_d_square:
        count_col_name = "1D² read count"
    else:
        count_col_name = "Read count"

    barcode_table = pd.DataFrame(
        {
            "Barcode arrangement (%)": read_count_sorted / sum(read_count_sorted) * 100,
            count_col_name: read_count_sorted,
            "Base count": base_count_sorted,
        }
    )
    if barcode_alias:
        barcode_table = barcode_table.rename(index=barcode_alias)

    barcode_table.sort_index(inplace=True)
    pd.options.display.float_format = percent_format_str.format
    barcode_table[count_col_name] = (
        barcode_table[count_col_name].astype(int).apply(lambda x: _format_int(x))
    )
    barcode_table["Base count"] = (
        barcode_table["Base count"].astype(int).apply(lambda x: _format_int(x))
    )
    table_html = _dataFrame_to_html(barcode_table)

    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def _read_length_distribution(
    graph_name: str,
    all_reads,
    pass_reads,
    fail_reads,
    all_color: str,
    pass_color: str,
    fail_color: str,
    xaxis_title: str,
    result_directory: str | None,
    bin_width: float | None = None,
) -> tuple[str, str | None, str | None, str]:
    npoints, sigma = interpolation_points(all_reads, "read_length_distribution")
    min_all_reads = min(all_reads)
    max_all_reads = max(all_reads)

    if bin_width is not None:
        read_span = max_all_reads - min_all_reads
        if read_span > 0:
            npoints = int(np.ceil(read_span / bin_width)) + 1
            npoints = max(2, npoints)

    count_x1, count_y1, cum_count_y1 = _smooth_data(
        npoints=npoints,
        sigma=sigma,
        data=all_reads,
        min_arg=min_all_reads,
        max_arg=max_all_reads,
    )
    count_x2, count_y2, cum_count_y2 = _smooth_data(
        npoints=npoints,
        sigma=sigma,
        data=pass_reads,
        min_arg=min_all_reads,
        max_arg=max_all_reads,
    )
    count_x3, count_y3, cum_count_y3 = _smooth_data(
        npoints=npoints,
        sigma=sigma,
        data=fail_reads,
        min_arg=min_all_reads,
        max_arg=max_all_reads,
    )

    sum_x1, sum_y1, cum_sum_y1 = _smooth_data(
        npoints=npoints,
        sigma=sigma,
        data=all_reads,
        weights=all_reads,
        min_arg=min_all_reads,
        max_arg=max_all_reads,
    )
    sum_x2, sum_y2, cum_sum_y2 = _smooth_data(
        npoints=npoints,
        sigma=sigma,
        data=pass_reads,
        weights=pass_reads,
        min_arg=min_all_reads,
        max_arg=max_all_reads,
    )
    sum_x3, sum_y3, cum_sum_y3 = _smooth_data(
        npoints=npoints,
        sigma=sigma,
        data=fail_reads,
        weights=fail_reads,
        min_arg=min_all_reads,
        max_arg=max_all_reads,
    )

    # Find 50 percentile for zoomed range on x axis
    max_x_range = np.percentile(all_reads, 99)

    coef = max_all_reads / npoints

    max_y = max(max(count_y1), max(count_y2), max(count_y3)) / coef
    max_sum_y = max(max(sum_y1), max(sum_y2), max(sum_y3)) / coef

    fig = go.Figure()

    # Read graphs
    fig.add_trace(
        go.Scatter(
            x=count_x1,
            y=count_y1 / coef,
            name="All reads",
            fill="tozeroy",
            marker_color=all_color,
            visible=True,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=count_x2,
            y=count_y2 / coef,
            name="Pass reads",
            fill="tozeroy",
            marker_color=pass_color,
            visible=True,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=count_x3,
            y=count_y3 / coef,
            name="Fail reads",
            fill="tozeroy",
            marker_color=fail_color,
            visible=True,
        )
    )

    # Threshold
    for p in [25, 50, 75]:
        x0 = np.percentile(all_reads, p)
        if p == 50:
            t = "median<br>all reads"
        else:
            t = str(p) + "%<br>all reads"
        fig.add_trace(
            go.Scatter(
                mode="lines+text",
                name="All reads",
                x=[x0, x0],
                y=[0, max_y],
                line=dict(color="gray", width=1, dash="dot"),
                text=["", t],
                textposition="top center",
                hoverinfo="skip",
                showlegend=False,
                visible=True,
            )
        )

    # Base plots
    # Read graphs
    fig.add_trace(
        go.Scatter(
            x=sum_x1,
            y=sum_y1 / coef,
            name="All reads",
            fill="tozeroy",
            marker_color=all_color,
            visible=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=sum_x2,
            y=sum_y2 / coef,
            name="Pass reads",
            fill="tozeroy",
            marker_color=pass_color,
            visible=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=sum_x3,
            y=sum_y3 / coef,
            name="Fail reads",
            fill="tozeroy",
            marker_color=fail_color,
            visible=False,
        )
    )

    # Threshold
    for p in [25, 50, 75]:
        x0 = np.percentile(all_reads, p)
        if p == 50:
            t = "median<br>all reads"
        else:
            t = str(p) + "%<br>all reads"
        fig.add_trace(
            go.Scatter(
                mode="lines+text",
                name="All reads",
                x=[x0, x0],
                y=[0, max_y],
                line=dict(color="gray", width=1, dash="dot"),
                text=["", t],
                textposition="top center",
                hoverinfo="skip",
                showlegend=False,
                visible=False,
            )
        )

    fig.update_layout(
        **_title(graph_name),
        **default_graph_layout,
        **_legend(args=dict(y=0.75)),
        hovermode="x",
        **_xaxis(xaxis_title, dict(range=[min_all_reads, max_x_range], type="linear")),
        **_yaxis("Read count", dict(range=[0, max_y * 1.10])),
    )

    # Add buttons
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="down",
                buttons=list(
                    [
                        dict(
                            args=[
                                {
                                    "visible": [
                                        True,
                                        True,
                                        True,
                                        True,
                                        True,
                                        True,
                                        False,
                                        False,
                                        False,
                                    ]
                                },
                                {
                                    "xaxis": {
                                        "type": "linear",
                                        "range": [min_all_reads, max_x_range],
                                    },
                                    "yaxis": {
                                        "title": "<b>Read count</b>",
                                        "range": [0, max_y * 1.10],
                                    },
                                },
                            ],
                            label="Reads linear",
                            method="update",
                        ),
                        dict(
                            args=[
                                {
                                    "visible": [
                                        True,
                                        True,
                                        True,
                                        True,
                                        True,
                                        True,
                                        False,
                                        False,
                                        False,
                                    ]
                                },
                                {
                                    "xaxis": {"type": "log"},
                                    "yaxis": {
                                        "title": "<b>Read count</b>",
                                        "range": [0, max_y * 1.10],
                                    },
                                },
                            ],
                            label="Reads log",
                            method="update",
                        ),
                        dict(
                            args=[
                                {
                                    "visible": [
                                        False,
                                        False,
                                        False,
                                        False,
                                        False,
                                        False,
                                        True,
                                        True,
                                        True,
                                    ]
                                },
                                {
                                    "xaxis": {
                                        "type": "linear",
                                        "range": [min_all_reads, max_x_range],
                                    },
                                    "yaxis": {
                                        "title": "<b>Base count</b>",
                                        "range": [0, max_sum_y * 1.10],
                                    },
                                },
                            ],
                            label="Bases linear",
                            method="update",
                        ),
                        dict(
                            args=[
                                {
                                    "visible": [
                                        False,
                                        False,
                                        False,
                                        False,
                                        False,
                                        False,
                                        True,
                                        True,
                                        True,
                                    ]
                                },
                                {
                                    "xaxis": {"type": "log"},
                                    "yaxis": {
                                        "title": "<b>Base count</b>",
                                        "range": [0, max_sum_y * 1.10],
                                    },
                                },
                            ],
                            label="Bases log",
                            method="update",
                        ),
                    ]
                ),
                pad={"r": 20, "t": 20, "l": 20, "b": 20},
                showactive=True,
                x=1.0,
                xanchor="left",
                y=1.25,
                yanchor="top",
            ),
        ]
    )

    # Create data for HTML table
    table_df = pd.concat(
        [pd.Series(all_reads), pass_reads, fail_reads],
        axis=1,
        keys=["All reads", "Pass reads", "Fail reads"],
    )
    table_html = _dataFrame_to_html(_make_describe_dataframe(table_df))

    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def _phred_score_density(
    graph_name: str,
    dataframe: pd.DataFrame,
    prefix: str,
    all_color: str,
    pass_color: str,
    fail_color: str,
    result_directory: str | None,
) -> tuple[str, str | None, str | None, str]:
    all_series = dataframe[prefix].dropna()
    pass_series = dataframe[prefix + " pass"].dropna()
    fail_series = dataframe[prefix + " fail"].dropna()

    npoints, sigma = interpolation_points(all_series, "phred_score_density")

    count_x2, count_y2, cum_count_y2 = _smooth_data(
        npoints=npoints,
        sigma=sigma,
        data=pass_series,
        min_arg=np.nanmin(all_series),
        max_arg=np.nanmax(all_series),
        density=True,
    )
    count_x3, count_y3, cum_count_y3 = _smooth_data(
        npoints=npoints,
        sigma=sigma,
        data=fail_series,
        min_arg=np.nanmin(all_series),
        max_arg=np.nanmax(all_series),
        density=True,
    )

    count_y2 = count_y2 / len(all_series)
    count_y3 = count_y3 / len(all_series)

    max_y = max(max(count_y2), max(count_y3))

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=count_x2,
            y=count_y2,
            name="Pass reads",
            fill="tozeroy",
            marker_color=pass_color,
            visible=True,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=count_x3,
            y=count_y3,
            name="Fail reads",
            fill="tozeroy",
            marker_color=fail_color,
            visible=True,
        )
    )

    # Threshold
    for p in [25, 50, 75]:
        x0 = np.percentile(pass_series, p)
        if p == 50:
            t = "median"
        else:
            t = str(p) + "%"
        fig.add_trace(
            go.Scatter(
                mode="lines+text",
                name="Pass read<br>percentiles",
                x=[x0, x0],
                y=[0, max_y],
                line=dict(color="gray", width=1, dash="dot"),
                text=["", t],
                textposition="top center",
                hoverinfo="skip",
                showlegend=(True if p == 50 else False),
                visible=True,
            )
        )

    fig.update_layout(
        **_title(graph_name),
        **_legend(),
        **default_graph_layout,
        hovermode="x",
        **_xaxis("PHRED score", dict(rangemode="tozero")),
        **_yaxis("Density probability", dict(rangemode="tozero")),
    )

    table_html = None
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def _quality_multiboxplot(
    graph_name: str,
    result_directory: str | None,
    df: pd.DataFrame,
    onedsquare: bool = False,
) -> tuple[str, str | None, str | None, str]:
    if onedsquare:
        prefix = "1D²"
    else:
        prefix = "1D"

    # If more than 10.000 reads, interpolate data
    npoints = interpolation_points(df[prefix], "phred_violin")[0]
    if len(df[prefix]) != npoints:
        violin_df = pd.DataFrame(
            {
                prefix: _interpolate(df[prefix], npoints),
                prefix + " pass": _interpolate(df[prefix + " pass"], npoints),
                prefix + " fail": _interpolate(df[prefix + " fail"], npoints),
            }
        )
    else:
        violin_df = df
    names = {
        prefix: "All reads",
        prefix + " pass": "Pass reads",
        prefix + " fail": "Fail reads",
    }

    colors = {
        prefix: toulligqc_colors["all"],
        prefix + " pass": toulligqc_colors["pass"],
        prefix + " fail": toulligqc_colors["fail"],
    }

    # Max yaxis value for displaying same scale between plots
    max_yaxis = (
        max(
            df.max(skipna=True, numeric_only=True).values.max(),
            violin_df.max(skipna=True, numeric_only=True).values.max(),
        )
        + 2.0
    )
    min_yaxis = (
        min(
            df.min(skipna=True, numeric_only=True).values.min(),
            violin_df.min(skipna=True, numeric_only=True).values.min(),
        )
        - 2.0
    )

    fig = go.Figure()

    for column in df.columns:
        d = _precompute_boxplot_values(df[column])
        fig.add_trace(
            go.Box(
                q1=[d["q1"]],
                median=[d["median"]],
                q3=[d["q3"]],
                lowerfence=[d["lowerfence"]],
                upperfence=[d["upperfence"]],
                name=names[column],
                x0=names[column],
                marker=dict(opacity=0.3, color=colors[column]),
                boxmean=False,
                showlegend=True,
            )
        )

        fig.add_trace(
            go.Violin(
                y=violin_df[column],
                name=names[column],
                meanline_visible=True,
                marker=dict(color=colors[column]),
                visible=False,
            )
        )

    fig.update_layout(
        **_title(graph_name),
        **default_graph_layout,
        **_legend(),
        hovermode="x",
        **_xaxis("Read type", dict(fixedrange=True)),
        **_yaxis("PHRED score", dict(range=[min_yaxis, max_yaxis])),
    )

    # Add buttons
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="left",
                buttons=list(
                    [
                        dict(
                            args=[{"visible": [True, False]}, {"hovermode": "x"}],
                            label="Boxplot",
                            method="update",
                        ),
                        dict(
                            args=[{"visible": [False, True]}, {"hovermode": False}],
                            label="Violin plot",
                            method="update",
                        ),
                    ]
                ),
                pad={"r": 20, "t": 20, "l": 20, "b": 20},
                showactive=True,
                x=1.0,
                xanchor="left",
                y=1.25,
                yanchor="top",
            ),
        ]
    )

    df = df[[prefix, prefix + " pass", prefix + " fail"]]
    df.columns = ["All reads", "Pass reads", "Fail reads"]
    table_html = _dataFrame_to_html(_make_describe_dataframe(df))

    div, output_file = _create_and_save_div(fig, result_directory, graph_name)
    return graph_name, output_file, table_html, div


def _twod_density_char(
    graph_name: str,
    dataframe_dict: dict,
    result_directory: str | None,
    onedsquare: bool = False,
) -> tuple[str, str | None, str | None, str]:
    read_pass_length = dataframe_dict["pass.reads.sequence.length"]
    read_pass_qscore = dataframe_dict["pass.reads.mean.qscore"]
    read_fail_length = dataframe_dict["fail.reads.sequence.length"]
    read_fail_qscore = dataframe_dict["fail.reads.mean.qscore"]

    all_length = dataframe_dict["all.reads.sequence.length"]
    all_qscore = dataframe_dict["all.reads.mean.qscore"]

    prefix = "1D² " if onedsquare else ""
    graph_name = prefix + graph_name
    npoint = 50000

    fig = go.Figure()

    idx_pass = np.random.choice(
        read_pass_length.index, min(npoint, len(read_pass_length)), replace=False
    )
    idx_fail = np.random.choice(
        read_fail_length.index, min(npoint, len(read_fail_length)), replace=False
    )
    idx_all = np.random.choice(
        all_length.index, min(npoint, len(all_length)), replace=False
    )

    pass_color = toulligqc_colors["pass"]
    fail_color = toulligqc_colors["fail"]
    all_color = toulligqc_colors["all"]

    empty_fail = len(read_fail_length) == 0

    fig.add_trace(
        go.Histogram2dContour(
            x=all_length[idx_all],
            y=all_qscore[idx_all],
            colorscale=[[0, "white"], [0.5, "khaki"], [1.0, all_color]],
            reversescale=False,
            xaxis="x",
            yaxis="y",
            colorbar=dict(title="<b>Legend</b>", len=0.5),
        )
    )

    fig.add_trace(
        go.Histogram2dContour(
            x=read_pass_length[idx_pass],
            y=read_pass_qscore[idx_pass],
            colorscale=[[0, "white"], [0.5, "honeydew"], [1.0, pass_color]],
            reversescale=False,
            xaxis="x",
            yaxis="y",
            colorbar=dict(
                # title = '<b>Legend</b>',
                len=0.4,
                title=dict(font=dict(size=9)),
            ),
            visible=False,
        )
    )

    fig.add_trace(
        go.Histogram2dContour(
            x=read_fail_length[idx_fail],
            y=read_fail_qscore[idx_fail],
            colorscale=[[0, "white"], [0.5, "coral"], [1.0, fail_color]],
            reversescale=False,
            xaxis="x",
            yaxis="y",
            colorbar=dict(
                title="<b>Legend</b>",
                len=0.4,
            ),
            visible=False,
        )
    )

    max_x_range = np.percentile(all_length, 99)
    max_y_range = np.percentile(all_qscore, 99.8)
    fig.update_xaxes(range=[min(all_length), max_x_range])
    fig.update_yaxes(range=[min(all_qscore), max_y_range])

    fig.add_trace(
        go.Histogram(
            y=all_qscore[idx_all], xaxis="x2", opacity=0.5, marker=dict(color=all_color)
        )
    )

    fig.add_trace(
        go.Histogram(
            x=all_length[idx_all], yaxis="y2", opacity=0.5, marker=dict(color=all_color)
        )
    )

    fig.add_trace(
        go.Histogram(
            y=read_pass_qscore[idx_pass],
            xaxis="x2",
            marker=dict(color=pass_color),
            visible=False,
        )
    )

    fig.add_trace(
        go.Histogram(
            x=read_pass_length[idx_pass],
            yaxis="y2",
            marker=dict(color=pass_color),
            visible=False,
        )
    )
    fig.add_trace(
        go.Histogram(
            y=read_fail_qscore[idx_fail],
            xaxis="x2",
            marker=dict(color=fail_color),
            visible=False,
        )
    )
    fig.add_trace(
        go.Histogram(
            x=read_fail_length[idx_fail],
            yaxis="y2",
            marker=dict(color=fail_color),
            visible=False,
        )
    )

    fig.update_layout(
        autosize=False,
        xaxis=dict(
            zeroline=False,
            domain=[0, 0.85],
            showgrid=False,
            title="<b>Sequence length</b>",
        ),
        yaxis=dict(
            zeroline=False, domain=[0, 0.85], showgrid=False, title="<b>PHRED score</b>"
        ),
        xaxis2=dict(zeroline=False, domain=[0.85, 1], showgrid=False, visible=False),
        yaxis2=dict(zeroline=False, domain=[0.85, 1], showgrid=False, visible=False),
        height=600,
        width=1000,
        bargap=0,
        paper_bgcolor="#FFFFFF",
        plot_bgcolor="#FFFFFF",
        hovermode="closest",
        showlegend=False,
    )

    # Add buttons
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="down",
                buttons=list(
                    [
                        dict(
                            args=[
                                {
                                    "visible": [
                                        True,
                                        False,
                                        False,
                                        True,
                                        True,
                                        False,
                                        False,
                                        False,
                                        False,
                                    ]
                                },
                                {"hovermode": False},
                            ],
                            label="All reads",
                            method="restyle",
                        ),
                        dict(
                            args=[
                                {
                                    "visible": [
                                        False,
                                        True,
                                        False,
                                        False,
                                        False,
                                        True,
                                        True,
                                        False,
                                        False,
                                    ]
                                },
                                {"hovermode": "x"},
                            ],
                            label="Pass reads",
                            method="restyle",
                        ),
                        dict(
                            args=[
                                {
                                    "visible": [
                                        False,
                                        False,
                                        True,
                                        False,
                                        False,
                                        False,
                                        False,
                                        True,
                                        True,
                                    ]
                                },
                                {"hovermode": False},
                            ],
                            label="Fail reads",
                            method="restyle",
                        )
                        if not empty_fail
                        else dict(),
                    ]
                ),
                pad={"r": 20, "t": 20, "l": 20, "b": 20},
                showactive=True,
                x=1.0,
                xanchor="left",
                y=1.25,
                yanchor="top",
            ),
        ]
    )

    fig.update_layout(
        **_title(graph_name),
        **_legend(),
        **default_graph_layout,
        xaxis_rangeselector_font_color="black",
        xaxis_rangeselector_activecolor="red",
        xaxis_rangeselector_bgcolor="green",
        # hovermode='x',
        # **_xaxis('PHRED score', dict(rangemode="tozero")),
        # **_yaxis('Density probability', dict(rangemode="tozero")),
    )

    table_html = None
    div, output_file = _create_and_save_div(fig, result_directory, graph_name)

    return graph_name, output_file, table_html, div


def interpolation_points(series, graph_name: str) -> tuple[int, int]:
    count = len(series)
    threshold, npoints, sigma = interpolation_point_count_dict[graph_name]

    if threshold is not None:
        if count > threshold:
            result = npoints
        else:
            result = count
    else:
        result = npoints

    return result, sigma
