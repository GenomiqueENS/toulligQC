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

# Generates a quality control report in HTML format including graphs and statistical tables

import base64
import datetime
import os
import pkgutil

from toulligqc.plotly_graph_common import (
    _format_int,
    figure_image_width,
    graph_font,
    help_html_link,
    title_size,
)


def html_report(config_dictionary, result_dict, graphs, substitutions=None):
    """
    Creation of a html report
    :param config_dictionary: dictionary containing file or directory paths
    :param result_dict: result dictionary containing all statistics
    :param graphs:
    :param substitutions: optional dict of {Measure: Value} to override in the Run statistics table
    """

    report_name = config_dictionary["report_name"]
    remove_image_files = (
        True if config_dictionary["images_directory"] is None else False
    )

    # Get report date
    report_date = _get_result_date_value(
        result_dict, "toulligqc.info.start.time", "Unknown"
    )

    # Get run date
    run_date = _get_result_date_value(
        result_dict, "sequencing.telemetry.extractor.exp.start.time", "Unknown"
    )

    sample_id = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.sample.id", "Unknown"
    )

    # Read CSS file resource
    css = pkgutil.get_data(__name__, "resources/toulligqc.css").decode("utf8")

    # Set CSS module class width to the width of the figures
    css = (
        css.replace("{figure_image_width}", str(figure_image_width) + "px")
        .replace("{title_size}", str(title_size))
        .replace("{graph_font}", str(graph_font))
    )

    # Read Plotly JavaScript code
    plotly_min_js = pkgutil.get_data(__name__, "resources/plotly-latest.min.js").decode(
        "utf8"
    )

    f = open(config_dictionary["html_report_path"], "w")

    # Create the report
    report = """<!doctype html>
<html>
  <head>
    <title>ToulligQC: {report_name} </title>
    <meta charset='UTF-8'>
    <script>{plotlyjs}</script>

    <!-- CSS stylesheet -->
    <style type="text/css">
    {css}
    </style>

  </head>

  <body>

    <!-- The banner -->
    <div id="banner">
      <div id="header_title"><img id="header_logo" alt="ToulligQC" src="{toulligqc_logo}"/>Report for {report_name}</div>
      <div id="header_filename">
        Sample ID: {sample_id} <br>
        Run date: {run_date} <br>
        Report date: {report_date} <br>
      </div>
    </div>

    <!-- The summary -->
    <div id="leftCol">
      <!--h2>Summary</h2-->
{summary_list}
    </div>

    <!-- Module results -->
    <div id="content">
{modules_report}
    </div> <!-- End of Content -->

    <!-- Footer -->
    <div id="footer"> Produced by <a href="{app_url}">{app_name}</a> (version {app_version})</div>
  </body>

</html>""".format(
        report_name=report_name,
        toulligqc_logo=_embedded_image("resources/toulligqc.png", True),
        plotlyjs=plotly_min_js,
        css=css,
        sample_id=sample_id,
        run_date=run_date,
        report_date=report_date,
        summary_list=_summary(graphs),
        modules_report=_modules_report(
            graphs,
            result_dict,
            sample_id,
            report_name,
            run_date,
            config_dictionary["app.version"],
            remove_image_files,
            substitutions=substitutions,
        ),
        app_url=config_dictionary["app.url"],
        app_name=config_dictionary["app.name"],
        app_version=config_dictionary["app.version"],
    )

    # Write the HTML page
    f.write(report)
    f.close()


def _summary(graphs):
    """
    Compose the summary section of the page
    :param graphs:
    :return: a string with HTML code for the module list
    """
    result = '        <ul class="menu-vertical">\n'
    result += (
        '          <li class="mv-item"><a href="#run_statistics'
        '">Run statistics</a></li>\n'
    )
    result += (
        '          <li class="mv-item"><a href="#software_info'
        '">Device and software</a></li>\n'
    )
    for i, t in enumerate(graphs):
        result += (
            '          <li class="mv-item"><a href="#M'
            + str(i)
            + '">'
            + t[0]
            + "</a></li>\n"
        )
    result += "        </ul>\n"
    return result


def _modules_report(
    graphs,
    result_dict,
    run_id,
    report_name,
    run_date,
    toulligqc_version,
    remove_image_files,
    substitutions=None,
):
    result = _basic_statistics_module_report(
        result_dict,
        run_id,
        report_name,
        run_date,
        toulligqc_version,
        substitutions=substitutions,
    )
    result += _other_module_reports(graphs, remove_image_files)
    return result


def _basic_statistics_module_report(
    result_dict,
    sample_id,
    report_name,
    run_date,
    toulligqc_version,
    substitutions=None,
):
    import sys as _sys

    minknow_version = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.minknow.version", "Unknown"
    )

    try:
        seconds = result_dict["basecaller.sequencing.summary.1d.extractor.run.time"]
        run_time = f"{seconds // 3600:.0f}h{(seconds % 3600) // 60:02.0f}m{seconds % 60:02.0f}s"
    except KeyError:
        run_time = "Unknown"

    read_count = result_dict["basecaller.sequencing.summary.1d.extractor.read.count"]
    run_yield = _format_int_with_prefix(
        result_dict["basecaller.sequencing.summary.1d.extractor.yield"]
    )
    n50 = result_dict["basecaller.sequencing.summary.1d.extractor.n50"]
    l50 = result_dict["basecaller.sequencing.summary.1d.extractor.l50"]

    # from telemetry file
    flow_cell_id = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.flowcell.id", "Unknown"
    )
    experiment_group = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.protocol.group.id", "Unknown"
    )
    run_id = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.run.id", "Unknown"
    )
    flowcell_version = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.flowcell.version", "Unknown"
    )
    kit_version = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.kit.version", "Unknown"
    )
    sequencing_kit_version = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.sequencing.kit.version", "Unknown"
    )
    barcode_kits_version = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.barcode.kits.version", "Unknown"
    )
    selected_speed_bases_per_second = _get_result_value(
        result_dict,
        "sequencing.telemetry.extractor.selected.speed.bases.per.second",
        "Unknown",
    )
    sample_frequency = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.sample.frequency", "Unknown"
    )
    basecaller_name = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.software.name", "Unknown"
    )
    basecaller_version = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.software.version", "Unknown"
    )
    basecaller_analysis = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.software.analysis", "Unknown"
    )
    hostname = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.hostname", "Unknown"
    )
    device_id = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.device.id", "Unknown"
    )
    device_type = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.device.type", "Unknown"
    )
    model_file = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.model.file", "Unknown"
    )
    min_qscore_threshold = _get_result_value(
        result_dict,
        "sequencing.telemetry.extractor.pass.threshold.qscore",
        value_type="float",
        default_value="Unknown",
    )

    distribution_version = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.distribution.version", "Unknown"
    )
    operating_system = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.operating.system", "Unknown"
    )
    flow_cell_product_code = _get_result_value(
        result_dict, "sequencing.telemetry.extractor.flow.cell.product.code", "Unknown"
    )
    basecalling_date = _get_result_date_value(
        result_dict, "sequencing.telemetry.extractor.basecalling.date", "Unknown"
    )

    # Build ordered list of (Measure, Value) rows for the Run statistics table
    # Using a list of tuples to preserve order and allow substitution by Measure name
    run_stats_rows = [
        ("Report name", report_name),
        ("Experiment group", experiment_group),
        ("Sample ID", sample_id),
        ("Run ID", run_id),
        ("Run date", run_date),
        ("Run duration", run_time),
        ("Flowcell ID", flow_cell_id),
        ("Flowcell product code", flow_cell_product_code),
        ("Flowcell version", flowcell_version),
        ("Kit", kit_version),
        ("Sequencing kit", sequencing_kit_version),
        ("Barcode kits", barcode_kits_version),
        ("Selected speed (bps)", selected_speed_bases_per_second),
        ("Sample frequency (Hz)", sample_frequency),
        ("Yield", run_yield),
        ("Read count", _format_int(read_count)),
        ("N50 (bp)", _format_int(int(n50))),
        ("L50", _format_int(int(l50))),
    ]

    # Apply --substitute overrides
    if substitutions:
        valid_measures = {row[0] for row in run_stats_rows}
        for measure in substitutions:
            if measure not in valid_measures:
                _sys.stderr.write(
                    f"\033[93mWarning:\033[0m --substitute measure '{measure}' "
                    f"not found in Run statistics table. "
                    f"Valid measures: {', '.join(sorted(valid_measures))}\n"
                )
        run_stats_rows = [
            (m, substitutions[m]) if m in substitutions else (m, v)
            for m, v in run_stats_rows
        ]

    from html import escape as _escape

    stats_tbody = ""
    for measure, value in run_stats_rows:
        stats_tbody += f"              <tr><th>{_escape(str(measure))}</th><td>{_escape(str(value))}</td></tr>\n"

    # Compose the Run statistics section
    result = f"""
      <div class="module" id="run_statistics">
            <h2>Run statistics {help_html_link("Run Statistics")}</h2>
            <table class="dataframe" border="">
              <thead><tr><th>Measure</th><th>Value</th></tr></thead>
              <tbody>
{stats_tbody}              </tbody>
            </table>
      </div> <!-- End of "Run-statistics" module -->
    """

    result += """
      <div class="module" id="software_info">
            <h2>Device and software {help_link}</h2>
            <table class="dataframe" border="">
                <thead><tr><th>Measure</th><th>Value</th></tr></thead>
                <tbody>
                <tr><th>Device type</th><td>{device_type}</td></tr>
                <tr><th>Device ID</th><td>{device_id}</td></tr>
                <tr><th>Device hostname</th><td>{hostname}</td></tr>
                <tr><th>Device OS</th><td>{operating_system}</td></tr>
                <tr><th>Distribution version</th><td>{distribution_version}</td></tr>
                <tr><th>MinKNOW version</th><td>{minknow_version}</td></tr>
                <tr><th>Basecaller name</th><td>{basecaller_name} </td></tr>
                <tr><th>Basecaller version</th><td>{basecaller_version}</td></tr>
                <tr><th>Basecaller analysis</th><td>{basecaller_analysis}</td></tr>
                <tr><th>Basecalling date</th><td>{basecalling_date}</td></tr>
                <tr><th>Model file</th><td>{model_file}</td></tr>
                <tr><th>Min qscore threshold</th><td>{min_qscore_threshold}</td></tr>
                <tr><th>ToulligQC version</th><td>{toulligqc_version}</td></tr>
                </tbody>
            </table>
      </div> <!-- End of "Software-info" module -->
    """.format(
        help_link=help_html_link("Software info"),
        minknow_version=minknow_version,
        basecaller_name=basecaller_name,
        basecaller_version=basecaller_version,
        basecaller_analysis=basecaller_analysis,
        basecalling_date=basecalling_date,
        toulligqc_version=toulligqc_version,
        hostname=hostname,
        operating_system=operating_system,
        distribution_version=distribution_version,
        device_type=device_type,
        device_id=device_id,
        model_file=model_file,
        min_qscore_threshold=min_qscore_threshold,
    )

    return result


def _other_module_reports(graphs, remove_image_files):
    result = ""

    for i, t in enumerate(graphs):
        if len(t) == 4:
            # Plotly Graph

            name, path, table, html = t

            # Plotly graph with table
            if table is not None:
                result += f"""
      <div class="module" id=M{i}>
        {html}
        {table}
      </div>
"""

            # Plotly graph without table
            else:
                result += f"""
      <div class="module" id=M{i}>
        {html}
      </div>
"""

        elif len(t) == 3:
            # image
            name, path, table = t

            # Image with table
            if table is not None:
                result += f"""
            <div class="module" id=M{i}>
              <h2>{name} {help_html_link(name)}</h2>
              <div class="box"><img src="{_embedded_image(path, remove=remove_image_files)}"/></div>
              {table}
            </div>
            """

            # Image without table
            else:
                result += f"""
            <div class="module" id=M{i}>
              <h2>{name} {help_html_link(name)}</h2>
              <div class="box"><img src="{_embedded_image(path, remove=remove_image_files)}"/></div>
            </div>
            """

    return result


def _embedded_image(image_path, resource=False, remove=False):
    """
    Embedded an image
    :param image_path: path of the image
    :return: a string with the image in base64
    """

    if resource:
        data = pkgutil.get_data(__name__, image_path)
    else:
        with open(image_path, "rb") as image_file:
            data = image_file.read()

    result = "data:image/png;base64," + base64.b64encode(data).decode("ascii")

    if remove:
        os.unlink(image_path)

    return result


def _get_result_value(result_dict, key, default_value="", value_type="str"):
    """
    Get the value of the result dictionary or a default value if the key does not exists.
    :param result_dict: result dictionary
    :param key: the key to use
    :param default_value: the default value
    :return: the value of key in the dictionary or the default value if the key does not exists in the dictionary
    """
    if key in result_dict:
        result = result_dict[key]
        if len(result) > 0:
            if value_type == "float":
                result = f"{float(result):.2f}"

            return result

    return default_value


def _get_result_date_value(result_dict, key, default_value=""):
    """
    Get a date value of the result dictionary and formot it. A default value is returned if the key does not exists.
    :param result_dict: result dictionary
    :param key: the key to use
    :param default_value: the default value
    :return: the value of key in the dictionary or the default value if the key does not exists in the dictionary
    """

    if key in result_dict:
        result = result_dict[key]
        if len(result) > 0:
            return _iso8601_to_formatted_date(result)

    return default_value


def _iso8601_to_formatted_date(date_string):
    """
    Format an ISO 8601 date.
    :param date_string: date to format
    :return: a formatted date
    """
    try:
        d = datetime.datetime.fromisoformat(date_string.replace("Z", "+00:00"))
    except ValueError:
        return date_string

    return d.strftime("%a %b %d %H:%M:%S %Z %Y")


def _format_int_with_prefix(i):
    for x in ((12, "T"), (9, "G"), (6, "M"), (3, "K")):
        if i / 10 ** x[0] > 1:
            return f"{float(i) / float(10 ** x[0]):.2f}{x[1]}"

    return i
