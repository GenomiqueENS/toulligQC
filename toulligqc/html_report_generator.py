# -*- coding: utf-8 -*-
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
#      https://github.com/GenomicParisCentre/toulligQC
#
# First author: Lionel Ferrato-Berberian
# Maintainer: Bérengère Laffay
# Since version 0.1

# Generates a quality control report in HTML format including graphs and statistical tables
import base64
import datetime
import os
import pkgutil

from toulligqc.plotly_graph_common import _format_int
from toulligqc.plotly_graph_common import figure_image_width
from toulligqc.plotly_graph_common import graph_font
from toulligqc.plotly_graph_common import help_html_link
from toulligqc.plotly_graph_common import title_size


def html_report(config_dictionary, result_dict, graphs):
    """
    Creation of a html report
    :param config_dictionary: dictionary containing file or directory paths
    :param result_dict: result dictionary containing all statistics
    :param graphs:
    """

    report_name = config_dictionary['report_name']
    remove_image_files = True if config_dictionary['images_directory'] is None else False

    # Get report date
    report_date = _get_result_date_value(result_dict, 'toulligqc.info.start.time', "Unknown")

    # Get run date
    run_date = _get_result_date_value(result_dict, 'sequencing.telemetry.extractor.exp.start.time', "Unknown")

    sample_id = _get_result_value(result_dict, 'sequencing.telemetry.extractor.sample.id', "Unknown")

    # Read CSS file resource
    css = pkgutil.get_data(__name__, "resources/toulligqc.css").decode('utf8')

    # Set CSS module class width to the width of the figures
    css = css.replace("{figure_image_width}", str(figure_image_width) + "px") \
        .replace("{title_size}", str(title_size)) \
        .replace("{graph_font}", str(graph_font))

    # Read Plotly JavaScript code
    plotly_min_js = pkgutil.get_data(__name__, "resources/plotly-latest.min.js").decode('utf8')

    f = open(config_dictionary['html_report_path'], 'w')

    # Create the report
    report = """<!doctype html>
<html>
  <head>
    <title>Report run MinION : {report_name} </title>
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
        Report date : {report_date} <br>
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

</html>""".format(report_name=report_name,
                  toulligqc_logo=_embedded_image("resources/toulligqc.png", True),
                  plotlyjs=plotly_min_js,
                  css=css,
                  sample_id=sample_id,
                  run_date=run_date,
                  report_date=report_date,
                  summary_list=_summary(graphs),
                  modules_report=_modules_report(graphs, result_dict, sample_id, report_name, run_date,
                                                 config_dictionary['app.version'], remove_image_files),
                  app_url=config_dictionary['app.url'],
                  app_name=config_dictionary['app.name'],
                  app_version=config_dictionary['app.version'])

    # Write the HTML page
    f.write(report)
    f.close()


def _summary(graphs):
    """
    Compose the summary section of the page
    :param graphs:
    :return: a string with HTML code for the module list
    """
    result = "        <ul class=\"menu-vertical\">\n"
    result += "          <li class=\"mv-item\"><a href=\"#run_statistics" "\">Run statistics</a></li>\n"
    result += "          <li class=\"mv-item\"><a href=\"#software_info" "\">Device and software</a></li>\n"
    for i, t in enumerate(graphs):
        result += "          <li class=\"mv-item\"><a href=\"#M" + str(i) + "\">" + t[0] + "</a></li>\n"
    result += "        </ul>\n"
    return result


def _modules_report(graphs, result_dict, run_id, report_name, run_date, toulligqc_version, remove_image_files):
    result = _basic_statistics_module_report(result_dict, run_id, report_name, run_date, toulligqc_version)
    result += _other_module_reports(graphs, remove_image_files)
    return result


def _basic_statistics_module_report(result_dict, sample_id, report_name, run_date, toulligqc_version):
    minknow_version = _get_result_value(result_dict, 'sequencing.telemetry.extractor.minknow.version', "Unknown")

    seconds = result_dict["basecaller.sequencing.summary.1d.extractor.run.time"]

    run_time = '%dh%02dm%02ds' % (seconds // 3600, (seconds % 3600) // 60, seconds % 60)

    read_count = result_dict["basecaller.sequencing.summary.1d.extractor.read.count"]
    run_yield = _format_int_with_prefix(result_dict["basecaller.sequencing.summary.1d.extractor.yield"])
    n50 = result_dict["basecaller.sequencing.summary.1d.extractor.n50"]
    l50 = result_dict["basecaller.sequencing.summary.1d.extractor.l50"]

    # from telemetry file
    flow_cell_id = _get_result_value(result_dict, 'sequencing.telemetry.extractor.flowcell.id', "Unknown")
    experiment_group = _get_result_value(result_dict, 'sequencing.telemetry.extractor.protocol.group.id', "Unknown")
    run_id = _get_result_value(result_dict, 'sequencing.telemetry.extractor.run.id', "Unknown")
    flowcell_version = _get_result_value(result_dict, 'sequencing.telemetry.extractor.flowcell.version', "Unknown")
    kit_version = _get_result_value(result_dict, 'sequencing.telemetry.extractor.kit.version', "Unknown")
    sequencing_kit_version = _get_result_value(result_dict, 'sequencing.telemetry.extractor.sequencing.kit.version', "Unknown")
    barcode_kits_version = _get_result_value(result_dict, 'sequencing.telemetry.extractor.barcode.kits.version', "Unknown")
    basecaller_name = _get_result_value(result_dict, 'sequencing.telemetry.extractor.software.name', "Unknown")
    basecaller_version = _get_result_value(result_dict, 'sequencing.telemetry.extractor.software.version', "Unknown")
    basecaller_analysis = _get_result_value(result_dict, 'sequencing.telemetry.extractor.software.analysis', "Unknown")
    hostname = _get_result_value(result_dict, 'sequencing.telemetry.extractor.hostname', "Unknown")
    device_id = _get_result_value(result_dict, 'sequencing.telemetry.extractor.device.id', "Unknown")
    device_type = _get_result_value(result_dict, 'sequencing.telemetry.extractor.device.type', "Unknown")
    model_file = _get_result_value(result_dict, 'sequencing.telemetry.extractor.model.file', "Unknown")
    min_qscore_threshold = _get_result_value(result_dict, 'sequencing.telemetry.extractor.pass.threshold.qscore',
                                             value_type='float', default_value="Unknown")

    distribution_version = _get_result_value(result_dict, 'sequencing.telemetry.extractor.distribution.version',
                                             "Unknown")
    operating_system = _get_result_value(result_dict, 'sequencing.telemetry.extractor.operating.system', "Unknown")
    flow_cell_product_code = _get_result_value(result_dict, 'sequencing.telemetry.extractor.flow.cell.product.code',
                                               "Unknown")
    basecalling_date = _get_result_date_value(result_dict, 'sequencing.telemetry.extractor.basecalling.date', "Unknown")

    # Compose the main of the page
    result = """
      <div class="module" id="run_statistics">
            <h2>Run statistics {help_link}</h2>
            <table class="dataframe" border="">
              <thead><tr><th>Measure</th><th>Value</th></tr></thead>
              <tbody>
              <tr><th>Report name </th><td>{report_name}</td></tr>
              <tr><th>Experiment group</th><td>{experiment_group}</td></tr>
              <tr><th>Sample ID</th><td>{sample_id}</td></tr>
              <tr><th>Run ID</th><td>{run_id}</td></tr>
              <tr><th>Run date</th><td>{run_date}</td></tr>
              <tr><th>Run duration </th><td>{run_time}</td></tr>
              <tr><th>Flowcell ID</th><td>{flow_cell_id}</td></tr>
              <tr><th>Flowcell product code</th><td>{flow_cell_product_code}</td></tr>
              <tr><th>Flowcell version</th><td>{flowcell_version}</td></tr>
              <tr><th>Kit</th><td>{kit_version}</td></tr>
              <tr><th>Sequencing kit</th><td>{sequencing_kit_version}</td></tr>
              <tr><th>Barcode kits</th><td>{barcode_kits_version}</td></tr>
              <tr><th>Yield</th><td>{run_yield}</td></tr>
              <tr><th>Read count</th><td>{read_count}</td></tr>
              <tr><th>N50 (bp)</th><td>{n50}</td></tr>
              <tr><th>L50</th><td>{l50}</td></tr>
              </tbody>
            </table>
      </div> <!-- End of "Run-statistics" module -->
    """.format(help_link=help_html_link("Run Statistics"),
               run_id=run_id,
               experiment_group=experiment_group,
               sample_id=sample_id,
               report_name=report_name,
               run_date=run_date,
               run_time=run_time,
               flow_cell_id=flow_cell_id,
               flow_cell_product_code=flow_cell_product_code,
               flowcell_version=flowcell_version,
               kit_version=kit_version,
               sequencing_kit_version=sequencing_kit_version,
               barcode_kits_version=barcode_kits_version,
               run_yield=run_yield,
               read_count=_format_int(read_count),
               n50=_format_int(int(n50)),
               l50=_format_int(int(l50)))

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
    """.format(help_link=help_html_link("Software info"),
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
               min_qscore_threshold=min_qscore_threshold)

    return result


def _other_module_reports(graphs, remove_image_files):
    result = ""

    for i, t in enumerate(graphs):

        if len(t) == 4:
            # Plotly Graph

            name, path, table, html = t

            # Plotly graph with table
            if table is not None:
                result += """
      <div class="module" id=M{i}>
        {html}
        {table}
      </div>
""".format(i=i, name=name, html=html, table=table)

            # Plotly graph without table
            else:
                result += """
      <div class="module" id=M{i}>
        {html}
      </div>
""".format(i=i, name=name, html=html, table=table)


        elif len(t) == 3:
            # image
            name, path, table = t

            # Image with table
            if table is not None:
                result += """
            <div class="module" id=M{i}>
              <h2>{name} {help_link}</h2>
              <div class="box"><img src="{image}"/></div>
              {table}
            </div>
            """.format(i=i, name=name, help_link=help_html_link(name),
                       image=_embedded_image(path, remove=remove_image_files), table=table)

            # Image without table
            else:
                result += """
            <div class="module" id=M{i}>
              <h2>{name} {help_link}</h2>
              <div class="box"><img src="{image}"/></div>
            </div>
            """.format(i=i, name=name, help_link=help_html_link(name),
                       image=_embedded_image(path, remove=remove_image_files))

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

    result = "data:image/png;base64," + base64.b64encode(data).decode('ascii')

    if remove:
        os.unlink(image_path)

    return result


def _get_result_value(result_dict, key, default_value="", value_type='str'):
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

            if value_type == 'float':
                result = '{:.2f}'.format(float(result))

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
        d = datetime.datetime.fromisoformat(date_string.replace('Z', '+00:00'))
    except ValueError:
        return date_string

    return d.strftime("%a %b %d %H:%M:%S %Z %Y")


def _format_int_with_prefix(i):
    for x in ((12, 'T'), (9, 'G'), (6, 'M'), (3, 'K')):
        if i / 10 ** x[0] > 1:
            return '{:.2f}{}'.format(float(i) / float(10 ** x[0]), x[1])

    return i
