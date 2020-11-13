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


def html_report(config_dictionary, result_dict, graphs):
    """
    Creation of a html report
    :param config_dictionary: dictionary containing file or directory paths
    :param result_dict: result dictionary containing all statistics
    :param graphs:
    """

    result_directory = config_dictionary['result_directory']
    report_name = config_dictionary['report_name']

    # from sequence summary file

    td = datetime.timedelta(hours=result_dict["basecaller.sequencing.summary.1d.extractor.run.time"])
    seconds = td.total_seconds()
    run_time = '%d:%02d:%02d' % (seconds / 3600, seconds / 60 % 60, seconds % 60)

    report_date = result_dict['toulligqc.info.start.time']

    # from Fast5 file
    run_date = result_dict['sequencing.telemetry.extractor.exp.start.time']
    flow_cell_id = result_dict['sequencing.telemetry.extractor.flowcell.id']
    run_id = result_dict['sequencing.telemetry.extractor.sample.id']
    minknow_version = result_dict['sequencing.telemetry.extractor.minknow.version']

    read_count = result_dict["basecaller.sequencing.summary.1d.extractor.read.count"]
    run_yield = round(result_dict["basecaller.sequencing.summary.1d.extractor.yield"]/1000000000, 2)
    n50 = result_dict["basecaller.sequencing.summary.1d.extractor.n50"]

    # from telemetry file
    flowcell_version = _get_result_value(result_dict, 'sequencing.telemetry.extractor.flowcell.version', "Unknown")
    kit_version = _get_result_value(result_dict, 'sequencing.telemetry.extractor.kit.version', "Unknown")
    basecaller_name = _get_result_value(result_dict, 'sequencing.telemetry.extractor.software.name', "Unknown")
    basecaller_version = _get_result_value(result_dict, 'sequencing.telemetry.extractor.software.version', "Unknown")
    basecaller_analysis = _get_result_value(result_dict, 'sequencing.telemetry.extractor.software.analysis', "Unknown")
    hostname = _get_result_value(result_dict, 'sequencing.telemetry.extractor.hostname', "Unknown")
    device_id = _get_result_value(result_dict, 'sequencing.telemetry.extractor.device.id', "Unknown")
    device_type = _get_result_value(result_dict, 'sequencing.telemetry.extractor.device.type', "Unknown")
    model_file = _get_result_value(result_dict, 'sequencing.telemetry.extractor.model.file', "Unknown")
    sample_id = _get_result_value(result_dict, 'sequencing.telemetry.extractor.sample.id', "Unknown")

    f = open(result_directory + 'report.html', 'w')

# Define the header of the page
    title = """<!doctype html>
<html>
  <head>
    <title>Report run MinION : {0} </title>
    <meta charset='UTF-8'>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style type="text/css">
    """.format(report_name)

    header = """

  @media screen {

    div.summary {
      width: 16em;
      position:fixed;
      top: 4.5em;
      margin:1em 0 0 1em;
    }

    div.main {
      display:block;
      position:absolute;
      overflow:auto;
      height:auto;
      width:auto;
      top:4.5em;
      bottom:2.3em;
      left:18em;
      right:0;
      border-left: 1px solid #CCC;
      padding:0 0 0 1em;
      background-color: white;
      z-index:1;
    }

    div.header {
      background-color: #EEE;
      border:0;
      margin:0;
      padding: 0.5em;
      font-size: 200%;
      font-weight: bold;
      position:fixed;
      width:100%;
      top:0;
      left:0;
      z-index:2;
    }

    div.footer {
      background-color: #EEE;
      border:0;
      margin:0;
      padding:0.5em;
      height: 1.3em;
      overflow:hidden;
      font-size: 100%;
      font-weight: bold;
      position:fixed;
      bottom:0;
      width:100%;
      z-index:2;
    }

    img.indented {
      margin-left: 3em;
    }
  }

  @media print {

    img {
      max-width:80% !important;
      page-break-inside: avoid;
    }

    h2, h3 {
      page-break-after: avoid;
    }

    div.header {
      background-color: #FFF;
    }

  }

  body {
    font-family: sans-serif;
    color: #000;
    background-color: #FFF;
    border: 0;
    margin: 0;
    padding: 0;
  }

  div.header {
    border:0;
    margin:0;
    padding: 0.5em;
    font-size: 200%;
    font-weight: bold;
    width:100%;
  }

  #header_title {
    display:inline-block;
    float:left;
    clear:left;
  }

  #header_filename {
    display:inline-block;
    float:right;
    clear:right;
    font-size: 50%;
    margin-right:2em;
    text-align: right;
  }

  div.header h3 {
    font-size: 50%;
    margin-bottom: 0;
  }

  div.summary ul {
    padding-left:0;
    list-style-type:none;
  }

  div.summary ul li img {
    margin-bottom:-0.5em;
    margin-top:0.5em;
  }

  div.main {
    background-color: white;
  }

  div.module {
    padding-bottom:1.5em;
    padding-top:1.5em;
  }
  
  .info-box {
    float:left;
    min-width: 400px;
    height: 350px;
    margin: 0em;
    padding-bottom:1.5em;
    padding-top:0em;
    top: 0 auto;
    bottom: 0 auto;
  }
  
  .info-box-left {
    float:left;
    min-width: 400px;
    height: 350px;
    margin: 0em;
    padding-bottom:1.5em;
    padding-top:0em;
    top: 0 auto;
    bottom: 0 auto;
  }
  
  .box {
    float:left;
    min-width: 1150px;
    margin: 0em;
    padding-bottom:1.5em;
    padding-top:1.8em;
    top: 0 auto;
    bottom: 0 auto;
  }
  
    .box-left {
    float:left;
    min-width: 350px;
    height: 350px;
    margin: 0em;
    padding-bottom:1.5em;
    padding-top:1.5em;
    top: 0.3px;
    bottom: 0 auto;
  }
  
  .after-box {
    clear: left;
    padding-bottom: 90px;
  }

  div.footer {
    background-color: #EEE;
    border:0;
    margin:0;
    padding: 0.5em;
    font-size: 100%;
    font-weight: bold;
    width:100%;
  }


  a {
    color: #000080;
  }

  a:hover {
    color: #800000;
  }

  h2 {
    color: #000000;
    padding-bottom: 0;
    margin-bottom: 0;
    clear:left;
  }

  table {
    margin-left: 3em;
    text-align: center;
    border-collapse:collapse;
  }

  th {
    text-align: center;
    background-color: #244C89;
    color: #FFF;
    padding: 0.4em;
  }

  td {
    font-family: monospace;
    text-align: right;
    background-color: #EFEFEF;
    color: #000;
    padding: 0.4em;
  }

  img {
    padding-top: 0;
    margin-top: 0;
    border-top: 0;
  }
  
  p {
    padding-top: 0;
    margin-top: 0;
  }

    </style>
  </head>
"""

    # Define the footer of the page
    footer = """
  <body>
    <div class="footer"> Produced by <a href="{0}">{1}</a> (version {2})</div>
  </body>

</html>
""".format(config_dictionary['app.url'], config_dictionary['app.name'], config_dictionary['app.version'])

    # Compose the Banner of the page
    banner = """
    <div class="header">
      <div id="header_title">ToulligQC report for {0} <br/></div>
      <div id="header_filename">
        Run id: {0} <br>
        Report name: {1} <br>
        Run date: {2} <br>
        Report date : {3} <br>
      </div>
    </div>
""".format(run_id, report_name, run_date, report_date)

    # Compose the summary section of the page
    summary = """
    <div class='summary'>
      <h2>Summary</h2>
      <ol>
"""
    summary += "<li><a href=\"#Basic-statistics" "\"> Basic Statistics </a></li>\n"
    for i, t in enumerate(graphs):
        summary += "        <li><a href=\"#M" + str(i) + "\">" + t[0] + "</a></li>\n"
    summary += """      </ol>
    </div>
"""
    # Compose the main of the page
    main_report = """
    <div class = 'main'>
      <div class=\"module\" id="Basic-statistics">
        <div class = "info-box"> 
            <h2 id=M{0}>Basic Statistics</h2>
            <h2 id=M{0}></h2>
            <h3 style="text-align:center">Run info</h3>
            <table class=\" dataframe \" border="1">
              <thead><tr><th>Measure</th><th>Value</th></tr></thead>      
              <tbody>
              <tr><th>Run id</th><td> {0} </td></tr>  
              <tr><th>Sample</th><td> {1} </td></tr>                      
              <tr><th>Report name </th><td> {2} </td></tr>
              <tr><th>Run date</th><td> {3} </td></tr>
              <tr><th>Run duration </th><td> {4} </td></tr>
              <tr><th>Flowcell id </th><td> {5} </td></tr>
              <tr><th>Flowcell version</th><td> {6} </td></tr>
              <tr><th>Kit</th><td> {7} </td></tr>
              <tr><th>Yield (Gbp)</th><td> {8} </td></tr>
              <tr><th>Read count</th><td> {9} </td></tr>
              <tr><th>N50 (bp)</th><td> {10} </td></tr>
              </tbody>
            </table> 
        </div> <!-- end .info-box -->
      </div>
    """.format(run_id,sample_id, report_name, run_date, run_time, flow_cell_id, flowcell_version, kit_version, run_yield, read_count, n50)

    main_report += """
      <div class=\"module\">
        <div class = "info-box-left">
            <h2 id=M{0}></h2>
            <h2 id=M{0}></h2>
            <h3 style="text-align:center">Software info</h3>
            <table class=\" dataframe \" border="1">
                <thead><tr><th>Measure</th><th>Value</th></tr></thead>
                <tbody>
                <tr><th>MinKNOW version </th><td> {0} </td></tr>
                <tr><th>Basecaller name</th><td> {1} </td></tr>
                <tr><th>Basecaller version</th><td> {2} </td></tr>
                <tr><th>Basecaller analysis</th><td> {3} </td></tr>
                <tr><th>ToulligQC version</th><td> {4} </td></tr>
                <tr><th>Hostname</th><td> {5} </td></tr>
                <tr><th>Device</th><td> {6} </td></tr>
                <tr><th>Device ID</th><td> {7} </td></tr>
                <tr><th>Model file</th><td> {8} </td></tr>
                </tbody>
            </table>
        </div> <!-- end .info-box-left -->
        <div class=\"after-box\"><p></p></div>
      </div>
    """.format(minknow_version,basecaller_name, basecaller_version, basecaller_analysis, config_dictionary['app.version'],hostname,device_type,device_id,model_file)

    for i, t in enumerate(graphs):
      if len(t)==5:
        main_report += "      <div class=\"module\" id=M{0}></div>".format(i)
        main_report += t[4]
        if t[2] is None:
          main_report += "      <div class=\"after-box\"></div>\n"
        else:
          main_report += "      <div class=\"box-left\">\n {} </div>\n".format(t[2])
          main_report += "      <div class=\"after-box\"><p></p></div>\n"

        
      else:
        main_report += "      <div class=\"module\" id=M{0}><h2> {1} <a title=\"<b>{4}</b>\">&#x1F263;</a></h2></div>" \
            .format(i, t[0], _embedded_image(t[1]), t[2], t[3])

        if t[2] is None:
            main_report += "      <div class=\"module\"><p><img src=\"{2}\" " \
                           "alt=\"{1} image\"></p></div>\n".format(i, t[0], _embedded_image(t[1]))
        else:
            main_report += "      <div class=\"box\"><img src=\"{2}\">" \
                           "</div>\n".format(i, t[0], _embedded_image(t[1]), t[2])

            main_report += "      <div class=\"box-left\">\n {3}" \
                           "</div>\n".format(i, t[0], _embedded_image(t[1]), t[2])

            main_report += "      <div class=\"after-box\"><p></p></div>\n"
    main_report += "    </div>\n"

    # Add all the element of the page
    report = title + header + banner + summary + main_report + footer

    # Write the HTML page
    f.write(report)
    f.close()


def _embedded_image(image_path):
    """
    Embedded an image
    :param image_path: path of the image
    :return: a string with the image in base64
    """
    with open(image_path, "rb") as image_file:
        result = "data:image/png;base64," + base64.b64encode(image_file.read()).decode('ascii')

    return result


def _get_result_value(result_dict, key , default_value = ""):
    """
    Get the value of the result dictionary or a default value if the key does not exists.
    :param result_dict: result dictionary
    :param key: the key to use
    :param default_value: the default value
    :return: the value of key in the dictionary or the default value if the key does not exists in the dictionary
    """
    if key in result_dict:
        return result_dict[key]
    else:
        return default_value
