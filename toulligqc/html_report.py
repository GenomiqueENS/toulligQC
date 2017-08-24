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
import re

def html_report(config_dictionary, result_dict):
    '''
    Creation of a html report
    :param config_dictionary: dictionary containing file or directory paths
    :param result_dict: result dictionary containing all statistics
    '''
    result_directory = config_dictionary['result_directory']
    sequence_length_template = result_dict['sequence_length_template']
    is_barcode = config_dictionary['barcoding']
    run_name = config_dictionary['run_name']
    flowcell_id = result_dict['flowcell_id']
    run_date = result_dict['run_date']
    date = re.search('(\d{4})(\d{2})(\d{2})', run_date).groups()
    run_date = '{}-{}-{}'.format(date[2], date[1], date[0])
    image_directory = result_directory + 'images/'
    f = open(result_directory + 'report.html', 'w')
    read_count = image_directory + "read_count_histogram.png"
    phred_score = image_directory + "read_quality_boxplot.png"
    channel_count = image_directory + "channel_count_histogram.png"
    read_number = image_directory + "read_number_run.png"
    read_length = image_directory + "read_length_histogram.png"
    channel_occupancy = image_directory + "channel_occupancy.png"
    barcode_pie_chart = image_directory + "barcode_percentage_pie_chart.png"
    barcode_length_boxplot = image_directory + 'barcode_length_boxplot.png'
    scatterplot = image_directory + "scatter_plot.png"
    phred_score_frequency = image_directory + "phred_score_frequency.png"
    barcode_phred_score_boxplot = image_directory + "barcode_phred_score_boxplot.png"

    number_of_read = len(sequence_length_template)
    if is_barcode:
        report = """<!DOCTYPE html>
        <html>
        <head>
            <title>Rapport run MinION</title>
            <meta charset='UTF-8'>
            <style>
            </style>
        </head>

        <body>
            <style>
            @media screen {{
                div.summary {{
                width: 16em;
                position:fixed;
                top: 3em;
                margin:1em 0 0 1em;
            }}

      div.main {{
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
      }}

      div.header {{
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
      }}

      div.footer {{
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
      }}

      img.indented {{
        margin-left: 3em;
      }}
     }}

     @media print {{
        img {{
            max-width:100% !important;
            page-break-inside: avoid;
        }}
        h2, h3 {{
            page-break-after: avoid;
        }}
        div.header {{
          background-color: #FFF;
        }}

      }}
                body {{
        font-family: sans-serif;
        color: #000;
        background-color: #FFF;
        border: 0;
        margin: 0;
        padding: 0;
        }}

        div.header {{
        border:0;
        margin:0;
        padding: 0.5em;
        font-size: 200%;
        font-weight: bold;
        width:100%;
        }}

        #header_title {{
        display:inline-block;
        float:left;
        clear:left;
        }}
        #header_filename {{
        display:inline-block;
        float:right;
        clear:right;
        font-size: 50%;
        margin-right:2em;
        text-align: right;
        }}

        div.header h3 {{
        font-size: 50%;
        margin-bottom: 0;
        }}

        div.summary ul {{
        padding-left:0;
        list-style-type:none;
        }}

        div.summary ul li img {{
        margin-bottom:-0.5em;
        margin-top:0.5em;
        }}

        div.main {{
        background-color: white;
        }}

        div.module {{
        padding-bottom:1.5em;
        padding-top:1.5em;
        }}

        div.footer {{
        background-color: #EEE;
        border:0;
        margin:0;
        padding: 0.5em;
        font-size: 100%;
        font-weight: bold;
        width:100%;
        }}


        a {{
        color: #000080;
        }}

        a:hover {{
        color: #800000;
        }}

        h2 {{
        color: #800000;
        padding-bottom: 0;
        margin-bottom: 0;
        clear:left;
        }}

        table {{
        margin-left: 3em;
        text-align: center;
        }}

        th {{
        text-align: center;
        background-color: #000080;
        color: #FFF;
        padding: 0.4em;
        }}

        td {{
        font-family: monospace;
        text-align: left;
        background-color: #EEE;
        color: #000;
        padding: 0.4em;
        }}

        img {{
        padding-top: 0;
        margin-top: 0;
        border-top: 0;
        }}


        p {{
        padding-top: 0;
        margin-top: 0;
        }}

        </style>

        <div class='summary'>
            <h2>Summary</h2>
            <ol>
              <li>
                <a href="#M0">Histogram of read count</a>
              </li>
              <li>
                <a href="#M1">Histogram of read length</a>
              </li>
              <li>
                <a href="#M2">Phred score according to the read type </a>
              </li>
              <li>
                <a href="#M3">Channel counts</a>
              </li>
              <li>
                <a href="#M4">Curve representing the reads number produced against the time</a>
              </li>
              <li>
                <a href="#M5">Channel occupancy</a>
              </li>
              <li>
                <a href="#M6">Phred score frequency</a>
              </li>
              <li>
                <a href="#M7">Relation between the sequence length template and the mean qscore template</a>
              </li>
              <li>
                <a href="#M8">Barcode pie chart</a>
              </li>
              <li>
                <a href="#M9">Boxplot of read length distribution for each barcode</a>
              </li>
              <li>
                <a href="#M10">Boxplot of phred score distribution by barcode</a>
              </li>
            </ol>
          </div>
          <div class="header">
            <div id="header_title">
              Run MinION report<br>

            </div>
            <div id="header_filename">
              Run name: {0}<br>
              Run date: {2}<br>
              Flowcell id: {1}
            </div>
          </div>

          <div class = 'main'>

            <div class="module">
            <p><b>Number of reads: {3}</b></p>


            <h2 id=M0>
              Histogram of read count
            </h2>
            <p>
               <img src={4} alt=read_count>
           </div>
          <div class="module">
            <h2 id=M1>
              Histogram of read length
            </h2>
            <p>
             <img src={5} alt=read_count>
            </p>
          </div>
          <div class="module">
            <h2 id=M2>
              Phred score according to the read type
            </h2>
            <p>
                <img src="{6}" alt=phred_score>
            </p>
          </div>

          <div class="module">
            <h2 id=M3>Channel counts</h2>
            <p>
                <img src="{7}" alt=channel_count>
            </p>
          </div>
            <div class="module">
              <h2 id=M4>Curve representing the reads number produced against the time</h2>
              <p>
                    <img src="{8}" alt=read_number, width=700, height=400>
              </p>
            </div>
            <div class="module">
              <h2 id="M5">Channel occupancy</h2>
              <p>
                <img src="{9}" alt=channel_occupancy>
              </p>
            </div>
            <div class="module">
              <h2 id="M6">Phred score frequency</h2>
              <p>
                <img src="{10}" alt=phred score frequency>
              </p>
            </div>

            <div class="module">
              <h2 id="M7">Relation between the sequence length template and the mean qscore template</h2>
              <p>
                <img src="{11}" alt=scatter plot>
              </p>
            </div>
     <div class="module">
              <h2 id="M8">Barcode pie chart</h2>
              <p>
                <img src="{12}" alt="Barcode pie chart", width=600, height=400>
              </p>
            </div>

            <div class="module">
              <h2 id="M9">Boxplot of read length distribution for each barcode</h2>
              <img src="{13}" alt="Barcode read length boxplot", width=800, height=700>
            </div>

            <div class="module">
                <h2 id="M10">Boxplot of phred score distribution by barcode</h2>
                <p>
                <img src="{14}" alt="Barcode read length boxplot", width=800, height=700>
                </p>
            </div>

        </div>
         <div class="footer">
            Produced by <a href="https://github.com/GenomicParisCentre/toulligQC">ToulligQC</a> (version 0.0.1)
            </div>
        </body>
        </html>""".format(run_name, flowcell_id, run_date, number_of_read,read_count, read_length, phred_score, \
                          channel_count, read_number, channel_occupancy, phred_score_frequency, scatterplot, barcode_pie_chart, \
                          barcode_length_boxplot, barcode_phred_score_boxplot)


    else:
        report = """
                 <style>
            @media screen {{
                div.summary {{
                width: 16em;
                position:fixed;
                top: 3em;
                margin:1em 0 0 1em;
            }}

  div.main {{
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
  }}

  div.header {{
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
  }}

  div.footer {{
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
  }}

  img.indented {{
    margin-left: 3em;
  }}
 }}

 @media print {{
	img {{
		max-width:100% !important;
		page-break-inside: avoid;
	}}
	h2, h3 {{
		page-break-after: avoid;
	}}
	div.header {{
      background-color: #FFF;
    }}

  }}
            body {{
    font-family: sans-serif;
    color: #000;
    background-color: #FFF;
    border: 0;
    margin: 0;
    padding: 0;
    }}

    div.header {{
    border:0;
    margin:0;
    padding: 0.5em;
    font-size: 200%;
    font-weight: bold;
    width:100%;
    }}

    #header_title {{
    display:inline-block;
    float:left;
    clear:left;
    }}
    #header_filename {{
    display:inline-block;
    float:right;
    clear:right;
    font-size: 50%;
    margin-right:2em;
    text-align: right;
    }}

    div.header h3 {{
    font-size: 50%;
    margin-bottom: 0;
    }}

    div.summary ul {{
    padding-left:0;
    list-style-type:none;
    }}

    div.summary ul li img {{
    margin-bottom:-0.5em;
    margin-top:0.5em;
    }}

    div.main {{
    background-color: white;
    }}

    div.module {{
    padding-bottom:1.5em;
    padding-top:1.5em;
    }}

    div.footer {{
    background-color: #EEE;
    border:0;
    margin:0;
    padding: 0.5em;
    font-size: 100%;
    font-weight: bold;
    width:100%;
    }}


    a {{
    color: #000080;
    }}

    a:hover {{
    color: #800000;
    }}

    h2 {{
    color: #800000;
    padding-bottom: 0;
    margin-bottom: 0;
    clear:left;
    }}

    table {{
    margin-left: 3em;
    text-align: center;
    }}

    th {{
    text-align: center;
    background-color: #000080;
    color: #FFF;
    padding: 0.4em;
    }}

    td {{
    font-family: monospace;
    text-align: left;
    background-color: #EEE;
    color: #000;
    padding: 0.4em;
    }}

    img {{
    padding-top: 0;
    margin-top: 0;
    border-top: 0;
    }}


    p {{
    padding-top: 0;
    margin-top: 0;
    }}

    </style>

    <div class='summary'>
        <h2>Summary</h2>
        <ol>
          <li>
            <a href="#M0">Histogram of read count</a>
          </li>
          <li>
            <a href="#M1">Histogram of read length</a>
          </li>
          <li>
            <a href="#M2">Phred score according to the read type </a>
          </li>
          <li>
            <a href="#M3">Channel counts</a>
          </li>
          <li>
            <a href="#M4">Curve representing the reads number produced against the time</a>
          </li>
          <li>
            <a href="#M5">Channel occupancy</a>
          </li>
          <li>
            <a href="#M6">Phred score frequency</a>
          </li>
          <li>
            <a href="#M7">Relation between the sequence length template and the mean qscore template</a>
          </li>
          <li>
            <a href="#M8">Barcode pie chart</a>
          </li>
          <li>
            <a href="#M9">Boxplot of read length distribution for each barcode</a>
          </li>
          <li>
            <a href="#M10">Boxplot of phred score distribution by barcode</a>
          </li>
        </ol>
      </div>
      <div class="header">
        <div id="header_title">
          Run MionION report
        </div>
        <div id="header_filename">
          Run name: {0},br>
          Run date: {1}<br>
          Flowcell id: {2}
        </div>
      </div>

      <div class = 'main'>


        <div class="module">
        <p><b>Number of reads: {3}</b></p>

        <h2 id=M0>
          Histogram of read count
        </h2>
        <p>
           <img src={4} alt=read_count>
       </div>
      <div class="module">
        <h2 id=M1>
          Histogram of read length
        </h2>
        <p>
         <img src={5} alt=read_length>
        </p>
      </div>
      <div class="module">
        <h2 id=M2>
          Phred score according to the read type
        </h2>
        <p>
            <img src="{6}" alt=phred_score>
        </p>
      </div>

      <div class="module">
        <h2 id=M3>Channel counts</h2>
        <p>
            <img src="{7}" alt=channel_count>
        </p>
      </div>
        <div class="module">
          <h2 id=M4>Curve representing the reads number produced against the time</h2>
          <p>
                <img src="{8}" alt=read_number, width=700, height=400>
          </p>
        </div>
        <div class="module">
          <h2 id="M5">Channel occupancy</h2>
          <p>
            <img src="{9}" alt=channel_occupancy>
          </p>
        </div>
        <div class="module">
          <h2 id="M6">Phred score frequency</h2>
          <p>
            <img src="{10}" alt=phred score frequency>
          </p>
        </div>

        <div class="module">
          <h2 id="M7">Relation between the sequence length template and the mean qscore template</h2>
          <p>
            <img src="{11}" alt=scatter plot>
          </p>
        </div>
        <div class="footer">
      Produced by <a href="https://github.com/GenomicParisCentre/toulligQC">ToulligQC</a> (version 0.0.1)
         </div>
    </div>
    </body>

    </html>""".format(run_name, run_date, flowcell_id, number_of_read, read_count, read_length, phred_score,
                      channel_count, read_number, channel_occupancy, phred_score_frequency, scatterplot)


    f.write(report)
    f.close()
