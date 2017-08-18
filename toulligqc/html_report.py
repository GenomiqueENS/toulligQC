def html_report(result_directory, run_date, flowcell_id, is_barcode, sequence_length_template):
    '''
    Creation of a html report
    :param result_directory: result directory
    :param run_date: run date
    :param flowcell_id: flowcell id
    :param is_barcode: boolean indicated if we use the barcodes
    :param sequence_length_template: sequence length template
    '''
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
    message1 = """
    <!DOCTYPE html>
    <html>
    <head>
            <title>Rapport run MinION</title>
            <meta charset='UTF-8'>
            <style>
            </style>
    </head>

    <body>
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
          Thu 20 Apr 2017
        </div>
      </div>

      <div class = 'main'>
        <div class="module">


        <h2 id=M0>
          Histogram of read count
        </h2>
        <p>
           <img src={3} alt=read_count>
       </div>
      <div class="module">
        <h2 id=M1>
          Histogram of read length
        </h2>
        <p>
         <img src={4} alt=read_count>
        </p>
      </div>
      <div class="module">
        <h2 id=M2>
          Phred score according to the read type
        </h2>
        <p>
            <img src="{5}" alt=phred_score>
        </p>
      </div>

      <div class="module">
        <h2 id=M3>Channel counts</h2>
        <p>
            <img src="{6}" alt=channel_count>
        </p>
      </div>
        <div class="module">
          <h2 id=M4>Curve representing the reads number produced against the time</h2>
          <p>
                <img src="{7}" alt=read_number, width=700, height=400>
          </p>
        </div>
        <div class="module">
          <h2 id="M5">Channel occupancy</h2>
          <p>
            <img src="{8}" alt=channel_occupancy>
          </p>
        </div>
        <div class="module">
          <h2 id="M6">Phred score frequency</h2>
          <p>
            <img src="{9}" alt=phred score frequency>
          </p>
        </div>

        <div class="module">
          <h2 id="M7">Relation between the sequence length template and the mean qscore template</h2>
          <p>
            <img src="{10}" alt=scatter plot>
          </p>
        </div>

         <div class="footer">
      Produced by <a href="https://github.com/GenomicParisCentre/toulligQC">ToulligQC</a> (version 0.0.1)
         </div>
      </div>""".format(flowcell_id, run_date, number_of_read, read_count, read_length, phred_score, channel_count, \
                       read_number, channel_occupancy, phred_score_frequency, scatterplot)

    if is_barcode:
        message2 = """
        <div class="module">
          <h2 id="M8">Barcode pie chart</h2>
          <p>
            <img src="/Users/lionelferrato/result_directory/FAF04250/images/barcode_percentage_pie_chart.png" alt="Barcode pie chart", width=600, height=400>
          </p>
        </div>

        <div class="module">
          <h2 id="M9">Boxplot of read length distribution for each barcode</h2>
          <img src="/Users/lionelferrato/result_directory/FAF04250/images/barcode_length_boxplot.png" alt="Barcode read length boxplot", width=800, height=700>
        </div>

        <div class="module">
            <h2 id="M10">Boxplot of phred score distribution by barcode</h2>
            <p>
            <img src="/Users/lionelferrato/result_directory/FAF04250/images/barcode_phred_score_boxplot.png" alt="Barcode read length boxplot", width=800, height=700>
            </p>
        </div>
        </div>
        <div class="footer">
        Produced by <a href="https://github.com/GenomicParisCentre/toulligQC">ToulligQC</a> (version 0.0.1)
        </div>
    </body>
    </html>""".format(barcode_pie_chart, barcode_length_boxplot, barcode_phred_score_boxplot)
    else:
        message2 = """</body>
                    </html>
                    """

    f.write(message1)
    f.write(message2)
    f.close()
