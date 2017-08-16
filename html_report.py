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
    f = open(result_directory+'report.html', 'w')
    read_count = image_directory+"read_count_histogram.png"
    phred_score = image_directory+"read_quality_boxplot.png"
    channel_count = image_directory+"channel_count_histogram.png"
    read_number = image_directory+"read_number_run.png"
    read_length = image_directory+"read_length_histogram.png"
    channel_occupancy = image_directory+"channel_occupancy.png"
    barcode_pie_chart = image_directory+"barcode_percentage_pie_chart.png"
    barcode_length_boxplot = image_directory+'barcode_length_boxplot.png'
    scatterplot = image_directory+"scatter_plot.png"
    phred_score_frequency = image_directory+"phred_score_frequency.png"
    barcode_phred_score_boxplot = image_directory+"barcode_phred_score_boxplot.png"
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
            <h1> RUN MinION report </h1>
            <p>Flowcell id:{0}</p>
            <p>Date:{1}</p>
            <p>Number of reads:{2}</p>
            <h3>Histogram of read count</h3>
            <img src={3} alt=read_count,width=500, height=400>
            <h3>Histogram of read length</h3>
            <img src={4} alt=read_count,width=1000, height=700>
            <h3>Phred score according to the read type</h3>
            <img src="{5}" alt=phred_score, width=800, height=700>
            <h3>Channel counts</h3>
            <img src="{6}" alt=channel_count,width=600, height=400>
            <h3>Curve representing the reads number produced against the time</h3>
            <img src="{7}" alt=read_number, width=700, height=400>
            <h3>Channel occupancy</h3>
            <img src="{8}" alt=channel_occupancy, width=1000,height=700>
            <h3>Phred score frequency</h3>
            <img src="{9}" alt=phred score frequency, width=800, height=700>
            <h3>Relation between the sequence length template and the mean qscore template</h3>
            <img src="{10}" alt=scatter plot, width=600, height=500>
            """.format(flowcell_id, run_date, number_of_read, read_count, read_length, phred_score, channel_count, read_number, channel_occupancy, phred_score_frequency, scatterplot)
    if is_barcode:
        message2 = """<h3>Barcode pie chart</h3>
            <img src="{0}" alt="Barcode pie chart", width=600, height=400>
            <h3>Boxplot of read length distribution for each barcode</h3>
            <img src="{1}" alt="Barcode read length boxplot", width=800, height=700>
            <h3>Boxplot of phred score distribution by barcode</h3>
            <img src="{2}" alt="Barcode read length boxplot", width=800, height=700>
        </body>
        </html>
        """.format(barcode_pie_chart, barcode_length_boxplot, barcode_phred_score_boxplot)

    else:
        message2 = """</body>
        </html>
        """

    f.write(message1)
    f.write(message2)
    f.close()
