def html_report(result_directory, run_date, flowcell_id, is_barcode):
    image_directory = result_directory + 'images/'
    f = open(result_directory+'report.html', 'w')
    read_count = image_directory+"read_count_histogram.png"
    phred_score = image_directory+"read_quality_boxplot.png"
    channel_count = image_directory+"channel_count_histogram.png"
    read_number = image_directory+"read_number_run.png"
    read_length = image_directory+"read_length_histogram.png"
    channel_occupancy = image_directory+"channel_occupancy.png"
    barcode_pie_chart = image_directory+"barcode_percentage_pie_chart.png"
    barcode_total = image_directory+"barcode_total.png"
    scatterplot = image_directory+"scatter_plot.png"
    phred_score_frequency = image_directory+"phred_score_frequency.png"

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
            <h3>Histogram of read count</h3>
            <img src={2} alt=read_count,width=600, height=300>
            <h3>Histogram of read length</h3>
            <img src={3} alt=read_count,width=1000, height=700>
            <h3>Phred score according to the read type</h3>
            <img src="{4}" alt=phred_score, width=600, height=300>
            <h3>Channel counts</h3>
            <img src="{5}" alt=channel_count,width=600, height=300>
            <h3>Curve representing the reads number produced against the time</h3>
            <img src="{6}" alt=read_number, width=500, height=300>
            <h3>Channel occupancy</h3>
            <img src="{7}" alt=channel_occupancy, width=1000,height=700>
            <h3>Phred score frequency</h3>
            <img src="{8}" alt=phred score frequency, width=600, height=300>
            <h3>Relation between the sequence length template and the mean qscore template</h3>
            <img src="{9}" alt=scatter plot, width=600, height=300>
            """.format(flowcell_id, run_date, read_count, read_length, phred_score, channel_count, read_number, channel_occupancy, phred_score_frequency, scatterplot)
    if is_barcode:
        message2 = """<h3>Barcode pie chart</h3>
            <img src="{}" alt="Barcode pie chart", width=600, height=400>
            <h3>Barcode read length boxplot</h3>
            <img src="{}" alt="Barcode read length boxplot", width=600, height=300>
        </body>
        </html>
        """.format(barcode_pie_chart,barcode_total)
    else:
        message2 = """</body>
        </html>
        """

    f.write(message1)
    f.write(message2)
    f.close()
