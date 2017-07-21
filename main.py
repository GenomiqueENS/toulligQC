import matplotlib
matplotlib.use('Agg')
import basecalling_stat_plotter1D
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import fast5_data_extractor
#import docxs
import log_file1D
import os
import parser
import html_report
import shutil
import sys

run_name, config_file, is_barcode, file_list = parser.get_args()
dico_path = parser.config_file_initialization(config_file, is_barcode, run_name, file_list)
if not dico_path:
    sys.exit("Error, dico_path is empty")

result_directory        = dico_path['result_directory']
basecall_log            = dico_path['basecall_log'] 
fastq_directory         = dico_path['fastq_directory'] 
fast5_directory         = dico_path['fast5_directory']

print(result_directory)

dico_extension = parser.extension(config_file, file_list)
fast5_file_extension = dico_extension['fast5_file_extension']
fastq_file_extension = dico_extension['fastq_file_extension']

if is_barcode:
    design_file_directory = dico_path['design_file_directory']
else:
    design_file_directory = ''

pdf_report = dico_path['result_directory']
basecall_log = dico_path['basecall_log'] +run_name+'/sequencing_summary.txt'

if os.path.isdir(result_directory):
    shutil.rmtree(result_directory, ignore_errors=True)
    os.makedirs(result_directory)
else:
    os.makedirs(result_directory)
pdf_report = result_directory+'Rapport_pdf.pdf'
pdf = PdfPages(pdf_report)


basecalling = basecalling_stat_plotter1D.basecalling_stat_plotter1D(basecall_log, pdf, is_barcode,result_directory, fastq_directory, dico_extension, design_file_directory,  file_list)
fast5_data = fast5_data_extractor.fast5_data_extractor(fast5_directory, result_directory, dico_extension)

#Date and flowcell id
flowcell_id, *_ = fast5_data

#Histogram of read counts according to the type read
basecalling.read_count_histogram()

#Phred score according to the read type
basecalling.read_quality_boxplot()

#Channel counts
basecalling.channel_count_histogram()

#Curve representing the number of reads produced along the runtime
basecalling.read_number_run()

if is_barcode:
    # Pie chart representing barcodes
    basecalling.barcode_percentage_pie_chart()
    basecalling.barcode_length_boxplot()
    basecalling.barcoded_phred_score_frequency()
basecalling.read_length_histogram()

#Representation of the channels occupation.
#The frame represents the flowcell containing 512 channels
#The layout is present in the pdf file result.pdf
channel_count = basecalling.channel
total_number_reads_per_pore = pd.value_counts(channel_count)
basecalling.plot_performance(total_number_reads_per_pore)
basecalling.occupancy_pore()

basecalling.phred_score_frequency()
basecalling.scatterplot()


pdf.close()

report_pdf_file = os.path.join(result_directory, 'Rapport_pdf.pdf')
html_report.html_report(result_directory, basecalling.run_date(), flowcell_id, is_barcode)

if is_barcode:
    #docxs.docxs(basecalling.barcode_selection,basecalling.run_date(), flowcell_id, is_barcode, result_directory, design_file_directory)
    basecalling.statistics_dataframe()
else:
    log_file1D.log_file1D(fast5_data, basecalling, result_directory)
    #docxs.docxs('', basecalling.run_date(), flowcell_id, is_barcode, result_directory)



