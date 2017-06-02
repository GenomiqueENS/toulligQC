import matplotlib
matplotlib.use('Agg')
import basecalling_stat_plotter1D
from matplotlib.backends.backend_pdf import PdfPages
import sys
from PyPDF2 import PdfFileMerger
import pandas as pd
import fast5_data_extractor
import docxs
import log_file1D
import os
import configparser

run_name = sys.argv[1]

configParser = configparser.ConfigParser()

#configParser.get('ferrato-config', 'fast5.directory')+'raw/'+run_name+'/0'
bz2_file_path = input('Path to bz2 fast5 files:')
barcode_present = input('Did you use barcodes ? Answer by y(yes) or n(no):')
question = input('Must the analysis performed on specific bz2 file ? Answer by y(yes) or n(no):')
if question == 'y':
    file_list = input('Enter your file (or file list) separated by a space:')
    file_list = file_list.split(" ")
else:
    file_list = 'None'
try:
    #In the docker image
    configFilePath = r'/configpass/docker_config.txt'
    basecall_log = '/log.file/' +run_name+'/sequencing_summary.txt'
    report_writing_directory = '/design.file.directory/'

except:
    configFilePath = r'config.txt'
    configParser.read(configFilePath)
    basecall_log = configParser.get('config', 'log.file') + run_name + '/sequencing_summary.txt'
    report_writing_directory = configParser.get('config', 'design.file.directory')
    
pdf_report = report_writing_directory+'Rapport_pdf.pdf'
pdf = PdfPages(pdf_report)

fast5_data = fast5_data_extractor.fast5_data_extractor(bz2_file_path)
basecalling = basecalling_stat_plotter1D.basecalling_stat_plotter1D(basecall_log,pdf, run_name, barcode_present, file_list)

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

if barcode_present == 'y':
    # Pie chart representing barcodes
    basecalling.barcode_percentage_pie_chart()

basecalling.read_length_histogram()

#Representation of the channels occupation.
#The frame represents the flowcell containing 512 channels
#The layout is present in the pdf file result.pdf
channel_count = basecalling.channel
total_number_reads_per_pore = pd.value_counts(channel_count)
basecalling.plot_performance(total_number_reads_per_pore)
basecalling.occupancy_pore()

pdf.close()


input1 = open(report_writing_directory, "rb")
input2 = open("layout.pdf", "rb")

pdfs = ["Rapport_pdf.pdf","layout.pdf"]

merger = PdfFileMerger()

for pdf in pdfs:
    merger.append(open(pdf, 'rb'))
    
result_pdf_path =os.path.join(report_writing_directory, 'result.pdf')

with open(report_writing_directory, 'wb') as fout:
    merger.write(fout)

if barcode_present == 'y':
    docxs.docxs(basecalling.selection,basecalling.run_date(), flowcell_id, barcode_present)
    basecalling.statistics_dataframe()
else:
    docxs.docxs('', basecalling.run_date(), flowcell_id, barcode_present)

if barcode_present == 'n':
    log_file1D.log_file1D(fast5_data, basecalling)
