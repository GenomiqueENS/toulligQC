from PyPDF2 import PdfFileMerger
import matplotlib
matplotlib.use('Agg')
import basecalling_stat_plotter1D
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import fast5_data_extractor
import docxs
import log_file1D
import os
import parser

run_name, selected_file, is_docker, is_barcode = parser.get_args()
 
dico_path = parser.file_path_initialization()
pdf_report = dico_path['result_directory']

fast5_directory = input('path to fast5 files:')
basecall_log = dico_path['basecall_log'] +run_name+'/sequencing_summary.txt'

#configParser.get('ferrato-config', 'fast5.directory')+'raw/'+run_name+'/0'
bz2_file_path = input('Path to bz2 fast5 files:')
###########barcode_present = input('Did you use barcodes ? Answer by y(yes) or n(no):')
###########question = input('Must the analysis performed on specific bz2 file ? Answer by y(yes) or n(no):')

report_writing_directory = dico_path['result_directory']
pdf_report = report_writing_directory+'Rapport_pdf.pdf'
pdf = PdfPages(pdf_report)

fast5_data = fast5_data_extractor.fast5_data_extractor(bz2_file_path)
basecalling = basecalling_stat_plotter1D.basecalling_stat_plotter1D(basecall_log, pdf, is_barcode, selected_file)

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

basecalling.read_length_histogram()

#Representation of the channels occupation.
#The frame represents the flowcell containing 512 channels
#The layout is present in the pdf file result.pdf
channel_count = basecalling.channel
total_number_reads_per_pore = pd.value_counts(channel_count)
basecalling.plot_performance(total_number_reads_per_pore)
basecalling.occupancy_pore()

pdf.close()


report_pdf_file = os.path.join(report_writing_directory, 'Rapport_pdf.pdf')

if is_docker:
    pdfs = [report_pdf_file,"/scripts/toulligQC/layout.pdf"]
else:
    pdfs = [report_pdf_file,"layout.pdf"]

merger = PdfFileMerger()

for pdf in pdfs:
    merger.append(open(pdf, 'rb'))
    
result_pdf_path =os.path.join(report_writing_directory, 'result.pdf')    

with open(result_pdf_path, 'wb') as fout:
    merger.write(fout)



if is_barcode:
    docxs.docxs(basecalling.barcode_selection,basecalling.run_date(), flowcell_id, is_barcode)
    basecalling.statistics_dataframe()
else:
    log_file1D.log_file1D(fast5_data, basecalling)
    docxs.docxs('', basecalling.run_date(), flowcell_id, is_barcode)



