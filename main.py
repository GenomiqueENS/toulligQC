import matplotlib
matplotlib.use('Agg')
import basecalling_stat_plotter1D
from matplotlib.backends.backend_pdf import PdfPages
import sys
from PyPDF2 import PdfFileMerger
import pandas as pd
import fast5_data_extractor
import docxs
import configparser

pdf = PdfPages('Rapport_pdf.pdf')

run_name = sys.argv[1]


configParser = configparser.ConfigParser()
configFilePath = r'config2.txt'
configParser.read(configFilePath)
#configParser.get('ferrato-config', 'fast5.directory')+'raw/'+run_name+'/0'
fast5_directory = input('path to fast5 files:')
basecall_log =configParser.get('config', 'log.file')+run_name+'/sequencing_summary.txt'


fast5_data = fast5_data_extractor.fast5_data_extractor(fast5_directory)
basecalling = basecalling_stat_plotter1D.basecalling_stat_plotter1D(basecall_log,pdf, run_name)


#Date and flowcell id

flowcell_id, *_ = fast5_data

#Histogram of read counts according to the type read

basecalling.read_count_histogram()

#Phred score according to the read type
basecalling.read_quality_boxplot()


#Channel counts
basecalling.channel_count_histogram()

#Pie chart representing barcodes
basecalling.barcode_percentage_pie_chart()

#Curve representing the number of reads produced along the runtime
basecalling.read_number_run()

#Representation of the channels occupation.
#The frame represents the flowcell containing 512 channels
#The layout is present in the pdf file result.pdf
channel_count = basecalling.channel
total_number_reads_per_pore = pd.value_counts(channel_count)
basecalling.plot_performance(total_number_reads_per_pore)
basecalling.occupancy_pore()
basecalling.reads_size_selection_barcode()
pdf.close()



input1 = open("Rapport_pdf.pdf", "rb")
input2 = open("layout.pdf", "rb")

pdfs = ["Rapport_pdf.pdf","layout.pdf"]

merger = PdfFileMerger()

for pdf in pdfs:
    merger.append(open(pdf, 'rb'))

with open('result.pdf', 'wb') as fout:
    merger.write(fout)


docxs.docxs(basecalling.selection,basecalling.date(), flowcell_id)
basecalling.statistics_dataframe()


#import log_file1D
#log_file1D.log_file1D(fast5_data, basecalling)



#/home/ferrato/shares-net/sequencages/nanopore/albacore-logs/FAF04250_20170328/sequencing_summary.txt
#/home/ferrato/shares-net/sequencages/nanopore/test_alabacore/save/workspace

#/home/ferrato/Documents/text.txt
#/home/ferrato/shares-net/sequencages/nanopore/fastq/FAF04250_20170328
#/home/ferrato/ownCloud/Documents/fast5_1D/design.csv
