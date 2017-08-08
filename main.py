import matplotlib
matplotlib.use('Agg')
import basecalling_stat_plotter1D
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import fast5_data_extractor
import log_file1D
import os
import html_report
import shutil
import sys
import csv
import re
import configparser
import argparse
from pathlib import Path

def get_args():
    '''
    Parsing the command line
    :return: different informations: run name,
    path towards files contained in the configuration file,
    boolean indicating if we use the barcode,
    path list used with the -f option
    '''
    parser = argparse.ArgumentParser()
    home_path = str(Path.home())
    parser.add_argument("-n", "--run_name", action='store', dest="run_name", help="Run name",default=True)
    parser.add_argument("-b","--barcode", action='store_true',dest='is_barcode',help="Barcode usage",default=False)
    parser.add_argument("-c", "--config_file", action='store', dest='config_file', help="Path to the configuration file", default=home_path+'/.toulligqc/config.txt')
    parser.add_argument('-f', '--arg', nargs='+', type=str, help="Path list of different files or directory needed when \
    we have not a config file. The paths indicated in the list must be in the same order that in the configuration file")

    argument_list = vars(parser.parse_args())
    file_list = argument_list["arg"]
    argument_value=parser.parse_args()
    run_name = argument_value.run_name
    config_file_path = argument_value.config_file
    is_barcode = argument_value.is_barcode

    print(config_file_path)
    return run_name, config_file_path, is_barcode, file_list


def config_file_initialization(config_file, is_barcode, run_name, list_file=''):
    """
    Creation of a dictionary of the path file contained in the configuration file

    :param config_file: configuration file path
    :param is_barcode: boolean
    :param run_name: Run name
    :param list_file: path list used with the -f option
    :return: Dictionary containing the file paths included into the configuration file
    """
    dico_path = {}
    configFilePath = config_file
    config = configparser.ConfigParser()
    config.read(configFilePath)
    print(config.get('config', 'result.directory'))

    if list_file:
        dico_path = {}
        dico_path['fast5_directory'] = list_file[0]
        dico_path['basecall_log'] = list_file[1]
        dico_path['fastq_directory'] = list_file[2]
        dico_path['result_directory'] = list_file[3]

        if is_barcode:
            dico_path['design_file_directory'] = list_file[4]

    elif config_file:
        dico_path['result_directory'] = config.get('config', 'result.directory')
        dico_path['basecall_log'] = config.get('config', 'log.file')
        dico_path['fastq_directory'] = config.get('config', 'fastq.directory')
        dico_path['fast5_directory'] = config.get('config', 'fast5.directory')

        if is_barcode:
            dico_path['design_file_directory'] = config.get('config', 'design.file.directory')

    else:
        print('Error, not a config file')
        return 0

    for key, value in dico_path.items():
        if value.endswith('/'):
            continue
        else:
            dico_path[key] = value + '/'

    dico_path['result_directory'] = dico_path['result_directory'] + run_name + '/'

    return dico_path


def extension(config_file, list_file, is_barcode):
    '''
    Creation of a dictionary containing the extension used for the fast5 and fastq files
    :param config_file: configuration file
    :param list_file: path list used with the -f option
    :return: Extension dictionary used for the fast5 and fastq files
    '''
    dico_extension = {}
    if list_file:
        if is_barcode:
            dico_extension['fast5_file_extension'] = list_file[5]
            dico_extension['fastq_file_extension'] = list_file[6]
        else:
            dico_extension['fast5_file_extension'] = list_file[4]
            dico_extension['fastq_file_extension'] = list_file[5]
    else:
        config = configparser.ConfigParser()
        config.read(config_file)
        dico_extension['fast5_file_extension'] = config.get('extension', 'fast5.file.extension')
        dico_extension['fastq_file_extension'] = config.get('extension', 'fastq.file.extension')

    return dico_extension

def graph_creation(basecalling, is_barcode):
    '''
    Creation of the different graphs produced by toulligQC
    :param basecalling: basecalling_stat_plotter instance
    :param is_barcode: boolean indicating if we use the barcodes or not
    '''

    basecalling.read_count_histogram()
    basecalling.read_quality_boxplot()
    basecalling.channel_count_histogram()
    basecalling.read_number_run()

    if is_barcode:
        basecalling.barcode_percentage_pie_chart()
        basecalling.barcode_length_boxplot()
        basecalling.barcoded_phred_score_frequency()
    basecalling.read_length_histogram()

    channel_count = basecalling.channel
    total_number_reads_per_pore = pd.value_counts(channel_count)
    basecalling.plot_performance(total_number_reads_per_pore)
    basecalling.occupancy_pore()

    basecalling.phred_score_frequency()
    basecalling.scatterplot()

def statistics_log_file(fast5_data, basecalling, result_directory,is_barcode):
    '''
    Production of statistics file in the form of a tsv file
    :param fast5_data: tuple containing the informations extracted from a fast5 file
    :param basecalling: basecalling_stat_plotter instance
    :param result_directory: result directory
    :param is_barcode: boolean indicating if we use the barcodes or not
    '''
    if is_barcode:
        basecalling.statistics_dataframe()
    else:
        log_file1D.log_file1D(fast5_data, basecalling, result_directory)
        log_file1D.log_file_tsv(fast5_data, basecalling, result_directory)

def get_barcode(design_file_directory):
    '''
    Get the barcode from a file given in input
    :param design_file_directory: sample sheet directory
    :return: sorted list containing the barcode indicated in the sample sheet
    '''
    barcode_file = design_file_directory+"design.csv"
    set_doublon = set()

    with open(barcode_file) as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t')

        for row in spamreader:
            pattern = re.search(r'BC(\d{2})', row[0])

            if pattern:
                barcode = 'barcode{}'.format(pattern.group(1))
                set_doublon.add(barcode)
    return sorted(set_doublon)

def main():

    #Initialization of the differents directories used by the program
    run_name, config_file, is_barcode, file_list = get_args()
    dico_path = config_file_initialization(config_file, is_barcode, run_name, file_list)
    if not dico_path:
        sys.exit("Error, dico_path is empty")
    result_directory = dico_path['result_directory']
    fastq_directory  = dico_path['fastq_directory']
    fast5_directory  = dico_path['fast5_directory']
            
    if is_barcode:
        design_file_directory = dico_path['design_file_directory']
        barcode_selection = get_barcode(design_file_directory)
    else:
        barcode_selection = ''

    basecall_log = dico_path['basecall_log'] +run_name+'/sequencing_summary.txt'


    if os.path.isdir(result_directory):
        shutil.rmtree(result_directory, ignore_errors=True)
        os.makedirs(result_directory)
        
    else:
        os.makedirs(result_directory)

    #Determination of fast5 and fastq files extension
    dico_extension = extension(config_file, file_list, is_barcode)
    fast5_file_extension = dico_extension['fast5_file_extension']
    fastq_file_extension = dico_extension['fastq_file_extension']


    fast5_data = fast5_data_extractor.fast5_data_extractor(fast5_directory, result_directory, fast5_file_extension, run_name)
    basecalling = basecalling_stat_plotter1D.basecalling_stat_plotter1D(basecall_log, is_barcode,result_directory, fastq_directory, fast5_data, run_name, is_barcode, fastq_file_extension, barcode_selection)

    #Date and flowcell id extracted from a FAST5 file
    flowcell_id, *_ = fast5_data

    #Graph creation
    graph_creation(basecalling, is_barcode)

    #Generation of the report
    html_report.html_report(result_directory, basecalling.run_date(), flowcell_id, is_barcode, basecalling.sequence_length_template)

    #Creation of the statistics files
    statistics_log_file(fast5_data, basecalling, result_directory, is_barcode)

main()