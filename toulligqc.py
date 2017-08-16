#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
import glob

__all__ = ['main']

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
    parser.add_argument("-b","--barcoding", action='store_true',dest='is_barcode',help="Barcode usage",default=False)
    parser.add_argument("-c", "--config_file", action='store', dest='config_file', help="Path to the configuration file", default=home_path+'/.toulligqc/config.txt')
    parser.add_argument('-f','--fast5-source', action='store', dest='fast5_source', help='Fast5 file source',default='')
    parser.add_argument('-a','--albacore-summary-source', action='store', dest='albacore_summary_source', help='Albacore summary source', default='')
    parser.add_argument('-q','--fastq-source', action='store', dest='fastq_source', help='fastq file source', default='')
    parser.add_argument('-o', '--output', action='store', dest='output_directory', help='output directory', default='')
    parser.add_argument('-s', '--sample-sheet-source', action='store', dest='sample_sheet_source', help='Sample sheet source')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    argument_value=parser.parse_args()
    fast5_source = argument_value.fast5_source
    albacore_summary_source = argument_value.albacore_summary_source
    fastq_source = argument_value.fastq_source
    run_name = argument_value.run_name

    if fast5_source and fastq_source:
        config_file_path = ''
    else:
        config_file_path = argument_value.config_file
    is_barcode = argument_value.is_barcode
    output_directory = argument_value.output_directory
    sample_sheet_source = argument_value.sample_sheet_source
    return run_name, config_file_path, is_barcode, fast5_source, fastq_source, albacore_summary_source, sample_sheet_source, output_directory


def config_file_initialization(is_barcode, run_name, config_file = '', fast5_source  = '', fastq_source = '', albacore_summary_source = '', sample_sheet_source = '', output_directory = '' ):
    """
    Creation of a dictionary of the path file contained in the configuration file

    :param config_file: configuration file path
    :param is_barcode: boolean
    :param run_name: Run name
    :param list_file: path list used with the -f option
    :return: Dictionary containing the file paths included into the configuration file
    """
    dico_path = {}

    if fast5_source and fastq_source and albacore_summary_source:
        dico_path = {}
        dico_path['fast5_source'] = fast5_source
        dico_path['basecall_log_source'] = albacore_summary_source
        dico_path['fastq_source'] = fastq_source
        dico_path['result_directory'] = output_directory

        if is_barcode:
            dico_path['design_file'] = sample_sheet_source


    elif config_file:
        configFilePath = config_file
        config = configparser.ConfigParser()
        config.read(configFilePath)
        dico_path['result_directory'] = config.get('config', 'result.directory')
        dico_path['basecall_log_source'] = config.get('config', 'albacore.summary.directory')
        dico_path['fastq_source'] = config.get('config', 'fastq.directory')
        dico_path['fast5_source'] = config.get('config', 'fast5.directory')

        if is_barcode:
            dico_path['design_file'] = config.get('config', 'design.file')

    else:
        print('Error, not a config file')
        return 0

    if dico_path['result_directory'] == '':
        dico_path['result_directory'] = os.getcwd()

    for key, value in dico_path.items():
        if value.endswith('/'):
            continue
        elif key == 'design_file':
            continue
        elif os.path.isfile(value):
            continue
        else:
            dico_path[key] = value + '/'

    dico_path['result_directory'] = dico_path['result_directory'] + run_name + '/'

    return dico_path


def extension(is_barcode, config_file = '', fast5_source = '', fastq_source = '', sample_sheet_source = '', albacore_summary_source = '', output_directory = ''):
    '''
    Creation of a dictionary containing the extension used for the fast5 and fastq files
    :param config_file: configuration file
    :param list_file: path list used with the -f option
    :return: Extension dictionary used for the fast5 and fastq files
    '''

    dico_extension = {}
    if fast5_source and fastq_source:
        if os.path.isdir(fast5_source):
            dico_extension['fast5_file_extension'] = 'fast5_directory'

        elif fast5_source.endswith('.fast5'):
            dico_extension['fast5_file_extension'] = 'fast5'

        elif fast5_source.endswith('.tar.bz2'):
            dico_extension['fast5_file_extension'] = 'tar.bz2'

        elif fast5_source.endswith('.tar.gz'):
            dico_extension['fast5_file_extension'] = 'tar.gz'
#Si c'est un répertoire que fait-t-on ?
            #A voir
        if os.path.isdir(fastq_source):
            fastq_directory = fastq_source
            if glob.glob(fastq_directory+'*.fastq'):
                dico_extension['fastq_file_extension'] = 'fastq'

            elif glob.glob(fastq_directory+'*.fq'):
                dico_extension['fastq_file_extension'] = 'fq'

            elif glob.glob(fastq_directory+'*.gz'):
                dico_extension['fastq_file_extension'] = 'gz'

            elif glob.glob(fastq_directory+'*.bz2'):
                dico_extension['fastq_file_extension'] = 'bz2'

        if fastq_source.endswith('.fastq'):
            dico_extension['fastq_file_extension'] = 'fastq'

        elif fastq_source.endswith('.fq'):
            dico_extension['fastq_file_extension'] = 'fq'

        elif fastq_source.endswith('.bz2') or fastq_source.endswith('.gz') or fastq_source.endswith('.zip'):
            pattern = '\.(gz|bz2|zip)$'
            if re.search(pattern, fastq_source):
                match = re.search(pattern, fastq_source)
                dico_extension['fastq_file_extension'] = match.groups()[0]
            
    elif config_file:
        config = configparser.ConfigParser()
        config.read(config_file)
        dico_extension['fast5_file_extension'] = config.get('extension', 'fast5.file.extension')
        dico_extension['fastq_file_extension'] = config.get('extension', 'fastq.file.extension')

    else:
        print('Pas de fichier source utilisé')
        return 0

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

def get_barcode(design_file):
    '''
    Get the barcode from a file given in input
    :param design_file: sample sheet directory
    :return: sorted list containing the barcode indicated in the sample sheet
    '''
    barcode_file = design_file
    set_doublon = set()

    with open(barcode_file) as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t')

        for row in spamreader:
            if row[0].startswith('#'):
                continue
            else:
                pattern = re.search(r'BC(\d{2})', row[0])

            if pattern:
                barcode = 'barcode{}'.format(pattern.group(1))
                set_doublon.add(barcode)
    return sorted(set_doublon)

def main():

    #Initialization of the differents directories used by the program
    run_name, config_file, is_barcode, fast5_source, fastq_source, albacore_summary_source,sample_sheet_source, output_directory  = get_args()
    dico_path = config_file_initialization(is_barcode, run_name, config_file, fast5_source, fastq_source, albacore_summary_source, sample_sheet_source, output_directory)
    if not dico_path:
        sys.exit("Error, dico_path is empty")


    result_directory = dico_path['result_directory']
    fastq_directory  = dico_path['fastq_source']
    fast5_directory  = dico_path['fast5_source']
    if is_barcode:
        design_file = dico_path['design_file']
        barcode_selection = get_barcode(design_file)
    else:
        barcode_selection = ''

    if os.path.isdir(dico_path['basecall_log_source']):
        albacore_summary_source  = dico_path['basecall_log_source'] +run_name+'/sequencing_summary.txt'


    if os.path.isdir(result_directory):
        shutil.rmtree(result_directory, ignore_errors=True)
        os.makedirs(result_directory)
        
    else:
        os.makedirs(result_directory)

    #Determination of fast5 and fastq files extension
    dico_extension = extension(is_barcode, config_file, fast5_source, fastq_source, sample_sheet_source, albacore_summary_source, output_directory)
    fast5_file_extension = dico_extension['fast5_file_extension']
    fastq_file_extension = dico_extension['fastq_file_extension']


    fast5_data = fast5_data_extractor.fast5_data_extractor(fast5_directory, result_directory, fast5_file_extension, run_name, config_file)
    basecalling = basecalling_stat_plotter1D.basecalling_stat_plotter1D(albacore_summary_source, is_barcode,result_directory, fastq_directory, fast5_data, run_name, is_barcode, fastq_file_extension, config_file, barcode_selection)

    #Date and flowcell id extracted from a FAST5 file
    flowcell_id, *_ = fast5_data

    #Graph creation
    graph_creation(basecalling, is_barcode)

    #Generation of the report
    html_report.html_report(result_directory, basecalling.run_date(), flowcell_id, is_barcode, basecalling.sequence_length_template)

    #Creation of the statistics files
    statistics_log_file(fast5_data, basecalling, result_directory, is_barcode)

main()
