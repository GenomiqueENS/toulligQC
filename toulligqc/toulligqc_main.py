#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')
import shutil
import sys
import csv
import re
import configparser
import argparse
import glob
import pandas as pd
import os
from toulligqc import basecalling_stat_plotter1D
from toulligqc import fast5_data_extractor
from toulligqc import log_file1D
from toulligqc import html_report
from toulligqc import version


def get_args():
    '''
    Parsing the command line
    :return: different informations: run name,
    path towards files contained in the configuration file,
    boolean indicating if we use the barcode,
    path list used with the -f option
    '''
    conf_parser = argparse.ArgumentParser(
        # Turn off help, so we print all options in response to -h
        add_help=False
    )
    conf_parser.add_argument("-c", "--conf_file",
                             help="Specify config file", metavar="FILE")
    args, remaining_argv = conf_parser.parse_known_args()

    defaults = {
        "fast5_source": "no fast5_source",
        "albacore_summary.source": "no albacore summary source",
        "fastq_source": "no fastq source",
        "output": "no output",
        "sample_sheet_file": "no sample sheet source"
    }

    if args.conf_file:
        config = configparser.SafeConfigParser()
        config.read([args.conf_file])
        defaults = dict(config.items("config"))

    # Don't surpress add_help here so it will handle -h
    parser = argparse.ArgumentParser(
        # Inherit options from config_parser
        parents=[conf_parser],
        # print script description with -h/--help
        description=__doc__,
        # Don't mess with format of description
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # print(remaining_argv)
    parser.set_defaults(**defaults)
    parser.add_argument("-n", "--run_name", action='store', dest="run_name", help="Run name", required=True)

    parser.add_argument('-f', '--fast5-source', action='store', dest='fast5_source', help='Fast5 file source')
    parser.add_argument('-a', '--albacore-summary-source', action='store', dest='albacore_summary_source',
                        help='Albacore summary source')
    parser.add_argument('-q', '--fastq-source', action='store', dest='fastq_source', help='Fastq file source')
    parser.add_argument('-o', '--output', action='store', dest='output', help='Output directory')
    parser.add_argument('-s', '--sample-sheet-file', action='store', dest='sample_sheet_file',
                        help='Path to sample sheet file')
    parser.add_argument("-b", "--barcoding", action='store_true', dest='is_barcode', help="Barcode usage",
                        default=False)
    parser.add_argument('--version', action='version', version=version.__version__)

    argument_value = parser.parse_args(remaining_argv)
    print(argument_value)
    fast5_source = argument_value.fast5_source
    albacore_summary_source = argument_value.albacore_summary_source
    fastq_source = argument_value.fastq_source
    run_name = argument_value.run_name

    is_barcode = argument_value.is_barcode
    output_directory = argument_value.output
    sample_sheet_file = argument_value.sample_sheet_file

    if not fast5_source:
        print('The fast5 source argument is empty')
        sys.exit(0)

    elif not fastq_source:
        print('The fastq source arugment is empty')
        sys.exit(0)

    elif not albacore_summary_source:
        print('The albacore summary source argument is empty')
        sys.exit(0)

    elif is_barcode:
        if not sample_sheet_file:
            print('The sample sheet source argument is empty')
            sys.exit(0)

    elif not output_directory:
        print('The output directory argument is empty')
        sys.exit(0)

    else:
        pass

    return run_name, is_barcode, fast5_source, fastq_source, albacore_summary_source, sample_sheet_file, output_directory


def config_file_initialization(is_barcode, run_name, fast5_source='', fastq_source='', albacore_summary_source='',
                               sample_sheet_file='', output_directory=''):
    """
    Creation of a dictionary of the path file contained in the configuration file

    :param is_barcode: boolean
    :param run_name: Run name
    :param list_file: path list used with the -f option
    :return: Dictionary containing the file paths included into the configuration file
    """
    dico_path = {}
    dico_path['fast5_source'] = fast5_source
    dico_path['basecall_log_source'] = albacore_summary_source
    dico_path['fastq_source'] = fastq_source
    dico_path['result_directory'] = output_directory

    if is_barcode:
        dico_path['design_file'] = sample_sheet_file

    if dico_path['result_directory'] == '':
        dico_path['result_directory'] = os.getcwd()

    for key, value in dico_path.items():
        if value.endswith('/'):
            continue
        elif os.path.isfile(value):
            continue
        else:
            dico_path[key] = value + '/'

    if not os.path.isdir(dico_path['result_directory']):
        os.makedirs(dico_path['result_directory'])

    if os.path.isdir(dico_path['result_directory'] + run_name):
        shutil.rmtree(dico_path['result_directory'] + run_name, ignore_errors=True)
        os.makedirs(dico_path['result_directory'] + run_name)

    else:
        os.makedirs(dico_path['result_directory'] + run_name)

    dico_path['result_directory'] = dico_path['result_directory'] + run_name + '/'

    return dico_path


def extension(run_name, is_barcode, fast5_source='', fastq_source='', sample_sheet_file='',
              albacore_summary_source='', output_directory=''):
    '''
    Creation of a dictionary containing the extension used for the fast5 and fastq files
    :param list_file: path list used with the -f option
    :return: Extension dictionary used for the fast5 and fastq files
    '''

    dico_extension = {}

    if os.path.isdir(fast5_source):
        dico_extension['fast5_file_extension'] = 'fast5_directory'

    elif fast5_source.endswith('.tar.gz'):
        dico_extension['fast5_file_extension'] = 'tar.gz'

    elif fast5_source.endswith('.fast5'):
        dico_extension['fast5_file_extension'] = 'fast5'

    elif fast5_source.endswith('.tar.bz2'):
        dico_extension['fast5_file_extension'] = 'tar.bz2'

    else:
        print('The fast5 extension is not supported (fast5, tar.bz2 or tar.gz format)')
        sys.exit(0)
    print(fastq_source)
    if os.path.isdir(fastq_source):
        fastq_directory = fastq_source + run_name + '/'

        if glob.glob(fastq_directory + '/*.fastq'):
            dico_extension['fastq_file_extension'] = 'fastq'

        elif glob.glob(fastq_directory + '/*.fq'):
            dico_extension['fastq_file_extension'] = 'fq'

        elif glob.glob(fastq_directory + '/*.gz'):
            dico_extension['fastq_file_extension'] = 'gz'

        elif glob.glob(fastq_directory + '/*.bz2'):
            dico_extension['fastq_file_extension'] = 'bz2'

        else:
            print('The fastq source extension is not supported (fast5, tar.bz2 or tar.gz format)')
            sys.exit(0)

    elif fastq_source.endswith('.fastq'):
        dico_extension['fastq_file_extension'] = 'fastq'

    elif fastq_source.endswith('.fq'):
        dico_extension['fastq_file_extension'] = 'fq'

    elif fastq_source.endswith('.bz2') or fastq_source.endswith('.gz') or fastq_source.endswith('.zip'):
        pattern = '\.(gz|bz2|zip)$'
        if re.search(pattern, fastq_source):
            match = re.search(pattern, fastq_source)
            dico_extension['fastq_file_extension'] = match.groups()[0]

    else:
        print('The fastq source extension is not supported (fast5, bz2 or gz format)')
        sys.exit(0)

    fast5_file_extension = dico_extension['fast5_file_extension']
    print(dico_extension['fast5_file_extension'])
    #if  fast5_file_extension != 'fast5' or fast5_file_extension != 'tar.bz2' or fast5_file_extension != 'tar.gz' or fast5_file_extension != 'fast5_directory':
     #   print('The fast5 source extension is not supported (fast5, tar.bz2 or tar.gz format) or delete the . in front of the extension name(bz2 and not .bz2)')
      #  sys.exit(0)

    #if dico_extension['fastq_file_extension'] != 'fastq' or dico_extension['fastq_file_extension'] != 'bz2' or \
     #               dico_extension['fastq_file_extension'] != 'gz':
     #   print('The fastq source extension is not supported (fast5, bz2 or gz format) or delete the . in \
      #        front of the extension name(bz2 and not .bz2)')
       # sys.exit(0)

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


def statistics_log_file(fast5_data, basecalling, result_directory, is_barcode):
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
    # Initialization of the differents directories used by the program
    run_name, is_barcode, fast5_source, fastq_source, albacore_summary_source, sample_sheet_file, output_directory = get_args()
    dico_path = config_file_initialization(is_barcode, run_name, fast5_source, fastq_source, albacore_summary_source,
                                           sample_sheet_file, output_directory)
    if not dico_path:
        sys.exit("Error, dico_path is empty")

    result_directory = dico_path['result_directory']
    fastq_directory = dico_path['fastq_source']
    fast5_directory = dico_path['fast5_source']

    if is_barcode:
        design_file = dico_path['design_file']
        barcode_selection = get_barcode(design_file)

        if barcode_selection == '':
            print("Sample sheet is empty")
            sys.exit(0)
    else:
        barcode_selection = ''

    if os.path.isdir(dico_path['basecall_log_source']):
        albacore_summary_source = dico_path['basecall_log_source'] + run_name + '/sequencing_summary.txt'


    # Determination of fast5 and fastq files extension
    dico_extension = extension(run_name, is_barcode, fast5_directory, fastq_directory, sample_sheet_file,
                               albacore_summary_source, output_directory)
    print(dico_extension)
    fast5_file_extension = dico_extension['fast5_file_extension']
    fastq_file_extension = dico_extension['fastq_file_extension']

    fast5_data = fast5_data_extractor.fast5_data_extractor(fast5_directory, result_directory, fast5_file_extension,
                                                           run_name)
    basecalling = basecalling_stat_plotter1D.basecalling_stat_plotter1D(albacore_summary_source, is_barcode,
                                                                        result_directory, fastq_directory, fast5_data,
                                                                        run_name, is_barcode, fastq_file_extension,
                                                                        barcode_selection)

    # Date and flowcell id extracted from a FAST5 file
    flowcell_id, *_ = fast5_data

    # Graph creation
    graph_creation(basecalling, is_barcode)

    # Generation of the report
    html_report.html_report(result_directory, basecalling.run_date(), flowcell_id, is_barcode,
                            basecalling.sequence_length_template)

    # Creation of the statistics files
    statistics_log_file(fast5_data, basecalling, result_directory, is_barcode)


main()
