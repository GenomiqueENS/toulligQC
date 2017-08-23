#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#                  ToulligQC development code
#
# This code may be freely distributed and modified under the
# terms of the GNU General Public License version 3 or later
# and CeCILL. This should be distributed with the code. If you
# do not have a copy, see:
#
#      http://www.gnu.org/licenses/gpl-3.0-standalone.html
#      http://www.cecill.info/licences/Licence_CeCILL_V2-en.html
#
# Copyright for this code is held jointly by the Genomic platform
# of the Institut de Biologie de l'École Normale Supérieure and
# the individual authors.

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
import basecalling_stat_plotter1D
import fastq_extractor
import fast5_extractor
import log_file1D
import html_report
import version
import albacore_stats_extractor
from pathlib import Path
import config


def parse_args(config_dictionary):
    '''
    Parsing the command line
    :return: different informations: run name,
    path towards files contained in the configuration file,
    boolean indicating if we use the barcode,
    path list used with the -f option
    '''

    home = str(Path.home())
    parser = argparse.ArgumentParser()
    # print(remaining_argv)
    #parser.set_defaults(**defaults)
    parser.add_argument("-c", "--conf-file",
                             help="Specify config file", metavar="FILE")
    parser.add_argument("-n", "--run-name", action='store', dest="run_name", help="Run name", required=True)

    parser.add_argument('-f', '--fast5-source', action='store', dest='fast5_source', help='Fast5 file source')
    parser.add_argument('-a', '--albacore-summary-source', action='store', dest='albacore_summary_source',
                        help='Albacore summary source')
    parser.add_argument('-q', '--fastq-source', action='store', dest='fastq_source', help='Fastq file source')
    parser.add_argument('-o', '--output', action='store', dest='output', help='Output directory')
    parser.add_argument('-s', '--samplesheet-file', action='store', dest='sample_sheet_file',
                        help='Path to sample sheet file')
    parser.add_argument("-b", "--barcoding", action='store_true', dest='is_barcode', help="Barcode usage",
                        default=False)
    parser.add_argument('--version', action='version', version=version.__version__)

    argument_value = parser.parse_args()
    conf_file = argument_value.conf_file
    fast5_source = argument_value.fast5_source
    albacore_summary_source = argument_value.albacore_summary_source
    fastq_source = argument_value.fastq_source
    run_name = argument_value.run_name
    is_barcode = argument_value.is_barcode
    result_directory = argument_value.output
    sample_sheet_file = argument_value.sample_sheet_file

    config_dictionary['run_name'] = run_name
    # Load confiuration file if possible
    if argument_value.conf_file:
        config_dictionary.load(conf_file)
    elif os.path.isfile(home + '/.toulligqc/config.txt'):
        config_dictionary.load(home + '/.toulligqc/config.txt')
    if fast5_source:
        config_dictionary['fast5_source'] = fast5_source

    if albacore_summary_source:
        config_dictionary['albacore_summary_source'] = albacore_summary_source


    if fastq_source:
        config_dictionary['fastq_source'] = fastq_source
    else:
        config_dictionary['fastq_source'] = config_dictionary['fastq_source']+'/'+run_name
    if result_directory:
        config_dictionary['result_directory'] = result_directory

    if sample_sheet_file:
        config_dictionary['sample_sheet_file'] = sample_sheet_file

    if is_barcode:
        config_dictionary['barcoding'] = True
    for key, value in config_dictionary.items():

        if type(value)==bool:
            continue
        elif value.endswith('/'):
            continue
        elif os.path.isfile(value):
            continue
        elif value == '':
            continue
        elif os.path.isdir(value):
            config_dictionary[key] = value + '/'
        else:
            continue

    return config_dictionary

def check_conf(config_dictionary):
    if not config_dictionary['fast5_source']:
        print('The fast5 source argument is empty')
        sys.exit(0)

    elif not config_dictionary['fastq_source']:
        print('The fastq source arugment is empty')
        sys.exit(0)

    elif not config_dictionary['albacore_summary_source']:
        print('The albacore summary source argument is empty')
        sys.exit(0)

    elif config_dictionary['barcoding']:
        if not config_dictionary['sample_sheet_file']:
            print('The sample sheet source argument is empty')
            sys.exit(0)

    elif not config_dictionary['result_directory']:
        print('The output directory argument is empty')
        sys.exit(0)

    else:
        pass
    if not os.path.isdir(config_dictionary['result_directory']):
        os.makedirs(config_dictionary['result_directory'])

    if os.path.isdir(config_dictionary['result_directory'] + config_dictionary['run_name']):
        shutil.rmtree(config_dictionary['result_directory'] + config_dictionary['run_name'], ignore_errors=True)
        os.makedirs(config_dictionary['result_directory'] + config_dictionary['run_name'])

    else:
        os.makedirs(config_dictionary['result_directory'] + config_dictionary['run_name'])

    config_dictionary['result_directory'] = config_dictionary['result_directory'] + config_dictionary['run_name'] + '/'


def statistics_log_file(config_dictionary,result_dict):
    '''
    Production of statistics file in the form of a tsv file
    :param fast5_data: tuple containing the informations extracted from a fast5 file
    :param basecalling: basecalling_stat_plotter instance
    :param result_directory: result directory
    :param is_barcode: boolean indicating if we use the barcodes or not
    '''

    log_file1D.log_file1D(config_dictionary,result_dict)
    #log_file1D.log_file_tsv(fast5_data, basecalling, result_directory)


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

def create_output_directories(config_dictionary):
    return


def main():
    # Initialization of the differents directories used by the program
    config_dictionary = config.toullig_conf()
    parse_args(config_dictionary)
    check_conf(config_dictionary)

   # dico_path = config_file_initialization(is_barcode, run_name, fast5_source, fastq_source, albacore_summary_source,
    #                                       sample_sheet_file, output_directory)
    if not config_dictionary:
        sys.exit("Error, dico_path is empty")

    if config_dictionary['barcoding']:
        sample_sheet_file = config_dictionary['sample_sheet_file']
        barcode_selection = get_barcode(sample_sheet_file)
        config_dictionary['barcode_selection'] = barcode_selection
        if barcode_selection == '':
            print("Sample sheet is empty")
            sys.exit(0)
    else:
        config_dictionary['barcode_selection'] = ''

    if os.path.isdir(config_dictionary['albacore_summary_source']):
        config_dictionary['albacore_summary_source'] = config_dictionary['albacore_summary_source'] + config_dictionary['run_name'] + '/sequencing_summary.txt'

    #Create extractors objects
    extractors = (fast5_extractor.fast5_extractor(config_dictionary), fastq_extractor.fastq_extractor(config_dictionary), albacore_stats_extractor.albacore_stats_extractor(config_dictionary))
    for extractor in extractors:
        extractor.check_conf()
        extractor.init()

    result_dict = {}
    graphs = []

    for extractor in extractors:
        print(type(extractor))
        extractor.extract(result_dict)
        graphs.extend(extractor.graph_generation())
        extractor.clean()

    # Generation of the report
    html_report.html_report(config_dictionary, result_dict, graphs)

    # Creation of the statistics files
    statistics_log_file(config_dictionary, result_dict)


if __name__ == "__main__":
    main()