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
import argparse
import os
from toulligqc import fastq_extractor
from toulligqc import fast5_extractor
from toulligqc import statistics_generator
from toulligqc import html_report
from toulligqc import version
from toulligqc import albacore_stats_extractor
from pathlib import Path
from toulligqc import toullig_conf


def parse_args(config_dictionary):
    '''
    Parsing the command line
    :return: config_dictionary containing the paths containing in the configuration file or specify by line arguments
    '''

    home = str(Path.home())
    parser = argparse.ArgumentParser()

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
    '''
    Check the configuration
    :param config_dictionary: configuration dictionary containing the file or directory paths
    '''
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

def get_barcode(samplesheet):
    '''
    Get the barcode from a file given in input
    :param samplesheet: sample sheet directory
    :return: sorted list containing the barcode indicated in the sample sheet
    '''
    barcode_file = samplesheet

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
    '''
    Main function creating graphs and statistics
    '''
    config_dictionary = toullig_conf.toullig_conf()
    parse_args(config_dictionary)
    check_conf(config_dictionary)


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

    html_report.html_report(config_dictionary, result_dict)

    statistics_generator.statistics_generator(config_dictionary, result_dict)


if __name__ == "__main__":
    main()