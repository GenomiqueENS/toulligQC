#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
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
#
# For more information on the ToulligQC project and its aims,
# visit the home page at:
#
#      https://github.com/GenomicParisCentre/toulligQC
#
# First author: Lionel Ferrato-Berberian
# Maintainer: Bérengère Laffay
# Since version 0.1
#
# Toulligqc.py: Main's constitution
# 1. It scans the user command line and creates the output directory
# 2. Depending on the options selected, it creates a list of the necessary extractors
# 3. For each extractor, it successively call the _check_conf, init, extract, graph_generation and clean methods
# to fill in the result_dict dictionary
# 4. In the case of barcoded sequencing, it searches all barcodes from the command line argument --barcodes
# 5. It uses all the information collected to generate a qc in the form of a htl-report and a report.data file

import matplotlib

matplotlib.use('Agg')
import shutil
import sys
import re
import argparse
import os
import time
import datetime

import warnings
from toulligqc import toulligqc_info_extractor
from toulligqc import report_data_file_generator
from toulligqc import html_report_generator
from toulligqc import version
from toulligqc import configuration
from toulligqc import fast5_extractor
from toulligqc import sequencing_summary_extractor
from toulligqc import sequencing_summary_onedsquare_extractor
from toulligqc import sequencing_telemetry_extractor
from toulligqc import common


def _parse_args(config_dictionary):
    """
    Parsing the command line
    :return: config_dictionary containing the paths specified by line arguments
    """

    parser = argparse.ArgumentParser(prog="ToulligQC V{0}".format(version.__version__), add_help=False)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    # Add all required arguments
    required.add_argument('-a', '--sequencing-summary-source', action='append', dest='sequencing_summary_source',
                          help='Basecaller sequencing summary source', metavar='SEQUENCING_SUMMARY_SOURCE',
                          required=True)
    required.add_argument('-t', '--telemetry-source', action='store', dest='telemetry_source',
                          help='Basecaller telemetry file source', default=False)

    required.add_argument('-f', '--fast5-source', action='store', dest='fast5_source',
                          help='Fast5 file source (necessary if no telemetry file)')

    # Add all optional arguments
    optional.add_argument("-n", "--report-name", action='store', dest="report_name", help="Report name", type=str)
    optional.add_argument('--output-directory', action='store', dest='output', help='Output directory')
    optional.add_argument('-o', '--html-report-path', action='store', dest='html_report_path',
                          help='Output HTML report')
    optional.add_argument('--data-report-path', action='store', dest='data_report_path', help='Output data report')
    optional.add_argument('--images-directory', action='store', dest='images_directory', help='Images directory')
    optional.add_argument('-d', '--sequencing-summary-1dsqr-source', action='append',
                          dest='sequencing_summary_1dsqr_source',
                          help='Basecaller 1dsq summary source')
    optional.add_argument("-b", "--barcoding", action='store_true', dest='is_barcode', help="Option for barcode usage",
                          default=False)
    optional.add_argument('-l', '--barcodes', action='store', default='', dest='barcodes',
                          help='Coma separated barcode list')
    optional.add_argument("--quiet", action='store_true', dest='is_quiet', help="Quiet mode",
                          default=False)
    optional.add_argument("--report-only", action='store_true', dest='report_only',
                          help=argparse.SUPPRESS,
                          default=False)
    optional.add_argument("--force", action='store_true', dest='force', help="Force overwriting of existing files",
                          default=False)
    optional.add_argument("--debug", action='store_true', dest='debug', help=argparse.SUPPRESS,
                          default=False)
    optional.add_argument("-h", "--help", action="help", help="Show this help message and exit")
    optional.add_argument('--version', action='version', version=version.__version__)

    # Parsing lone arguments and assign each argument value to a variable
    args = parser.parse_args()
    report_name = args.report_name
    is_barcode = args.is_barcode
    barcodes = args.barcodes

    # If a barcode list is provided, automatically add --barcoding argument
    if len(barcodes) > 0:
        is_barcode = True

    # If no report_name specified, create default one : ToulligQC-report-YYYYMMDD_HHMMSS
    if not report_name:
        timestamp = datetime.datetime.now()
        config_dictionary['report_name'] = "Toulligqc-report-" + str((timestamp.strftime("%Y-%m-%d-%H%M%S")))
    else:
        config_dictionary['report_name'] = report_name

    # Rewrite the configuration file value if argument option is present
    args_dict = {
        ('fast5_source', args.fast5_source),
        ('sequencing_summary_source', _join_parameter_arguments(args.sequencing_summary_source)),
        ('sequencing_summary_1dsqr_source', _join_parameter_arguments(args.sequencing_summary_1dsqr_source)),
        ('sequencing_telemetry_source', args.telemetry_source),
        ('result_directory', args.output),
        ('html_report_path', args.html_report_path),
        ('data_report_path', args.data_report_path),
        ('images_directory', args.images_directory),
        ('barcoding', is_barcode),
        ('barcodes', barcodes),
        ('quiet', args.is_quiet),
        ('report_only', args.report_only),
        ('force', args.force),
        ('debug', args.debug)
    }

    # Put arguments values in configuration object
    for key, value in args_dict:
        if value:
            config_dictionary[key] = value

    # Directory paths must ends with '/'
    for key, value in config_dictionary.items():
        if type(value) == str and (key.endswith('_source') or key.endswith('_directory')) and os.path.isdir(value) and (
                not value.endswith('/')):
            config_dictionary[key] = value + '/'

    # Convert all configuration values in strings
    for key, value in config_dictionary.items():
        config_dictionary[key] = str(value)

    return config_dictionary


def _check_conf(config_dictionary):
    """
    Check the configuration
    :param config_dictionary: configuration dictionary containing the file or directory paths
    """

    force = True if config_dictionary.get('force', 'False').lower() == 'true' else False

    if ('sequencing_summary_source' not in config_dictionary or not config_dictionary['sequencing_summary_source']) and \
            ('sequencing_telemetry_source' not in config_dictionary or not config_dictionary[
                'sequencing_telemetry_source']):
        argparse.ArgumentParser.print_help

    if 'sequencing_summary_source' not in config_dictionary or not config_dictionary['sequencing_summary_source']:
        sys.exit('ERROR: The sequencing summary file argument is empty')

    if 'html_report_path' not in config_dictionary or not config_dictionary['html_report_path']:

        # If no --output argument provided, create output folder in current directory
        if 'result_directory' not in config_dictionary or not config_dictionary['result_directory']:
            current_directory = os.getcwd()
            config_dictionary['result_directory'] = current_directory + '/'

        # Create the root output directory if not exists
        if not os.path.isdir(config_dictionary['result_directory']):
            os.makedirs(config_dictionary['result_directory'])

        # Define the output directory
        config_dictionary['result_directory'] = \
            config_dictionary['result_directory'] + config_dictionary['report_name'] + '/'

        _check_if_dir_exists(config_dictionary['result_directory'], force)

        # Define the output paths
        config_dictionary['images_directory'] = config_dictionary['result_directory'] + 'images/'
        config_dictionary['html_report_path'] = config_dictionary['result_directory'] + 'report.html'
        config_dictionary['data_report_path'] = config_dictionary['result_directory'] + 'report.data'
        del config_dictionary['result_directory']

    if 'images_directory' not in config_dictionary:
        config_dictionary['images_directory'] = None

    if 'data_report_path' not in config_dictionary:
        config_dictionary['data_report_path'] = None

    _check_if_dir_exists(config_dictionary['images_directory'], force)
    _check_if_file_exists(config_dictionary['html_report_path'], force)
    _check_if_file_exists(config_dictionary['data_report_path'], force)

    print(config_dictionary['html_report_path'])


def _check_if_dir_exists(dir, force):
    if dir is None:
        return

    if os.path.isdir(dir):
        if not force:
            sys.exit("Error directory already exists: " + dir)
        else:
            shutil.rmtree(dir, ignore_errors=True)
    os.makedirs(dir)


def _check_if_file_exists(path, force):
    if path is None:
        return

    if os.path.isfile(path):

        if not force:
            sys.exit("Error file already exists: " + path)
        else:
            os.remove(path)


def _welcome(config_dictionary):
    """
    Print welcome message
    """
    _show(config_dictionary, "ToulligQC version " + config_dictionary['app.version'])


def _show(config_dictionary, msg):
    """
    Print a message on the screen
    :param config_dictionary: configuration dictionary
    :param msg: message to print
    """
    if 'quiet' not in config_dictionary or config_dictionary['quiet'].lower() != 'true':
        print(msg)


def _join_parameter_arguments(arg):
    """
    Join parameter arguments
    :param arg: argument to join
    :return: a string with arguments separated by tab character or None if the input parameter is None
    """

    if (arg is None):
        return None
    return '\t'.join(arg)


def _create_extractor_list(config_dictionary):
    result = []

    if 'sequencing_telemetry_source' in config_dictionary and \
            config_dictionary['sequencing_telemetry_source']:
        result.append(sequencing_telemetry_extractor.SequencingTelemetryExtractor(config_dictionary))

    if 'fast5_source' in config_dictionary and config_dictionary['fast5_source']:
        result.append(fast5_extractor.Fast5Extractor(config_dictionary))

    if 'sequencing_summary_1dsqr_source' in config_dictionary and \
            config_dictionary['sequencing_summary_1dsqr_source']:
        result.append(sequencing_summary_onedsquare_extractor.
                      OneDSquareSequencingSummaryExtractor(config_dictionary))
    else:
        result.append(sequencing_summary_extractor.SequencingSummaryExtractor(config_dictionary))

    result.insert(0, toulligqc_info_extractor.ToulligqcInfoExtractor(config_dictionary, result))

    return result


def main():
    """
    Main function creating graphs and statistics
    """
    config_dictionary = configuration.ToulligqcConf()
    _parse_args(config_dictionary)
    _check_conf(config_dictionary)

    warnings.simplefilter('ignore')

    if not config_dictionary:
        sys.exit("ERROR: dico_path is empty")

    # Get barcode selection
    if config_dictionary['barcoding'].lower() == 'true':
        config_dictionary['barcode_selection'] = []

        if 'barcodes' in config_dictionary:
            barcode_set = set()
            for b in config_dictionary['barcodes'].strip().split(','):
                pattern = re.search(r'(BC|RB|NB|BARCODE)(\d{2})', b.strip().upper())
                if pattern:
                    barcode = 'barcode{}'.format(pattern.group(2))
                    barcode_set.add(barcode)
            barcode_selection = sorted(barcode_set)

            if len(barcode_selection) == 0:
                sys.exit("ERROR: No known barcode found in provided list of barcodes")
            config_dictionary['barcode_selection'] = barcode_selection
    else:
        config_dictionary['barcode_selection'] = ''

    # Print welcome message
    _welcome(config_dictionary)

    # Configuration checking and initialisation of the extractors
    _show(config_dictionary, "* Initialize extractors")

    # Create the list of extractors to execute
    extractors_list = _create_extractor_list(config_dictionary)

    # Check extractor configuration
    for extractor in extractors_list:
        (check_result, error_message) = extractor.check_conf()
        if not check_result:
            sys.exit("ERROR: Error while checking " + extractor.get_name() + " configuration: " + error_message)

    result_dict = {}
    graphs = []
    qc_start = time.time()

    # Information extraction about statistics and generation of the graphs
    for extractor in extractors_list:
        _show(config_dictionary, "* Start {0} extractor".format(extractor.get_name()))
        extractor_start = time.time()

        # Execute extractor
        extractor.init()
        extractor.extract(result_dict)
        graphs.extend(extractor.graph_generation(result_dict))
        extractor.clean(result_dict)

        extractor_end = time.time()
        extract_time = extractor_end - extractor_start
        result_dict['{}.duration'.format(extractor.get_report_data_file_id())] = round(extract_time, 2)

        _show(config_dictionary, "* End of {0} extractor (done in {1})".format(extractor.get_name(),
                                                                               common.format_duration(extract_time)))

    # HTML report and report.data file generation
    _show(config_dictionary, "* Write HTML report")
    html_report_generator.html_report(config_dictionary, result_dict, graphs)

    qc_end = time.time()
    result_dict['toulligqc.info.execution.duration'] = round((qc_end - qc_start), 2)

    if config_dictionary['report_only'].lower() != 'true':
        _show(config_dictionary, "* Write statistics files")
        report_data_file_generator.statistics_generator(config_dictionary, result_dict)
    _show(config_dictionary, "* End of the QC extractor (done in {})".format(common.format_duration(qc_end - qc_start)))


if __name__ == "__main__":
    main()
