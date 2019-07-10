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
#
# For more information on the ToulligQC project and its aims,
# visit the home page at:
#
#      https://github.com/GenomicParisCentre/toulligQC
#
# First author: Lionel Ferrato-Berberian
# Maintainer: Bérengère Laffay
# Since version 0.1

# Extraction of statistics about the FASTQ files

import os
import glob
from collections import Counter
import pandas as pd
import bz2
import io
import gzip
import sys
import re


class FastqExtractor:
    """
    Uses fastq files to infer the nucleotide rate per read.
    Also gives the nucleotide rate per read per sample in case of barcoded samples
    param Fastq_source: Path to the Fastq files
    return: result_dict dictionary filled
    """
    def __init__(self, config_dictionary):
        self.config_dictionary = config_dictionary
        self.global_dico = {}
        self.global_length_array = []
        self.report_name = config_dictionary['report_name']
        self.is_barcode = config_dictionary['barcoding']
        self.get_report_data_file_id()

        if self.is_barcode == 'True':
            self.is_barcode = True
        elif self.is_barcode == 'False':
            self.is_barcode = False

        self.fastq_file = ''
        self.result_directory = config_dictionary['result_directory']
        self.barcode_selection = config_dictionary['barcode_selection']
        self.fastq_source = config_dictionary['fastq_source']
        self.selection_global = []

        self.fastq_file_extension = ''

    def check_conf(self):
        """Configuration checking"""

        if os.path.isdir(self.fastq_source):
            if glob.glob(self.fastq_source + '/*.fastq') \
                    or glob.glob(self.fastq_source + self.report_name + '/*.fastq'):
                self.fastq_file_extension = 'fastq'

            elif glob.glob(self.fastq_source + '/*.fq') or glob.glob(self.fastq_source + self.report_name + '/*.fq'):
                self.fastq_file_extension = 'fq'

            elif glob.glob(self.fastq_source + '/*.gz') or glob.glob(self.fastq_source + self.report_name + '/*.gz'):
                self.fastq_file_extension = 'gz'

            elif glob.glob(self.fastq_source + '/*.bz2') or glob.glob(self.fastq_source + self.report_name + '/*.bz2'):
                self.fastq_file_extension = 'bz2'

            elif glob.glob(self.fastq_source + self.report_name + '/*.fastq'):
                self.fastq_source = self.fastq_source + self.report_name + '/'
                self.fastq_file_extension = 'fastq'

            elif glob.glob(self.fastq_source + self.report_name + '/*.fq'):
                self.fastq_source = self.fastq_source + self.report_name + '/'
                self.fastq_file_extension = 'fq'

            elif glob.glob(self.fastq_source + self.report_name + '/*.gz'):
                self.fastq_source = self.fastq_source + self.report_name + '/'
                self.fastq_file_extension = 'gz'

            elif glob.glob(self.fastq_source + self.report_name + '/*.bz2'):
                self.fastq_source = self.fastq_source + self.report_name + '/'
                self.fastq_file_extension = 'bz2'

            else:
                return False, 'The fastq source extension is not supported (fast5, tar.bz2 or tar.gz format)'

        elif self.fastq_source.endswith('.fastq'):
            self.fastq_file_extension = 'fastq'

        elif self.fastq_source.endswith('.fq'):
            self.fastq_file_extension = 'fq'

        elif self.fastq_source.endswith('.bz2') or self.fastq_source.endswith('.gz') or self.fastq_source.endswith(
                '.zip'):
            pattern = '\.(gz|bz2|zip)$'
            if re.search(pattern, self.fastq_source):
                match = re.search(pattern, self.fastq_source)
                self.fastq_file_extension = match.groups()[0]
        else:
            return False, 'The fastq source extension is not supported (fast5, bz2 or gz format)'

        return True, ""

    def init(self):
        """
        Initialization and determination of the fastq file extension
        :return:
        """
        return

    @staticmethod
    def get_name():
        """
        Get the name of the extractor.
        :return: the name of the extractor
        """
        return 'FASTQ'

    @staticmethod
    def get_report_data_file_id():
        """
        Get the report.data id of the extractor.
        :return: the report.data id
        """
        return 'fastq.extractor'

    def add_key_to_result_dict(self, key):
        """
        Adds a key to the result_dict dictionary with the module name as a prefix
        :param key: key suffix
        :return: result_dict entry (string)
        """
        return '{0}.{1}'.format(self.get_report_data_file_id(), key)

    @staticmethod
    def add_value_to_unwritten_key(result_dict, value):
        """
        Adds a key name to the list of dictionary entries that will not be in the report.data file
        :param result_dict:
        :param value: result_dict entry (string)
        :return:
        """
        return result_dict['unwritten.keys'].append(value)

    def extract(self, result_dict):
        """
        Extraction of different information from the fastq file
        :param result_dict: dictionary which gathers all the extracted
        information that will be reported in the report.data file
        :return: result_dict
        """
        result_dict[self.add_key_to_result_dict('source')] = self.fastq_source

        if self.is_barcode:
            self._read_fastq_without_barcode()
            self._read_fastq_barcoded()
            result_dict.update(self.global_dico)

            for selected_barcode in self.barcode_selection:
                key = self.add_key_to_result_dict('length.') + selected_barcode
                self.add_value_to_unwritten_key(result_dict, key)

        else:
            self._read_fastq_without_barcode()
            result_dict.update(self.global_dico)
            self.add_value_to_unwritten_key(result_dict, '')

    @staticmethod
    def graph_generation(result_dict):
        """
        Generation of graphs for fastq files
        :return:
        """
        return []

    def clean(self, result_dict):
        """
        Removing dictionary entries that will not be kept in the report.data file
        :param result_dict: dictionary which gathers all the extracted
        information that will be reported in the report.data file
        :return:
        """

        keys = ['nucleotide.counter']
        key_list = []
        for key in keys:
            key_list.append(self.add_key_to_result_dict(key))
        result_dict['unwritten.keys'].extend(key_list)

    def _fastq_metrics(self):
        """
        Determination of different metrics
        :return:
        total_nucs_template: number of total nucleotides
        self.global_length_array: sequence length contained in the fastq file,
        barcode_length_array: sequence length for each barcode sample,
        template_nucleotide_counter: counting of the nucleotide present in each barcode
        """
        barcode_length_array = []
        template_nucleotide_counter = Counter()
        read_count = 0

        with self._open_compressed_file(self.fastq_file, self.fastq_file_extension) as input_file:
            with io.TextIOWrapper(input_file, encoding='utf-8') as fastq_file:
                entry_line = 0
                for line in fastq_file:
                    line = line.strip()
                    entry_line += 1

                    if entry_line == 2:
                        read_count += 1
                        template_nucleotide_counter.update(line)
                        self.global_length_array.append(len(line))

                        if self.is_barcode:
                            barcode_length_array.append(len(line))

                    if entry_line == 4:
                        entry_line = 0

        total_nucs_template = sum(template_nucleotide_counter.values())
        return (total_nucs_template, self.global_length_array, barcode_length_array,
                template_nucleotide_counter, read_count)

    @staticmethod
    def _open_compressed_file(file_path, file_extension):
        """
        Open a compressed file or not
        :param file_path: file path
        :param file_extension: file compressed file format in a string
        """

        if file_extension == 'bz2':
            return bz2.BZ2File(file_path, 'rb')

        elif file_extension == 'gz':
            return gzip.open(file_path, 'rb')

        else:
            return open(file_path, 'rb')

    def _barcoded_fastq_informations(self, selected_barcode=''):
        """
        Get different information about fastq files
        :param selected_barcode: barcode selection taken from the samplesheet file
        """
        total_nucs_template, self.global_length_array, barcode_length_array, template_nucleotide_counter, read_count = \
            self._fastq_metrics()

        series_read_size = pd.Series(barcode_length_array)
        selected_barcode_fastq_size_statistics = pd.Series.describe(series_read_size)
        self.global_dico[self.add_key_to_result_dict(selected_barcode + '.total.nucleotide')] = total_nucs_template

        for nucleotide, count in template_nucleotide_counter.items():
            self.global_dico[self.add_key_to_result_dict(selected_barcode + '.total.nucleotide.') + nucleotide] = \
                template_nucleotide_counter[nucleotide]

        for key in dict(selected_barcode_fastq_size_statistics).keys():
            self.global_dico[self.add_key_to_result_dict(selected_barcode + '.length.') + key] = \
                dict(selected_barcode_fastq_size_statistics)[key]

    def _read_fastq_barcoded(self):
        """
        Get information about the barcoded fastq sequence
        """
        self.init()
        if os.path.isfile(self.fastq_source):
            self.fastq_file = self.fastq_source
            for selected_barcode in self.barcode_selection:
                self._barcoded_fastq_informations(selected_barcode)

        if self.fastq_file_extension == 'bz2':
            for bz2_fastq_file in glob.glob("{}/*.bz2".format(self.fastq_source)):
                for selected_barcode in self.barcode_selection:
                    if selected_barcode in bz2_fastq_file:
                        self.fastq_file = bz2_fastq_file
                        self._barcoded_fastq_informations(selected_barcode)
                    else:
                        self.fastq_file = bz2_fastq_file
                        self._barcoded_fastq_informations('unclassified')

        elif self.fastq_file_extension == 'gz':
            for gz_fastq_file in glob.glob("{}/*.gz".format(self.fastq_source)):
                for selected_barcode in self.barcode_selection:
                    if selected_barcode in gz_fastq_file:
                        self.fastq_file = gz_fastq_file
                        self._barcoded_fastq_informations(selected_barcode)

        elif self.fastq_file_extension == 'fq':
            for fastq_file in glob.glob("{}/*.fq".format(self.fastq_source)):
                for selected_barcode in self.barcode_selection:
                    if selected_barcode in fastq_file:
                        self.fastq_file = fastq_file
                        self._barcoded_fastq_informations(selected_barcode)
        else:
            for fastq_files in glob.glob("{}/*.fastq".format(self.fastq_source)):
                for selected_barcode in self.barcode_selection:
                    if selected_barcode in fastq_files:
                        self.fastq_file = open(fastq_files, 'r')
                        self._barcoded_fastq_informations(selected_barcode)

    def _read_fastq_without_barcode(self):
        """
        Get information about the fastq sequence not barcoded
        """

        self.global_dico[self.add_key_to_result_dict('nucleotide.count')] = 0
        self.global_dico[self.add_key_to_result_dict('read.count')] = 0
        self.global_dico[self.add_key_to_result_dict('nucleotide.counter')] = Counter()

        self.init()
        if os.path.isfile(self.fastq_source):
            self.fastq_file = self.fastq_source
            self._fastq_without_barcode_information()

        elif self.fastq_file_extension == 'bz2':
            for bz2_fastq_file in glob.glob("{}/*.bz2".format(self.fastq_source)):
                self.fastq_file = bz2_fastq_file
                self._fastq_without_barcode_information()

        elif self.fastq_file_extension == 'gz':
            for bz2_fastq_file in glob.glob("{}/*.gz".format(self.fastq_source)):
                self.fastq_file = bz2_fastq_file
                self._fastq_without_barcode_information()

        elif self.fastq_file_extension == 'fq':
            for fastq_files in glob.glob("{}/*.fq".format(self.fastq_source)):
                self.fastq_file = fastq_files
                self._fastq_without_barcode_information()

        else:
            for fastq_files in glob.glob("{}/*.fastq".format(self.fastq_source)):
                self.fastq_file = fastq_files
                self._fastq_without_barcode_information()

        # Fill the result_dict dictionary
        self.global_dico[self.add_key_to_result_dict('mean.nucleotide.count.per.read')] = \
            self.global_dico[self.add_key_to_result_dict('nucleotide.count')] / \
            self.global_dico[self.add_key_to_result_dict('read.count')]

        for nucleotide, count in self.global_dico[self.add_key_to_result_dict('nucleotide.counter')].items():
            self.global_dico[self.add_key_to_result_dict('mean.nucleotide.count.') + nucleotide + '.per.read'] = \
                float(self.global_dico[self.add_key_to_result_dict('nucleotide.counter')][nucleotide]) / \
                float(self.global_dico[self.add_key_to_result_dict('read.count')])

            self.global_dico[self.add_key_to_result_dict('mean.nucleotide.ratio.') + nucleotide + '.per.read'] = \
                float(self.global_dico[self.add_key_to_result_dict('mean.nucleotide.count.') +
                                       nucleotide + '.per.read']) / \
                float(self.global_dico[self.add_key_to_result_dict('mean.nucleotide.count.per.read')])

            self.global_dico[self.add_key_to_result_dict('mean.nucleotide.frequency.') + nucleotide + '.per.read'] = \
                float(self.global_dico[self.add_key_to_result_dict('mean.nucleotide.ratio.')
                                       + nucleotide + '.per.read']) * 100

    def _fastq_without_barcode_information(self):
        """
        Get different information about all fastq files
        """
        total_nucs_template, self.global_length_array, _, template_nucleotide_counter, read_count = \
            self._fastq_metrics()

        # Fill the result_dict dictionary
        self.global_dico[self.add_key_to_result_dict('read.count')] += read_count
        self.global_dico[self.add_key_to_result_dict('nucleotide.count')] += total_nucs_template
        self.global_dico[self.add_key_to_result_dict('nucleotide.counter')] += template_nucleotide_counter
