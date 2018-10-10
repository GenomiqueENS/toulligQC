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

# Extraction of statistics from sequencing_summary.txt file

import glob
import sys
import os
import tarfile
import shutil
import tempfile
import re  # python 3.5 package


class albacore_log_extractor():
    '''
    Extraction of informations from piepline.log file
    :param result_dict:
    :param config_dictionary:
    '''

    def __init__(self, config_dictionary):
        self.config_file_dictionary = config_dictionary
        self.pipeline_source = config_dictionary['albacore_pipeline_source']
        self.result_directory = config_dictionary['result_directory']
        self.pipeline_file = ''
        self.my_dpi = int(config_dictionary['dpi'])
        self.pipeline_dict = {}
        self.get_report_data_file_id()
        self.dict={}

    def check_conf(self):
        '''
        Configuration checking
        :return:
        '''
        return

    def init(self):
        '''
        Determination of the pipeline.log file extension
        '''
        if os.path.isdir(self.pipeline_source):
            self.pipeline_file = self.pipeline_source + "/pipeline.log"
        else:
            self.pipeline_file = self.pipeline_source

    def get_name(self):
        '''
        Get the name of the extractor.
        :return: the name of the extractor
        '''
        return 'ALBACORE PIPELINE LOG'

    def get_report_data_file_id(self):
        '''
        Get the report.data id of the extractor.
        :return: the report.data id
        '''
        return 'albacore.log.extractor'

    def add_key_to_result_dict(self, key):
        return '{0}.{1}'.format(self.get_report_data_file_id(), key)

    def _is_in_result_dict(self,result_dict, dict_key, default_value):
        if dict_key not in result_dict or not result_dict[dict_key]:
            result_dict[dict_key] = default_value
        return result_dict[dict_key]

    def extract(self, result_dict):
        '''
        Extraction of the different informations about the fast5 files
        :param result_dict:
        :return: result_dict
        '''

        result_dict[self.add_key_to_result_dict('source')] = self.pipeline_source

        with open(self.pipeline_file, 'r') as pipeline_file:
            result_dict[self.add_key_to_result_dict('fast5.files.failed.toload.key')] = 0
            result_dict[self.add_key_to_result_dict('fast5.files.failed.count')] = 0
            result_dict[self.add_key_to_result_dict('fast5.files.processed')] = 0
            result_dict[self.add_key_to_result_dict('fast5.files.submitted')] = 0

            for line in pipeline_file:
                if re.compile("(version)\s(\d+\.)(\d+\.)(\d)").search(line):
                    result_dict[self.add_key_to_result_dict('albacore.version')] = re.compile("(\d+\.)(\d+\.)(\d)").search(line).group(0)

                if re.compile("(SQK)\-([A-Z]{3})([0-9]{3})").search(line):
                    result_dict[self.add_key_to_result_dict('kit.version')] = re.compile("(SQK)\-([A-Z]{3})([0-9]{3})").search(line).group(0)

                if re.compile("(FLO)\-([A-Z]{3})([0-9]{3})").search(line):
                    result_dict[self.add_key_to_result_dict('flowcell.version')] = re.compile("(FLO)\-([A-Z]{3})([0-9]{3})").search(line).group(0)

                if re.compile("(key\:)\s('(sequence)\_(length)\_(template)')").search(line):
                    result_dict[self.add_key_to_result_dict('fast5.files.failed.toload.key')] += 1

                if re.compile('(ERROR)\s(inserting)\s(read)').search(line):
                    result_dict[self.add_key_to_result_dict('fast5.files.failed.count')] += 1

                if re.compile('(Finished)').search(line):
                    result_dict[self.add_key_to_result_dict('fast5.files.processed')] += 1

                if re.compile('(Submitting)').search(line):
                    result_dict[self.add_key_to_result_dict('fast5.files.submitted')] += 1

        pipeline_file.close()

        result_dict[self.add_key_to_result_dict('raw.fast5.files.not.processed')] = result_dict[self.add_key_to_result_dict('fast5.files.submitted')] - result_dict[self.add_key_to_result_dict('fast5.files.processed')]
        result_dict[self.add_key_to_result_dict('fast5.files.basecalled.error.count')] = result_dict[self.add_key_to_result_dict('raw.fast5.files.not.processed')] + result_dict[self.add_key_to_result_dict('fast5.files.failed.toload.key')] + result_dict[self.add_key_to_result_dict('fast5.files.failed.count')]
        result_dict[self.add_key_to_result_dict('fast5.files.ratio')] = result_dict[self.add_key_to_result_dict('fast5.files.submitted')]/result_dict[self.add_key_to_result_dict('fast5.files.submitted')]
        result_dict[self.add_key_to_result_dict('fast5.files.basecalled.error.ratio')] = result_dict[self.add_key_to_result_dict('fast5.files.basecalled.error.count')]/result_dict[self.add_key_to_result_dict('fast5.files.submitted')]
        result_dict[self.add_key_to_result_dict('fast5.files.frequency')] = result_dict[self.add_key_to_result_dict('fast5.files.submitted')]/result_dict[self.add_key_to_result_dict('fast5.files.submitted')]*100
        result_dict[self.add_key_to_result_dict('fast5.files.basecalled.error.frequency')] = result_dict[self.add_key_to_result_dict('fast5.files.basecalled.error.count')]/result_dict[self.add_key_to_result_dict('fast5.files.submitted')]*100

        result_dict[self.add_key_to_result_dict('flowcell.version')] = self._is_in_result_dict(result_dict,'albacore.log.extractor.flowcell.version', "unknown")
        result_dict[self.add_key_to_result_dict('kit.version')] = self._is_in_result_dict(result_dict, 'albacore.log.extractor.kit.version', "unknown")
        result_dict[self.add_key_to_result_dict('albacore.version')] = self._is_in_result_dict(result_dict, 'albacore.log.extractor.albacore.version', "unknown")

    def graph_generation(self, result_dict):
        '''
        Graph generaiton
        :return:
        '''
        return []

    def clean(self, result_dict):
        '''
        Cleaning
        :return:
        '''
        keys = []
        key_list = []
        for key in keys:
            key_list.extend(self.add_key_to_result_dict(key))
        result_dict['unwritten.keys'].extend(key_list)
