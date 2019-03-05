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
# First author: Bérengère Laffay
# Maintainer: Bérengère Laffay
# Since version 0.10 

# Extraction of information from the pipeline.log file


import os
import re


class AlbacorePipelineLogExtractor:
    """
    Extraction of information from the pipeline.log file
    like the flowcell version.
    param pipeline_source: path to the pipeline.log file (config_dictionary)
    :return result_dict filled
    """

    def __init__(self, config_dictionary):
        self.config_file_dictionary = config_dictionary
        self.pipeline_source = config_dictionary['albacore_pipeline_source']
        self.result_directory = config_dictionary['result_directory']
        self.pipeline_file = ''
        self.my_dpi = int(config_dictionary['dpi'])
        self.pipeline_dict = {}
        self.get_report_data_file_id()
        self.dict = {}

        if os.path.isdir(self.pipeline_source):
            self.pipeline_file = self.pipeline_source + "/pipeline.log"
        else:
            self.pipeline_file = self.pipeline_source

    def check_conf(self):
        """
        Configuration checking
        :return:
        """

        if not os.path.isfile(self.pipeline_file):
            return False, "Pipeline log file does not exists: " + self.pipeline_file

        return True, ""

    def init(self):
        """
        Determination of the pipeline.log file extension
        """
        return

    @staticmethod
    def get_name():
        """
        Get the name of the extractor.
        :return: the name of the extractor
        """
        return 'Albacore pipeline log'

    @staticmethod
    def get_report_data_file_id():
        """
        Get the report.data id of the extractor.
        :return: the report.data id
        """
        return 'albacore.log.extractor'

    def add_key_to_result_dict(self, key):
        """
        Adding a key to the result_dict dictionary with the module name as a prefix
        :param key: key suffix
        :return: result_dict entry (string)
        """
        return '{0}.{1}'.format(self.get_report_data_file_id(), key)

    @staticmethod
    def _is_in_result_dict(result_dict, dict_key, default_value):
        """
        Global function to check for the presence of an entry in a dictionary
        and give it a default value.
        :param result_dict: result_dict dictionary
        :param dict_key: entry (string)
        :param default_value:
        :return:
        """
        if dict_key not in result_dict or not result_dict[dict_key]:
            result_dict[dict_key] = default_value
        return result_dict[dict_key]

    def extract(self, result_dict):
        """
        Extraction of the different information about the pipeline.log file
        :param result_dict: result_dict dictionary
        :return: result_dict filled
        """

        result_dict[self.add_key_to_result_dict('source')] = self.pipeline_source

        albacore_version_regex = re.compile("(version)\s(\d+\.)(\d+\.)(\d)")
        kit_version_regex = re.compile("(SQK)" + re.escape('-') + "([A-Z]{3})([0-9]{3})")
        flowcell_version_regex = re.compile("(FLO)" + re.escape('-') + "([A-Z]{3})([0-9]{3})")
        fail_to_load_regex = re.compile("(key)" + re.escape(':') + " \s('(sequence) " + re.escape('_')
                                      + "(length)" + re.escape('_') + "(template)')")
        error_regex = re.compile('(ERROR)\s(inserting)\s(read)')
        finished_regex = re.compile('(Finished)')
        submitting_regex = re.compile('(Submitting)')

        with open(self.pipeline_file, 'r') as pipeline_file:
            result_dict[self.add_key_to_result_dict('fast5.files.failed.to.load.key')] = 0
            result_dict[self.add_key_to_result_dict('fast5.files.failed.count')] = 0
            result_dict[self.add_key_to_result_dict('fast5.files.processed')] = 0
            result_dict[self.add_key_to_result_dict('fast5.files.submitted')] = 0

            for line in pipeline_file:
                if albacore_version_regex.search(line):
                    result_dict['sequencing.telemetry.extractor.software.version'] = \
                        re.compile("(\d+\.)(\d+\.)(\d)").search(line).group(0)
                    result_dict['sequencing.telemetry.extractor.software.name'] = "albacore-basecalling"

                if kit_version_regex.search(line):
                    result_dict['sequencing.telemetry.extractor.kit.version'] = \
                        re.compile("(SQK)" + re.escape('-') + "([A-Z]{3})([0-9]{3})").search(line).group(0)

                if flowcell_version_regex.search(line):
                    result_dict['sequencing.telemetry.extractor.flowcell.version'] = \
                        re.compile("(FLO)" + re.escape('-') + "([A-Z]{3})([0-9]{3})").search(line).group(0)

                if fail_to_load_regex.search(line):
                    result_dict[self.add_key_to_result_dict('fast5.files.failed.to.load.key')] += 1

                if error_regex.search(line):
                    result_dict[self.add_key_to_result_dict('fast5.files.failed.count')] += 1

                if finished_regex.search(line):
                    result_dict[self.add_key_to_result_dict('fast5.files.processed')] += 1

                if submitting_regex.search(line):
                    result_dict[self.add_key_to_result_dict('fast5.files.submitted')] += 1

        pipeline_file.close()

        # Count of the read number (Fast5 files) that have not been converted into Fastq
        result_dict[self.add_key_to_result_dict('raw.fast5.files.not.processed')] = \
            result_dict[self.add_key_to_result_dict('fast5.files.submitted')] \
            - result_dict[self.add_key_to_result_dict('fast5.files.processed')]

        result_dict[self.add_key_to_result_dict('fast5.files.basecalled.error.count')] = \
            result_dict[self.add_key_to_result_dict('raw.fast5.files.not.processed')] \
            + result_dict[self.add_key_to_result_dict('fast5.files.failed.to.load.key')] \
            + result_dict[self.add_key_to_result_dict('fast5.files.failed.count')]

        # Ration and frequency deduction

        result_dict[self.add_key_to_result_dict('fast5.files.ratio')] = \
            result_dict[self.add_key_to_result_dict('fast5.files.submitted')] / \
            result_dict[self.add_key_to_result_dict('fast5.files.submitted')]

        result_dict[self.add_key_to_result_dict('fast5.files.basecalled.error.ratio')] = \
            result_dict[self.add_key_to_result_dict('fast5.files.basecalled.error.count')] / \
            result_dict[self.add_key_to_result_dict('fast5.files.submitted')]

        result_dict[self.add_key_to_result_dict('fast5.files.frequency')] = \
            result_dict[self.add_key_to_result_dict('fast5.files.submitted')] / \
            result_dict[self.add_key_to_result_dict('fast5.files.submitted')]*100

        result_dict[self.add_key_to_result_dict('fast5.files.basecalled.error.frequency')] = \
            result_dict[self.add_key_to_result_dict('fast5.files.basecalled.error.count')] / \
            result_dict[self.add_key_to_result_dict('fast5.files.submitted')]*100

        result_dict[self.add_key_to_result_dict('flowcell.version')] = \
            self._is_in_result_dict(result_dict, self.add_key_to_result_dict('flowcell.version'), "unknown")
        result_dict[self.add_key_to_result_dict('kit.version')] = \
            self._is_in_result_dict(result_dict, self.add_key_to_result_dict('kit.version'), "unknown")
        result_dict[self.add_key_to_result_dict('albacore.version')] = \
            self._is_in_result_dict(result_dict, self.add_key_to_result_dict('albacore.version'), "unknown")

    @staticmethod
    def graph_generation(result_dict):
        """
        Graph generaiton
        :return: nothing
        """
        return []

    def clean(self, result_dict):
        """
        Removing dictionary entries that will not be kept in the report.data file
        :param result_dict: result_dict dictionary
        :return:
        """
        keys = []
        key_list = []
        for key in keys:
            key_list.extend(self.add_key_to_result_dict(key))
        result_dict['unwritten.keys'].extend(key_list)
