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
# First author: Laurent Jourdren
# Maintainer: Bérengère Laffay
# Since version 1.1

# Extraction of run information from the sequencing_telemetry.js file

import json
import os.path


class SequencingTelemetryExtractor:

    def __init__(self, config_dictionary):
        self.config_file_dictionary = config_dictionary
        self.telemetry_source = config_dictionary['sequencing_telemetry_source']

        if os.path.isdir(self.telemetry_source):
            self.telemetry_file = self.telemetry_source + "/sequencing_telemetry.js"
        else:
            self.telemetry_file = self.telemetry_source

    def check_conf(self):
        """
        Configuration checking
        :return: nothing
        """

        if not os.path.isfile(self.telemetry_file):
            return False, "Telemetry file does not exists: " + self.telemetry_file

        return True, ""

    def init(self):
        """
        Determination of the telemetry file
        """
        return

    @staticmethod
    def get_name():
        """
        Get the name of the extractor.
        :return: the name of the extractor
        """
        return 'Sequencing telemetry'

    @staticmethod
    def get_report_data_file_id():
        """
        Get the report.data id of the extractor.
        :return: the report.data id
        """
        return 'sequencing.telemetry.extractor'

    def graph_generation(self, result_dict):
        """
        Graph generation
        :return: nothing
        """
        return []

    def clean(self, result_dict):
        """
        :return: nothing
        """
        return

    def extract(self, result_dict):
        """
        Extraction of data from the sequencing_telemetry.js file
        :param result_dict: Dictionary which gathers all the extracted
        information that will be reported in the report.data file
        :return: result_dict filled
        """

        with open(self.telemetry_file, 'r') as f:
            array = json.load(f)

            result_dict[self.get_report_data_file_id() + '.source'] = self.telemetry_file
            result_dict[self.get_report_data_file_id() + '.flowcell.id'] = array[0]['tracking_id']['flow_cell_id']
            result_dict[self.get_report_data_file_id() + '.minknow.version'] = array[0]['tracking_id']['version']
            result_dict[self.get_report_data_file_id() + '.hostname'] = array[0]['tracking_id']['hostname']
            result_dict[self.get_report_data_file_id() + '.operating.system'] = array[0]['tracking_id']['operating_system']
            result_dict[self.get_report_data_file_id() + '.protocol.run.id'] = array[0]['tracking_id']['protocol_run_id']
            result_dict[self.get_report_data_file_id() + '.sample.id'] = array[0]['tracking_id']['sample_id']
            result_dict[self.get_report_data_file_id() + '.exp.start.time'] = array[0]['tracking_id']['exp_start_time']
            result_dict[self.get_report_data_file_id() + '.device.id'] = array[0]['tracking_id']['device_id']
            result_dict[self.get_report_data_file_id() + '.software.name'] = array[0]['software']['name']
            result_dict[self.get_report_data_file_id() + '.software.version'] = array[0]['software']['version']
            result_dict[self.get_report_data_file_id() + '.software.analysis'] = array[0]['software']['analysis']

            if 'albacore_opts' in array[0]:
                result_dict[self.get_report_data_file_id() + '.kit.version'] = array[0]['albacore_opts']['kit']
                result_dict[self.get_report_data_file_id() + '.flowcell.version'] = array[0]['albacore_opts']['flowcell']
                result_dict[self.get_report_data_file_id() + '.model.file'] = array[0]['context_tags']['local_bc_temp_model']
            else:
                result_dict[self.get_report_data_file_id() + '.kit.version'] = array[0]['opts']['kit']
                result_dict[self.get_report_data_file_id() + '.flowcell.version'] = array[0]['opts']['flowcell']
                result_dict[self.get_report_data_file_id() + '.model.file'] = array[0]['opts']['model_file']
                if 'device_type' in array[0]['tracking_id']:
                    result_dict[self.get_report_data_file_id() + '.device.type'] = array[0]['tracking_id']['device_type']

        return
