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
        self.telemetry_file = ''

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

    def init(self):
        """
        Determination of the telemetry file
        """
        if os.path.isdir(self.telemetry_source):
            self.telemetry_file = self.telemetry_source + "/sequencing_telemetry.js"
        else:
            telemetry_file = self.telemetry_source

    def check_conf(self):
        """
        Configuration checking
        :return: nothing
        """
        return

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

        with open(self.telemetry_source, 'r') as f:
            array = json.load(f)

            result_dict['fast5.extractor.source'] = self.telemetry_source
            result_dict['fast5.extractor.flowcell.id'] = array[0]['tracking_id']['flow_cell_id']
            result_dict['fast5.extractor.minknow.version'] = array[0]['tracking_id']['version']
            result_dict['fast5.extractor.hostname'] = array[0]['tracking_id']['hostname']
            result_dict['fast5.extractor.operating.system'] = array[0]['tracking_id']['operating_system']
            result_dict['fast5.extractor.device.id'] = array[0]['tracking_id']['hostname']
            result_dict['fast5.extractor.protocol.run.id'] = array[0]['tracking_id']['protocol_run_id']
            result_dict['fast5.extractor.sample.id'] = array[0]['tracking_id']['sample_id']
            result_dict['fast5.extractor.exp.start.time'] = array[0]['tracking_id']['exp_start_time']

            # TODO Add software and protocol name
            #result_dict['albacore.log.extractor.albacore.version'] = array[0]['software']['analysis']
            #result_dict['albacore.log.extractor.albacore.version'] = array[0]['software']['name']
            result_dict['albacore.log.extractor.albacore.version'] = array[0]['software']['version']

            if 'albacore_opts' in array[0]:
                result_dict['albacore.log.extractor.kit.version'] = array[0]['albacore_opts']['kit']
                result_dict['albacore.log.extractor.flowcell.version'] = array[0]['albacore_opts']['flowcell']
            else:
                result_dict['albacore.log.extractor.kit.version'] = array[0]['opts']['kit']
                result_dict['albacore.log.extractor.flowcell.version'] = array[0]['opts']['flowcell']

        return
