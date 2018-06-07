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
from toulligqc import graph_generator



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


    def extract(self, result_dict):
        '''
        Extraction of the different informations about the fast5 files
        :param result_dict:
        :return: result_dict
       '''

        with open(self.pipeline_file, 'r') as pipeline_file:
            self.pipeline_dict['Fast5_submitted'] = 0
            self.pipeline_dict['Fast5_failed_to_load_key'] = 0
            self.pipeline_dict['Fast5_failed_count'] = 0
            self.pipeline_dict['Fast5_processed'] = 0

            for line in pipeline_file:
                if re.compile("(version)\s(\d+\.)(\d+\.)(\d)").search(line):
                    self.pipeline_dict['albacore_version'] = re.compile("\s(\d+\.)(\d+\.)(\d)").search(line).group(0)

                if re.compile("(SQK)\-([A-Z]{3})([0-9]{3})").search(line):
                    self.pipeline_dict['kit_version'] = re.compile("(SQK)\-([A-Z]{3})([0-9]{3})").search(line).group(0)

                if re.compile("(FLO)\-([A-Z]{3})([0-9]{3})").search(line):
                    self.pipeline_dict['flowcell_version'] = re.compile("(FLO)\-([A-Z]{3})([0-9]{3})").search(line).group(0)

                if re.compile('(sequence)\_(length)\_(template)').search(line):
                    self.pipeline_dict['Fast5_failed_to_load_key'] += 1

                if re.compile('(ERROR)\s(inserting)\s(read)').search(line):
                    self.pipeline_dict['Fast5_failed_count'] += 1

                if re.compile('(Finished)').search(line):
                    self.pipeline_dict['Fast5_processed'] += 1

                if re.compile('(Submitting)').search(line):
                    self.pipeline_dict['Fast5_submitted'] += 1

        pipeline_file.close()

        result_dict['albacore_version'] = self.pipeline_dict['albacore_version']
        result_dict['kit_version'] = self.pipeline_dict['kit_version']
        result_dict['flowcell_version'] = self.pipeline_dict['flowcell_version']
        result_dict['Fast5_failed_to_load_key'] = self.pipeline_dict['Fast5_failed_to_load_key']
        result_dict['Fast5_failed_count'] = self.pipeline_dict['Fast5_failed_count']
        result_dict['Fast5_processed'] = self.pipeline_dict['Fast5_processed']
        result_dict['Fast5_submitted'] = self.pipeline_dict['Fast5_submitted']

    def graph_generation(self):
        '''
        Graph generaiton
        :return:
        '''
        images_directory = self.result_directory + '/images'
        images = []
        images.append(graph_generator.log_count_histogram(self.pipeline_dict, 'About Albacore log', self.my_dpi,
                                                           images_directory,
                                                           "Number of reads submitted (Fast 5 in blue), proccesed (Fast 5 in blue) and with Error load key (Fast 5 in blue) or Error inserting file (Fast 5 in blue)."))
        return images

    def clean(self):
        '''
        Cleaning
        :return:
        '''
        return
