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

        result_dict['albacore_version'] = "Unknown"
        result_dict['kit_version'] = "Unknown"
        result_dict['flowcell_version'] = "Unknown"
        self.pipeline_dict['fast5_submitted'] = 0
        self.pipeline_dict['fast5_failed_to_load_key'] = 0
        self.pipeline_dict['fast5_failed_count'] = 0
        self.pipeline_dict['fast5_processed'] = 0

        with open(self.pipeline_file, 'r') as pipeline_file:


            for line in pipeline_file:
                if re.compile("(version)\s(\d+\.)(\d+\.)(\d)").search(line):
                    self.pipeline_dict['albacore_version'] = re.compile("\s(\d+\.)(\d+\.)(\d)").search(line).group(0)
                    result_dict['albacore_version'] = self.pipeline_dict['albacore_version']

                if re.compile("(SQK)\-([A-Z]{3})([0-9]{3})").search(line):
                    self.pipeline_dict['kit_version'] = re.compile("(SQK)\-([A-Z]{3})([0-9]{3})").search(line).group(0)
                    result_dict['kit_version'] = self.pipeline_dict['kit_version']

                if re.compile("(FLO)\-([A-Z]{3})([0-9]{3})").search(line):
                    self.pipeline_dict['flowcell_version'] = re.compile("(FLO)\-([A-Z]{3})([0-9]{3})").search(line).group(0)
                    result_dict['flowcell_version'] = self.pipeline_dict['flowcell_version']

                if re.compile("(key\:)\s('(sequence)\_(length)\_(template)')").search(line):
                    self.pipeline_dict['fast5_failed_to_load_key'] += 1

                if re.compile('(ERROR)\s(inserting)\s(read)').search(line):
                    self.pipeline_dict['fast5_failed_count'] += 1

                if re.compile('(Finished)').search(line):
                    self.pipeline_dict['fast5_processed'] += 1

                if re.compile('(Submitting)').search(line):
                    self.pipeline_dict['fast5_submitted'] += 1

        pipeline_file.close()

        result_dict['raw_fast5'] = self.pipeline_dict['fast5_submitted']
        result_dict['fast5_failed_to_load_key'] = self.pipeline_dict['fast5_failed_to_load_key']
        result_dict['fast5_failed_count'] = self.pipeline_dict['fast5_failed_count']
        result_dict['fast5_processed'] = self.pipeline_dict['fast5_processed']

        result_dict['raw_fast5_no_processed']=self.pipeline_dict['fast5_submitted'] - self.pipeline_dict['fast5_processed']
        result_dict['basecalled_error_count']= result_dict['raw_fast5_no_processed'] + result_dict['fast5_failed_to_load_key'] + result_dict['fast5_failed_count']



    def graph_generation(self,result_dict):
        '''
        Graph generaiton
        :return:
        '''
        return []

    def clean(self):
        '''
        Cleaning
        :return:
        '''
        return
