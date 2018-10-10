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

from toulligqc import fastq_extractor
from toulligqc import fast5_extractor
from toulligqc import albacore_1dsqr_stats_generator
from toulligqc import albacore_stats_extractor
from toulligqc import pipeline_log_extractor


class toulligqc_extractor():
    '''
    Extraction of the extractor selected list
    :param config_dict: toulligqc config dictionnary
    :param extractors_list: extractors list
    :return: extractors list
    '''

    def __init__(self, config_dictionary):
        self.get_name()

    def get_name(self):
        '''
        Get the name of the extractor.
        :return: the name of the extractor
        '''
        return 'Toulligqc extractors'

    def init(self):
        '''
        Initialisation
        :return:
        '''
        return

    def check_conf(self):
        '''Configuration checking'''
        return

    def extract(config_dictionary, extractors, result_dict):
        '''
        Extraction of the different informations about the config dictionary
        :param extractors_list:
        :param result_dict:
        :return: extractors_list
       '''
        extractors.append(fast5_extractor.fast5_extractor(config_dictionary))
        result_dict['toulligqc.info.extractors'] = ["fast5.extractor"]

        if 'albacore_pipeline_source' in config_dictionary and config_dictionary['albacore_pipeline_source']:
            extractors.append(pipeline_log_extractor.albacore_log_extractor(config_dictionary))
            result_dict['toulligqc.info.extractors'].append("albacore.log.extractor")

        if 'fastq_source' in config_dictionary and config_dictionary['fastq_source']:
            extractors.append(fastq_extractor.fastq_extractor(config_dictionary))
            result_dict['toulligqc.info.extractors'].append("fastq.extractor")

        if 'albacore_1dsqr_summary_source' in config_dictionary and config_dictionary['albacore_1dsqr_summary_source']:
            extractors.append(albacore_1dsqr_stats_generator.albacore_1dsqr_stats_extractor(config_dictionary))
            result_dict['toulligqc.info.extractors'].append("albacore.1dsqr.stats.extractor")
        else:
            extractors.append(albacore_stats_extractor.albacore_stats_extractor(config_dictionary))
            result_dict['toulligqc.info.extractors'].append("albacore.stats.extractor")



