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
from toulligqc import albacore_1dsqr_sequencing_summary_extractor
from toulligqc import albacore_sequencing_summary_extractor
from toulligqc import albacore_pipeline_log_extractor


class ToulligqcExtractor:
    """
    Extraction of the extractor selected list
    :return: extractors list : list of modules (extractors) that are called by the user from the command line.
    It depends on the parse_args function (toulligqc.py) that fills the config_dictionary.
    """

    def __init__(self, config_dictionary):
        self._get_name()
        self._config_dictionary = config_dictionary

    def __getitem__(self, item):
        return self._config_dictionary[item]

    @staticmethod
    def _get_name():
        """
        Get the name of the extractor.
        :return: the name of the extractor
        """
        return 'Toulligqc extractors'

    def init(self):
        """
        Initialisation
        :return:
        """
        return

    def check_conf(self):
        """Configuration checking"""
        return

    def extract(config_dictionary, extractors_list, result_dict):
        """
        Extraction of the different details about the config dictionary
        to create the extractors list that will be stored in the result_dict dictionary
        :param config_dictionary:
        :param extractors_list:
        :param result_dict:
        :return: extractors_list
        """
        extractors_list.append(fast5_extractor.Fast5Extractor(config_dictionary))
        result_dict['toulligqc.info.extractors'] = ["fast5.extractor"]

        if 'albacore_pipeline_source' in config_dictionary and config_dictionary['albacore_pipeline_source']:
            extractors_list.append(albacore_pipeline_log_extractor.AlbacorePipelineLogExtractor(config_dictionary))
            result_dict['toulligqc.info.extractors'].append("albacore.log.extractor")

        if 'fastq_source' in config_dictionary and config_dictionary['fastq_source']:
            extractors_list.append(fastq_extractor.FastqExtractor(config_dictionary))
            result_dict['toulligqc.info.extractors'].append("fastq.extractor")

        if 'albacore_1dsqr_summary_source' in config_dictionary and config_dictionary['albacore_1dsqr_summary_source']:
            extractors_list.append(albacore_1dsqr_sequencing_summary_extractor.
                                   Albacore1DsqrSequencingSummaryExtractor(config_dictionary))
            result_dict['toulligqc.info.extractors'].append("albacore.1dsqr.stats.extractor")
        else:
            extractors_list.append(albacore_sequencing_summary_extractor.
                                   AlbacoreSequencingSummaryExtractor(config_dictionary))
            result_dict['toulligqc.info.extractors'].append("albacore.stats.extractor")
