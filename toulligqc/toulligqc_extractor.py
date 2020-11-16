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


from toulligqc import fastq_extractor
from toulligqc import fast5_extractor
from toulligqc import onedsquare_sequencing_summary_extractor
from toulligqc import sequencing_summary_extractor
from toulligqc import sequencing_telemetry_extractor


class ToulligqcExtractor:
    """
    Extraction of the extractor selected list.
    param config_dictionary: (configuration.py)
    param result_dict: dictionary which gathers all the extracted
    information that will be reported in the report.data file
    :return: extractors list : list of modules (extractors)
    that are called by the user from the command line.
    It depends on the parse_args function (toulligqc.py)
    that fills the config_dictionary.
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

    def check_conf(self):
        """Configuration checking"""
        return True, ""

    def init(self):
        """
        Initialisation
        :return:
        """
        return

    def extract(config_dictionary, extractors_list, result_dict):
        """
        Extraction of the different details about the config dictionary
        to create the extractors list that will be stored in the result_dict dictionary
        :param config_dictionary: configuration.py source
        :param extractors_list: list of modules that will be executed
        :param result_dict: result_dict dictionary
        :return: extractors_list
        """

        # Initialize the list of extractors
        result_dict['toulligqc.info.extractors'] = []

        if 'sequencing_telemetry_source' in config_dictionary and config_dictionary['sequencing_telemetry_source']:
            extractors_list.append(sequencing_telemetry_extractor.SequencingTelemetryExtractor(config_dictionary))
            result_dict['toulligqc.info.extractors'].append("sequencing.telemetry.extractor")

        if 'fast5_source' in config_dictionary and config_dictionary['fast5_source']:
            extractors_list.append(fast5_extractor.Fast5Extractor(config_dictionary))
            result_dict['toulligqc.info.extractors'].append("fast5.extractor")

        if 'sequencing_summary_1dsqr_source' in config_dictionary and config_dictionary['sequencing_summary_1dsqr_source']:
            extractors_list.append(onedsquare_sequencing_summary_extractor.
                                   OneDSquareSequencingSummaryExtractor(config_dictionary))
            result_dict['toulligqc.info.extractors'].append("sequencing.summary.1dsqr.extractor")
        else:
            extractors_list.append(sequencing_summary_extractor.
                                   SequencingSummaryExtractor(config_dictionary))
            result_dict['toulligqc.info.extractors'].append("sequencing.summary.extractor")
