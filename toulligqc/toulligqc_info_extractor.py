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

import datetime
import os
import platform as pf
import sys
import tempfile as tp


class ToulligqcInfoExtractor:
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

    def __init__(self, config_dictionary, extractors_list):
        self._config_dictionary = config_dictionary
        self._extractors_list = extractors_list
        self._debug = False

        if 'debug' in config_dictionary and config_dictionary['debug'].lower() == 'true':
            self._debug = True

    @staticmethod
    def get_name() -> str:
        """
        Get the name of the extractor.
        :return: the name of the extractor
        """
        return 'Toulligqc info'

    @staticmethod
    def get_report_data_file_id():
        """
        Get the report.data id of the extractor
        :return: the report.data id
        """
        return 'toulligqc.info.extractor'

    def check_conf(self):
        """Configuration checking"""
        return True, ""

    def init(self):
        """
        Initialisation
        :return:
        """
        return

    def extract(self, result_dict):
        """
        Extraction of the different details about the config dictionary
        to create the extractors list that will be stored in the result_dict dictionary
        :param config_dictionary: configuration.py source
        :param extractors_list: list of modules that will be executed
        :param result_dict: result_dict dictionary
        :return: extractors_list

        """

        result_dict['unwritten.keys'] = ['unwritten.keys']

        # Add ToulligQC info
        self._toulligqc_info(result_dict)

        # Add system and python information
        if self._debug:
            self._system_and_python_info(result_dict)

        # Add QC info
        self._qc_info(result_dict)

        # Add the list of used extractors
        result_dict['toulligqc.info.extractors'] = []
        for e in self._extractors_list:
            result_dict['toulligqc.info.extractors'].append(e.get_report_data_file_id())

    def graph_generation(self, result_dict):
        """
        Graph generation
        :return: nothing
        """
        return []

    def clean(self, result_dict):
        return

    @staticmethod
    def _system_and_python_info(result_dict):
        """
        Initialization of the result_dict with the OS parameters and the environment variables
        :param config_dictionary: details from command user line
        :param result_dict: Dictionary which gathers all the extracted
        information that will be reported in the report.data file
        :return: result_dict dictionary and extractors list
        """
        result_dict['toulligqc.info.system.hostname'] = os.uname()[1]
        result_dict['toulligqc.info.system.username'] = os.environ.get('USERNAME')
        result_dict['toulligqc.info.system.user.home'] = os.environ['HOME']
        result_dict['toulligqc.info.system.temporary.directory'] = tp.gettempdir()
        result_dict['toulligqc.info.system.operating.system'] = pf.processor()

        # Environment variables
        for name, value in os.environ.items():
            result_dict['toulligqc.info.system.env.' + name] = value

        # Python info
        result_dict['toulligqc.info.python.version'] = pf.python_version()
        result_dict['toulligqc.info.python.implementation'] = pf.python_implementation()

        # Python dependencies versions
        for name, module in sorted(sys.modules.items()):
            if hasattr(module, '__version__'):
                result_dict['toulligqc.info.python.dependancy.' + name + '.version'] = module.__version__
            elif hasattr(module, 'VERSION'):
                result_dict['toulligqc.info.python.dependancy.' + name + '.version'] = module.VERSION

        return result_dict

    def _toulligqc_info(self, result_dict):

        result_dict['toulligqc.info.version'] = self._config_dictionary['app.version']
        result_dict['toulligqc.info.start.time'] = datetime.datetime.now().astimezone().replace(
            microsecond=0).isoformat()
        result_dict['toulligqc.info.report.name'] = self._config_dictionary['report_name']
        result_dict['toulligqc.info.executable.path'] = sys.argv[0]
        result_dict['toulligqc.info.command.line'] = sys.argv

    def _qc_info(self, result_dict):

        result_dict['toulligqc.info.html.report.path'] = self._config_dictionary.get('html_report_path', 'Undefined')
        result_dict['toulligqc.info.data.report.path'] = self._config_dictionary.get('data_report_path', 'Undefined')
        result_dict['toulligqc.info.image.directory'] = self._config_dictionary.get('images_directory', 'Undefined')
        result_dict['toulligqc.info.barcode.option'] = "False"

        if self._config_dictionary['barcoding'].lower() == 'true':
            result_dict['toulligqc.info.barcode.option'] = "True"
            result_dict['toulligqc.info.barcode.selection'] = self._config_dictionary['barcode_selection']
