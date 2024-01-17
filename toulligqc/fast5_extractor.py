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
# First author: Lionel Ferrato-Berberian
# Maintainer: Bérengère Laffay
# Since version 0.1

# Extraction of the information about the FAST5 files

import os
import shutil
import sys
import tarfile
import tempfile
import h5py

from toulligqc.common import find_file_in_directory
from toulligqc.common import set_result_dict_value

class Fast5Extractor:
    """
    Extraction of different information from a FAST5 file
    param fast5_source: FAST5 file directory
    param result_directory: dictionary which gathers all the extracted
    information that will be reported in the report.data file
    param fast5_file_extension: extension used for the storage of the set of FAST5 files if there's one
    param report_name: report name
    return: a tuple containing the information about a FAST5 file
    """

    def __init__(self, config_dictionary):
        self.config_file_dictionary = config_dictionary
        self.fast5_source = config_dictionary['fast5_source']
        self.file_to_process = None
        self.report_name = config_dictionary['report_name']
        self.fast5_file_extension = ''
        self.fast5_file = ''
        self.get_report_data_file_id()

    def check_conf(self):
        """
        Configuration checking
        :return:
        """

        if not os.path.exists(self.fast5_source):
            return False, 'The input file or directory for Fast5 file does not exists: ' + self.fast5_source

        if os.path.isdir(self.fast5_source):
            file_found = find_file_in_directory(self.fast5_source, "fast5")
            if file_found is None:
                return False, 'No Fast5 file found in directory: ' + self.fast5_source
            self.file_to_process = file_found
        else:
            self.file_to_process = self.fast5_source

        if self.file_to_process.endswith('.tar'):
            self.fast5_file_extension = 'tar'

        elif self.file_to_process.endswith('.tar.gz'):
            self.fast5_file_extension = 'tar.gz'

        elif self.file_to_process.endswith('.tar.bz2'):
            self.fast5_file_extension = 'tar.bz2'

        elif self.file_to_process.endswith('.fast5'):
            self.fast5_file_extension = 'fast5'

        else:
            return (False, 'The file extension for FAST5 input is not supported '
                           '(only .fast5, .tar, .tar.gz or .tar.bz2 are supported): ' + self.fast5_source)

        return True, ""


    def init(self):
        """
        Determination of the fast5 file extension
        """
        return


    @staticmethod
    def get_name():
        """
        Get the name of the extractor.
        :return: the name of the extractor
        """
        return 'Fast5'


    @staticmethod
    def get_report_data_file_id():
        """
        Get the report.data id of the extractor
        :return: the report.data id
        """
        return 'fast5.extractor'


    def extract(self, result_dict):
        """
        Extraction of the different information about the fast5 files
        :param result_dict: Dictionary which gathers all the extracted
        information that will be reported in the report.data file
        :return: result_dict filled
        """
        h5py_file = self._read_fast5()
        tracking_id_dict = self._get_fast5_items(h5py_file, 'tracking_id')

        if len(tracking_id_dict) == 0:
            return

        prefix = 'sequencing.telemetry.extractor'
        result_dict[prefix + '.source'] = self.fast5_source
        set_result_dict_value(result_dict, prefix + '.flowcell.id', tracking_id_dict, 'flow_cell_id')
        set_result_dict_value(result_dict, prefix + '.minknow.version', tracking_id_dict, 'version')
        set_result_dict_value(result_dict, prefix + '.hostname', tracking_id_dict, 'hostname')
        set_result_dict_value(result_dict, prefix + '.operating.system', tracking_id_dict, 'operating_system')
        set_result_dict_value(result_dict, prefix + '.run.id', tracking_id_dict, 'run_id')
        set_result_dict_value(result_dict, prefix + '.protocol.run.id', tracking_id_dict, 'protocol_run_id')
        set_result_dict_value(result_dict, prefix + '.protocol.group.id', tracking_id_dict, 'protocol_group_id')
        set_result_dict_value(result_dict, prefix + '.sample.id', tracking_id_dict, 'sample_id')
        set_result_dict_value(result_dict, prefix + '.exp.start.time', tracking_id_dict, 'exp_start_time')
        set_result_dict_value(result_dict, prefix + '.device.id', tracking_id_dict, 'device_id')
        set_result_dict_value(result_dict, prefix + '.device.type', tracking_id_dict, 'device_type')
        set_result_dict_value(result_dict, prefix + '.distribution.version', tracking_id_dict, 'distribution_version')
        set_result_dict_value(result_dict, prefix + '.flow.cell.product.code', tracking_id_dict,
                               'flow_cell_product_code')

        context_tags_dict = self._get_fast5_items(h5py_file, 'context_tags')
        if len(context_tags_dict) != 0:
            set_result_dict_value(result_dict, prefix + '.selected.speed.bases.per.second', context_tags_dict, 'selected_speed_bases_per_second')
            set_result_dict_value(result_dict, prefix + '.sample.frequency', context_tags_dict, 'sample_frequency')
            set_result_dict_value(result_dict, prefix + '.sequencing.kit.version', context_tags_dict, 'sequencing_kit')


    def graph_generation(self, result_dict):
        """
        Graph generation
        :return: nothing
        """
        return []


    def clean(self, result_dict):
        """
        Deleting the temporary fast5 file extracted from the tar archive if used
        and removing dictionary entries that will not be kept in the report.data file
        :param result_dict: dictionary which gathers all the extracted
        information that will be reported in the report.data file
        :return:
        """
        if self.temporary_directory:
            shutil.rmtree(self.temporary_directory, ignore_errors=True)


    def _fast5_tar_extraction(self, tar_file, extension, output_directory):
        """
        Extraction of a FAST5 file stored in a tar file,
        :param tar_file: tar file containing the set of the raw FAST5 files
        :param extension of the file
        :param output_directory: dictionary which gathers all the extracted
        information that will be reported in the report.data file
        :return: a FAST5 file
        """

        if extension == 'tar.gz':
            tar_mode = 'r:gz'
        elif extension == 'tar.bz2':
            tar_mode = 'r:bz2'
        else:
            tar_mode = 'r'

        tf = tarfile.open(name=tar_file, mode=tar_mode)
        while True:
            member = tf.next()
            if member.name.endswith('.fast5'):
                tf.extract(member, path=output_directory)
                break
        return output_directory + '/' + member.name


    def _read_fast5(self):
        """
        Extraction of one fast5 file from the archive and stores
        it in a h5py object for next retrieving information
        :return: h5py_file: h5py file
        """
        self.temporary_directory = tempfile.mkdtemp()
        if self.fast5_file_extension == 'tar' or \
                self.fast5_file_extension == 'tar.gz' or \
                self.fast5_file_extension == 'tar.bz2':
            self.fast5_file = self._fast5_tar_extraction(self.file_to_process, self.fast5_file_extension,
                                                         self.temporary_directory)
        elif self.fast5_file_extension == 'fast5' or self.fast5_file_extension == '.fast5':
            self.fast5_file = self.file_to_process
        else:
            err_msg = 'There is a problem with the fast5 file or the tar file'
            sys.exit(err_msg)
        h5py_file = h5py.File(self.fast5_file)

        return h5py_file


    def _get_fast5_items(self, h5py_file, group):
        """
        Global function to extract run information stores in h5py format
        :param h5py_file: fast5 file store in a h5py object
        :param key:  required h5py attributes
        :return: h5py value, for example flow_cell_id : FAE22827
        """

        for k in h5py_file['/'].keys():
            new_group = '/' + k + '/' + group
            if new_group in h5py_file:
                tracking_id_items = list(h5py_file[new_group].attrs.items())
                tracking_id_dict = {key: value.decode('utf-8') for key, value in tracking_id_items}
                return tracking_id_dict

        return {}
