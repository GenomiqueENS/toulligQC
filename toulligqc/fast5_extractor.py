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

import h5py
import glob
import sys
import os
import tarfile
import shutil
import tempfile
import dateutil


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
        self.result_directory = config_dictionary['result_directory']
        self.report_name = config_dictionary['report_name']
        self.fast5_file_extension = ''
        self.fast5_file = ''
        self.get_report_data_file_id()


    def check_conf(self):
        """
        Configuration checking
        :return:
        """

        if os.path.isdir(self.fast5_source):
            self.fast5_file_extension = 'fast5_directory'

        elif self.fast5_source.endswith('.tar.gz'):
            self.fast5_file_extension = 'tar.gz'

        elif self.fast5_source.endswith('.fast5'):
            self.fast5_file_extension = 'fast5'

        elif self.fast5_source.endswith('.tar.bz2'):
            self.fast5_file_extension = 'tar.bz2'

        else:
            return False, 'The fast5 extension is not supported (fast5, tar.bz2 or tar.gz format)'

        if self.fast5_file_extension != 'fast5_directory' and not os.path.isfile(self.fast5_source):
            return False, "The Fast5 source does not exists: " + self.fast5_source

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
        result_dict['sequencing.telemetry.extractor.source'] = self.fast5_source
        result_dict['sequencing.telemetry.extractor.flowcell.id'] = self._get_fast5_items(h5py_file, 'flow_cell_id')
        result_dict['sequencing.telemetry.extractor.minknow.version'] = self._get_fast5_items(h5py_file, 'version')
        result_dict['sequencing.telemetry.extractor.hostname'] = self._get_fast5_items(h5py_file, 'hostname')

        result_dict['sequencing.telemetry.extractor.operating.system'] = \
            self._get_fast5_items(h5py_file, 'operating_system')
        result_dict['sequencing.telemetry.extractor.device.id'] = self._get_fast5_items(h5py_file, 'device_id')

        result_dict['sequencing.telemetry.extractor.protocol.run.id'] = \
            self._get_fast5_items(h5py_file, 'protocol_run_id')
        result_dict['sequencing.telemetry.extractor.sample.id'] = self._get_fast5_items(h5py_file, 'sample_id')

        run_date = self._get_fast5_items(h5py_file, 'exp_start_time')
        exp_start_time = dateutil.parser.parse(run_date)
        result_dict['sequencing.telemetry.extractor.exp.start.time'] = exp_start_time.strftime("%x %X %Z")


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

    def _fast5_tar_bz2_extraction(self, tar_bz2_file, result_directory):
        """
        Extraction of the FAST5 file stored in tar_bz2 format
        :param tar_bz2_file: tar bz2 file containing the set of the raw FAST5 files
        :param result_directory:dictionary which gathers all the extracted
        information that will be reported in the report.data file
        :return: a FAST5 file
        """
        tar_bz2 = tarfile.open(tar_bz2_file, 'r:bz2')
        while True:
            member = tar_bz2.next()
            if member.name.endswith('.fast5'):
                tar_bz2.extract(member, path=result_directory)
                break
        return member.name

    def _fast5_tar_gz_extraction(self, tar_gz_file, result_directory):
        """
        Extraction of a FAST5 file stored in tar_gz format
        :param tar_gz_file: tar gz file containing the set of the raw FAST5 files
        :param result_directory: dictionary which gathers all the extracted
        information that will be reported in the report.data file
        :return: a FAST5 file
        """
        tar_gz = tarfile.open(self, tar_gz_file, 'r:gz')
        while True:
            member = tar_gz.next()
            if member.name.endswith('.fast5'):
                tar_gz.extract(member, path=result_directory)
                break
        return member.name

    def _read_fast5(self):
        """
        Extraction of one fast5 file from the archive and stores
        it in a h5py object for next retrieving information
        :return: h5py_file: h5py file
        """
        self.temporary_directory = tempfile.mkdtemp(dir=self.result_directory)
        if self.fast5_file_extension == 'tar.bz2':
            tar_bz2_file = self.fast5_source
            self.fast5_file = \
                self.temporary_directory + '/' + self._fast5_tar_bz2_extraction(tar_bz2_file, self.temporary_directory)

        elif self.fast5_file_extension == 'tar.gz':
            tar_gz_file = self.fast5_source
            self.fast5_file = \
                self.temporary_directory + '/' + self._fast5_tar_gz_extraction(tar_gz_file, self.temporary_directory)

        elif self.fast5_file_extension == 'fast5' or self.fast5_file_extension == '.fast5':
            self.fast5_file = self.fast5_source

        elif self.fast5_file_extension == 'fast5_directory':

            if glob.glob(self.fast5_source+'/*.fast5'):
                self.fast5_file = self.fast5_source+os.listdir(self.fast5_source)[0]

            elif glob.glob(self.fast5_source + '/*.tar.bz2'):
                tar_bz2_file = self.fast5_source+self.report_name+'.tar.bz2'
                self.fast5_file = \
                    self.temporary_directory + '/' + self._fast5_tar_bz2_extraction(tar_bz2_file,
                                                                                    self.temporary_directory)
            elif glob.glob(self.fast5_source + '/*.tar.gz'):
                tar_gz_file = self.fast5_source+self.report_name + '.tar.gz'
                self.fast5_file = self.temporary_directory + '/' + self._fast5_tar_gz_extraction(tar_gz_file,
                                                                                                 self.result_directory)
        else:
            err_msg = 'There is a problem with the fast5 file or the tar file'
            sys.exit(err_msg)
        h5py_file = h5py.File(self.fast5_file)

        return h5py_file

    def _get_fast5_items(self, h5py_file, params):
        """
        Global function to extract run information stores in h5py format
        :param h5py_file: fast5 file store in a h5py object
        :param params:  required h5py attributes
        :return: h5py value, for example flow_cell_id : FAE22827
        """
        if [key for key in h5py_file['/'].keys()][0] == 'Raw':
            tracking_id_items = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
            tracking_id_dict = {key: value.decode('utf-8') for key, value in tracking_id_items}

        else:
            tracking_id_items = list(h5py_file[[key for key in h5py_file['/'].keys()][0]+'/tracking_id'].attrs.items())
            tracking_id_dict = {key: value.decode('utf-8') for key, value in tracking_id_items}

        return tracking_id_dict[params]
