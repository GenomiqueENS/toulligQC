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

import h5py
import glob
import sys
import os
import tarfile
import shutil
import tempfile
import pathlib

class fast5_extractor():
    '''
       Extraction of different informations from a FAST5 file
       :param fast5_source: FAST5 file directory
       :param result_directory: result directory
       :param fast5_file_extension: extension used for the storage of the set of FAST5 files if there's one
       :param run_name: run name
       :return: a tuple containing the informations about a FAST5 file
       '''

    def __init__(self, config_dictionary):
        self.config_file_dictionary =config_dictionary
        self.fast5_source = config_dictionary['fast5_source']
        self.result_directory = config_dictionary['result_directory']
        self.run_name = config_dictionary['run_name']
        self.fast5_file_extension = ''
        self.fast5_file = ''

    def init(self):

        if os.path.isdir(self.fast5_source):
           self.fast5_file_extension = 'fast5_directory'

        elif self.fast5_source.endswith('.tar.gz'):
            self.fast5_file_extension = 'tar.gz'

        elif self.fast5_source.endswith('.fast5'):
            self.fast5_file_extension = 'fast5'

        elif self.fast5_source.endswith('.tar.bz2'):
            self.fast5_file_extension = 'tar.bz2'

        else:
            print('The fast5 extension is not supported (fast5, tar.bz2 or tar.gz format)')
            sys.exit(0)

    def extract(self, result_dict):
        h5py_file = self._read_fast5()
        result_dict['flowcell_id'] = self._get_flowcell_id(h5py_file)
        result_dict['minknow_version'] = self._get_minknow_version(h5py_file)
        result_dict['hostname'] = self._get_hostname(h5py_file)
        result_dict['minion_run_id'] = self._get_minion_run_id(h5py_file)
        result_dict['protocol_run_id'] = self._get_protocol_run_id(h5py_file)
        return result_dict

    def check_conf(self):
        return

    def graph_generation(self):
        return []

    def clean(self):
        if self.temporary_directory:
            shutil.rmtree(self.temporary_directory, ignore_errors=True)
        else:
            return

    def _fast5_tar_bz2_extraction(self, tar_bz2_file, result_directory):
        '''
        Extraction of a FAST5 file from a set of FAST5 files
        :param tar_bz2_file: tar bz2 file containing the set of the raw FAST5 files
        :param result_directory: result directory
        :return: a FAST5 file
        '''
        tar_bz2 = tarfile.open(tar_bz2_file, 'r:bz2')
        while True:
            member = tar_bz2.next()
            if member.name.endswith('.fast5'):
                tar_bz2.extract(member, path=result_directory)
                break
        return member.name

    def _fast5_tar_gz_extraction(self, tar_gz_file, result_directory):
        '''
        Extraction of a FAST5 file from a set of FAST5 files
        :param tar_gz_file: tar gz file containing the set of the raw FAST5 files
        :param result_directory: result directory
        :return: a FAST5 file
        '''
        tar_gz = tarfile.open(self, tar_gz_file, 'r:gz')
        while True:
            member = tar_gz.next()
            if member.name.endswith('.fast5'):
                tar_gz.extract(member, path=result_directory)
                break
        return member.name

    def _read_fast5(self):

        self.temporary_directory = tempfile.mkdtemp(dir=self.result_directory)
        if self.fast5_file_extension == 'tar.bz2':
            tar_bz2_file = self.fast5_source
            self.fast5_file = self.temporary_directory + '/' + self._fast5_tar_bz2_extraction(tar_bz2_file, self.temporary_directory)

        elif self.fast5_file_extension == 'tar.gz':
            tar_gz_file = self.fast5_source
            self.fast5_file = self.temporary_directory + '/' + self._fast5_tar_gz_extraction(tar_gz_file, self.temporary_directory)

        elif self.fast5_file_extension == 'fast5' or self.fast5_file_extension == '.fast5':
            self.fast5_file = self.fast5_source

        elif self.fast5_file_extension == 'fast5_directory':

            if glob.glob(self.fast5_source+self.run_name+'/*.fast5'):
                self.fast5_file = self.fast5_source+self.run_name+'.fast5'

            elif glob.glob(self.fast5_source + '/*.tar.bz2'):
                tar_bz2_file = self.fast5_source+self.run_name+'.tar.bz2'
                self.fast5_file = self.temporary_directory + '/' + self._fast5_tar_bz2_extraction(tar_bz2_file, self.temporary_directory)
            elif glob.glob(self.fast5_source + '/*.tar.gz'):
                tar_gz_file = self.fast5_source+self.run_name+ '.tar.gz'
                self.fast5_file = self.temporary_directory + '/' + self._fast5_tar_gz_extraction(tar_gz_file, self.result_directory)



        else:
            print('There is a problem with the fast5 file or the tar file')
            sys.exit(0)
        h5py_file = h5py.File(self.fast5_file)

        return h5py_file

    def _get_minknow_version(self,h5py_file):
        """
        Get the Minknow version from fast5 file
        """
        version = list(h5py_file['/UniqueGlobalKey/tracking_id'].attrs.items())
        version_d = {key: value.decode('utf-8') for key, value in version}
        return version_d['version']

    def _get_flowcell_id(self,h5py_file):
        """
        Get the flowcell id from fast5 file
        """
        flowcell_id = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
        flowcell_id_dico = {key: value.decode('utf-7') for key, value in flowcell_id}
        return flowcell_id_dico['flow_cell_id']

    def _get_hostname(self,h5py_file):
        """
        Get the hostname from fast5 file
        """
        host_name = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
        host_name_dico = {key: value.decode('utf-8') for key, value in host_name}
        return host_name_dico['hostname']

    def _get_minion_run_id(self,h5py_file):
        """
        Get the number of Minion run
        """
        numMinION = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
        minion_run_id_dico = {key: value.decode('utf-8') for key, value in numMinION}
        return minion_run_id_dico['device_id']

    def _get_protocol_run_id(self, h5py_file):
        """
        Get the run id protocol from fast 5 file
        """
        protocol_run_id =  list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
        protocol_run_id_dico = {key: value.decode('utf-8') for key, value in protocol_run_id}
        return protocol_run_id_dico['protocol_run_id']



