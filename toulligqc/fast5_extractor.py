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
#

#Extraction of the informations about the FAST5 files

import h5py
import glob
import sys
import os
import tarfile
import shutil
import tempfile

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

    def get_name(self):
        '''
        Get the name of the extractor.
        :return: the name of the extractor
        '''
        return 'FAST5'


    def init(self):
        '''
        Determination of the fast5 file extension
        '''
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
        '''
        Extraction of the different informations about the fast5 files
        :param result_dict:
        :return: result_dict
       '''
        h5py_file = self._read_fast5()
        result_dict['flow_cell_id'] = self._get_fast5_items(h5py_file,'flow_cell_id')
        result_dict['minknow_version'] = self._get_fast5_items(h5py_file,'version')
        result_dict['hostname'] = self._get_fast5_items(h5py_file,'hostname')
        result_dict['minion_run_id'] = self._get_fast5_items(h5py_file,'device_id')
        result_dict['protocol_run_id'] = self._get_fast5_items(h5py_file,'protocol_run_id')
        result_dict['exp_start_time'] = self._get_fast5_items(h5py_file,'exp_start_time')
        result_dict['sample_id'] = self._get_fast5_items(h5py_file,'sample_id')

    def check_conf(self):
        '''
        Configuration checking
        :return:
        '''
        return

    def graph_generation(self):
        '''
        Graph generaiton
        :return:
        '''
        return []

    def clean(self):
        '''
        Deleting the temporary fast5 file extracted from the tar archive if used
        :return:
        '''
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
        '''
        Extraction of one fast5 file from the archive and stores it in a h5py object for next retrieving informations
        :return: h5py_file: h5py file
        '''
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

    def _get_fast5_items(self,h5py_file,params):
        '''
        Global function to exctract run informations stores in h5py format
        :param h5py_file: fast5 file store in a h5py object
        :param params:  required h5py attributes
        :return: h5py value, for example flow_cell_id : FAE22827
        '''
        tracking_id_items = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
        tracking_id_dict = {key: value.decode('utf-8') for key, value in tracking_id_items}
        return tracking_id_dict[params]



