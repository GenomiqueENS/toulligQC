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

import os
import glob
from collections import Counter
import pandas as pd
import bz2
import io
import gzip
import sys
import re

class fastq_extractor():
    def __init__(self, config_dictionary):
        self.config_dictionary = config_dictionary
        self.global_dico = {}
        self.global_length_array = []
        self.run_name = config_dictionary['run_name']
        self.is_barcode = config_dictionary['barcoding']
        self.fastq_file = ''
        self.result_directory = config_dictionary['result_directory']
        self.barcode_selection = config_dictionary['barcode_selection']
        print(self.barcode_selection)
        self.fastq_source = config_dictionary['fastq_source']
       # if os.path.isdir(self.fastq_source):
        #    if self.run_name in self.fastq_source:
         #       self.fastq_source = config_dictionary['fastq_source']
          #  else:
           #     self.fastq_source = config_dictionary['fastq_source']+self.run_name

        self.image_directory = self.result_directory+'images/'
        self.statistic_directory = self.result_directory+'statistics/'
        self.selection_global = []

        self.fastq_file_extension = ''

    def init(self):
        if os.path.isdir(self.fastq_source):
            if glob.glob(self.fastq_source + '/*.fastq') or glob.glob(self.fastq_source + self.run_name+ '/*.fastq'):
                self.fastq_file_extension = 'fastq'
            
            elif glob.glob(self.fastq_source + '/*.fq') or glob.glob(self.fastq_source + self.run_name+ '/*.fq'):
                self.fastq_file_extension = 'fq'
            
            elif glob.glob(self.fastq_source + '/*.gz') or glob.glob(self.fastq_source + self.run_name+ '/*.gz'):
                self.fastq_file_extension = 'gz'
            
            elif glob.glob(self.fastq_source + '/*.bz2') or glob.glob(self.fastq_source + self.run_name+ '/*.bz2'):
                self.fastq_file_extension = 'bz2'

            elif glob.glob(self.fastq_source + self.run_name+ '/*.fastq'):
                self.fastq_source = self.fastq_source+self.run_name+'/'
                self.fastq_file_extension = 'fastq'

            elif glob.glob(self.fastq_source + self.run_name+ '/*.fq'):
                self.fastq_source = self.fastq_source+self.run_name+'/'
                self.fastq_file_extension = 'fq'

            elif glob.glob(self.fastq_source + self.run_name+ '/*.gz'):
                self.fastq_source = self.fastq_source+self.run_name+'/'
                self.fastq_file_extension = 'gz'

            elif glob.glob(self.fastq_source + self.run_name+ '/*.bz2'):
                self.fastq_source = self.fastq_source+self.run_name+'/'
                self.fastq_file_extension = 'bz2'

            else:
                print('The fastq source extension is not supported (fast5, tar.bz2 or tar.gz format)')
                sys.exit(0)


        elif self.fastq_source.endswith('.fastq'):
            self.fastq_file_extension = 'fastq'

        elif self.fastq_source.endswith('.fq'):
            self.fastq_file_extension= 'fq'

        elif self.fastq_source.endswith('.bz2') or self.fastq_source.endswith('.gz') or self.fastq_source.endswith('.zip'):
            pattern = '\.(gz|bz2|zip)$'
            if re.search(pattern, self.fastq_source):
                match = re.search(pattern, self.fastq_source)
                self.fastq_file_extension = match.groups()[0]
        else:
            print('The fastq source extension is not supported (fast5, bz2 or gz format)')
            sys.exit(0)

    def check_conf(self):
        return


    def extract(self,result_dict):
        if self.is_barcode:
            self._read_fastq_barcoded()

            result_dict.update(self.global_dico)
        else:
            self._read_fastq_without_barcode()
            result_dict.update(self.global_dico)
        
        return result_dict

    def graph_generation(self):
        return []

    def clean(self):
        return

    def _get_fastq_configuration(self):

        os.makedirs(self.image_directory)
        os.makedirs(self.statistic_directory)

    def _fastq_metrics(self):
        '''
        Determination of different metrics
        :return: total nucleotides,
        sequence length contained in the fastq file,
        sequence length for each barcode sample,
        counting of the nucleotide present in each barcode
        '''
        counter = 0
        barcode_length_array = []
        sequence = ''
        if self.fastq_file_extension == 'bz2':
            with bz2.BZ2File(self.fastq_file, 'rb') as inputo:
                with io.TextIOWrapper(inputo, encoding='utf-8') as bz2_fastq_file:
                    for line in bz2_fastq_file:
                        counter += 1

                        if counter == 2:
                            sequence += line.strip()
                            self.global_length_array.append(len(line))

                            if self.is_barcode:
                                barcode_length_array.append(len(line))

                        if counter == 4:
                            counter = 0

        elif self.fastq_file_extension == 'gz':
            with gzip.open(self.fastq_file, 'rb') as input_file:
                with io.TextIOWrapper(input_file, encoding='utf-8') as bz2_fastq_file:
                    for line in bz2_fastq_file:
                        counter += 1

                        if counter == 2:
                            sequence += line.strip()
                            self.global_length_array.append(len(line))

                            if self.is_barcode:
                                barcode_length_array.append(len(line))

                        if counter == 4:
                            counter = 0

        else:
            with open(self.fastq_file, 'r') as fastq_file:
                for line in fastq_file:
                    counter += 1

                    if counter == 2:
                        sequence += line.strip()
                        self.global_length_array.append(len(line))

                    if self.is_barcode:
                        barcode_length_array.append(len(line))

                    if counter == 4:
                        counter = 0

        template_nucleotide_counter = Counter(sequence)
        total_nucs_template = len(sequence)
        return total_nucs_template, self.global_length_array, barcode_length_array, template_nucleotide_counter

    def _barcoded_fastq_informations(self, selected_barcode= ''):
        '''
        Get different information about fastq files
        :param selected_barcode: barcode selection
        '''
        total_nucs_template, self.global_length_array, barcode_length_array, template_nucleotide_counter = self._fastq_metrics()
        series_read_size = pd.Series(barcode_length_array)
        selected_barcode_fastq_size_statistics = pd.Series.describe(series_read_size)
        self.global_dico['nucleotide_count_'+selected_barcode] = template_nucleotide_counter
        self.global_dico['total_nucleotide_'+selected_barcode] = total_nucs_template
        self.global_dico['fastq_length_'+selected_barcode] = selected_barcode_fastq_size_statistics


    def _read_fastq_barcoded(self):
        '''
        Get informations about the barcoded fastq sequence
        :param selection: barcode selection
        :return: length of all sequences of a fastq barcoded file,
        dictionary containing different information such as the nucleotide counting and other
        '''
        self._get_fastq_configuration()
        self.init()
        if os.path.isfile(self.fastq_source):
            self.fastq_file = self.fastq_source
            for selected_barcode in self.barcode_selection:
                self._barcoded_fastq_informations(selected_barcode)

        if self.fastq_file_extension == 'bz2':
            for bz2_fastq_file in glob.glob("{}/*.bz2".format(self.fastq_source)):
                for selected_barcode in self.barcode_selection:
                    if selected_barcode in bz2_fastq_file:
                        self.fastq_file = bz2_fastq_file
                        self._barcoded_fastq_informations(selected_barcode)

        elif self.fastq_file_extension == 'gz':
            for gz_fastq_file in glob.glob("{}/*.gz".format(self.fastq_source)):
                for selected_barcode in self.barcode_selection:
                    if selected_barcode in gz_fastq_file:
                        self.fastq_file = gz_fastq_file
                        self._barcoded_fastq_informations(selected_barcode)

        elif self.fastq_file_extension == 'fq':
            for fastq_file in glob.glob("{}/*.fq".format(self.fastq_source)):
                for selected_barcode in self.barcode_selection:
                    if selected_barcode in fastq_file:
                        self.fastq_file = fastq_file
                        self._barcoded_fastq_informations(selected_barcode)
        else:
            for fastq_files in glob.glob("{}/*.fastq".format(self.fastq_source)):
                for selected_barcode in self.barcode_selection:
                    if selected_barcode in fastq_files:
                        self.fastq_file = open(fastq_files, 'r')
                        self._barcoded_fastq_informations(selected_barcode)

        #self.global_dico['global_length_array'] = self.global_length_array

                #return self.global_length_array, self.global_dico


    def _read_fastq_without_barcode(self):
        '''
        Gets informations about the fastq sequence not barcoded
        :return: total nucleotides,
        sequence length contained in the fastq file,
        counting of the different nucleotide for the entire sample
        '''
        self._get_fastq_configuration()
        self.init()
        if os.path.isfile(self.fastq_source):
            self.fastq_file = self.fastq_source
            total_nucs_template, self.global_length_array, _, template_nucleotide_counter = self._fastq_metrics()

        elif self.fastq_file_extension == 'bz2':
            for bz2_fastq_file in glob.glob("{}/*.bz2".format(self.fastq_source)):
                self.fastq_file = bz2_fastq_file
                total_nucs_template, self.global_length_array, _, template_nucleotide_counter = self._fastq_metrics()

        elif self.fastq_file_extension == 'gz':
            for bz2_fastq_file in glob.glob("{}/*.gz".format(self.fastq_source)):
                self.fastq_file = bz2_fastq_file
                total_nucs_template, self.global_length_array, _, template_nucleotide_counter = self._fastq_metrics()

        elif self.fastq_file_extension == 'fq':
            for fastq_files in glob.glob("{}/*.fq".format(self.fastq_source)):
                self.fastq_file = fastq_files
                total_nucs_template, self.global_length_array, _, template_nucleotide_counter = self._fastq_metrics()

        else:
            for fastq_files in glob.glob("{}/*.fastq".format(self.fastq_source)):
                self.fastq_file = fastq_files
                total_nucs_template, self.global_length_array, _, template_nucleotide_counter = self._fastq_metrics()

        fastq_length_statistics = pd.Series(self.global_length_array).describe()
        #self.global_dico['global_length_array'] = fastq_length_statistics
        self.global_dico['total_nucleotide'] = total_nucs_template
        self.global_dico['nucleotide_count'] = template_nucleotide_counter

#return total_nucs_template, self.global_length_array, _, template_nucleotide_counter


