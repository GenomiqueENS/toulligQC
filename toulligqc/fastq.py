import os
import glob
from collections import Counter
import pandas as pd
import bz2
import io
import numpy as np
import gzip

class fastq():
    def __init__(self, result_directory, fastq_directory, run_name, is_barcode, fastq_file_extension):
        self.global_dico = {}
        self.global_length_array = []
        self.run_name,self.is_barcode = run_name, is_barcode
        self.fastq_file = ''
        self.result_directory = result_directory
        if os.path.isdir(fastq_directory):
            self.fastq_directory = fastq_directory+self.run_name

        self.image_directory = self.result_directory+'images/'
        self.statistic_directory = self.result_directory+'statistics/'
        self.selection_global = []
        self.fastq_file_extension = fastq_file_extension

    def get_fastq_configuration(self):

        os.makedirs(self.image_directory)
        os.makedirs(self.statistic_directory)

    def fastq_metrics(self):
        '''
        Determination of different metrics
        :return: total nucleotides,
        sequence length contained in the fastq file,
        sequence length for each barcode sample,
        counting of the nucleotide present in each barcode
        '''
        counter = 0
        barcode_length_array = []
        variable = ''
        if self.fastq_file_extension == 'bz2':
            print(self.fastq_file)
            with bz2.BZ2File(self.fastq_file, 'rb') as inputo:
                with io.TextIOWrapper(inputo, encoding='utf-8') as bz2_fastq_file:
                    for line in bz2_fastq_file:
                        counter += 1

                        if counter == 2:
                            variable += line.strip()
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
                            variable += line.strip()
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
                        variable += line.strip()
                        self.global_length_array.append(len(line))

                    if self.is_barcode:
                        barcode_length_array.append(len(line))

                    if counter == 4:
                        counter = 0

        template_nucleotide_counter = Counter(variable)
        total_nucs_template = len(variable)

        return total_nucs_template, self.global_length_array,barcode_length_array, template_nucleotide_counter

    def barcoded_fastq_informations(self, selected_barcode= ''):
        '''
        Get different information about fastq files
        :param selected_barcode: barcode selection
        '''
        total_nucs_template, self.global_length_array, barcode_length_array, template_nucleotide_counter = self.fastq_metrics()
        completeName = os.path.join(self.statistic_directory, selected_barcode)
        barcode_file = open(completeName, 'w')
        series_read_size = pd.Series(barcode_length_array)
        selected_barcode_fastq_size_statistics = pd.Series.describe(series_read_size)
        self.global_dico['fastq_size_'+selected_barcode] = barcode_length_array
        self.global_dico['nucleotide_count_'+selected_barcode] = template_nucleotide_counter
        self.global_dico['total_nucleotide_'+selected_barcode] = total_nucs_template
        for index, value in selected_barcode_fastq_size_statistics.iteritems():
            print(type(value))
            barcode_file.write("Read.fastq.length.{}={}\n".format(index, np.round(value, decimals=2)))
        for nucleotide, count in template_nucleotide_counter.items():
            barcode_file.write("nucleotide.{}.template={}\n".format(nucleotide, np.round(value, decimals=2)))
            if nucleotide == 'total':
                continue
            calcul = float(count) / float(total_nucs_template)
            calcul = calcul*100
            barcode_file.write("nucleotide.{}.proportion={}\n".format(nucleotide,  np.round(calcul, decimals=2)))
        barcode_file.close()
        #self.dico[selected_barcode] = barcode_length_array

    def get_fastq_barcoded(self, selection):
        '''
        Get informations about the barcoded fastq sequence
        :param selection: barcode selection
        :return: length of all sequences of a fastq barcoded file,
        dictionary containing different information such as the nucleotide counting and other
        '''
        self.get_fastq_configuration()

        if os.path.isfile(self.fastq_directory):
            self.fastq_file = self.fastq_directory
            for selected_barcode in selection:
                self.barcoded_fastq_informations(selected_barcode)

        if self.fastq_file_extension == 'bz2':
            for bz2_fastq_file in glob.glob("{}/*.bz2".format(self.fastq_directory)):
                for selected_barcode in selection:
                    if selected_barcode in bz2_fastq_file:
                        self.fastq_file = bz2_fastq_file
                        self.barcoded_fastq_informations(selected_barcode)

        elif self.fastq_file_extension == 'gz':
            for gz_fastq_file in glob.glob("{}/*.gz".format(self.fastq_directory)):
                for selected_barcode in selection:
                    if selected_barcode in gz_fastq_file:
                        self.fastq_file = gz_fastq_file
                        self.barcoded_fastq_informations(selected_barcode)

        elif self.fastq_file_extension == 'fq':
            for fastq_file in glob.glob("{}/*.fq".format(self.fastq_directory)):
                for selected_barcode in selection:
                    if selected_barcode in fastq_file:
                        self.fastq_file = fastq_file
                        self.barcoded_fastq_informations(selected_barcode)
        else:
            for fastq_files in glob.glob("{}/*.fastq".format(self.fastq_directory)):
                for selected_barcode in selection:
                    if selected_barcode in fastq_files:
                        self.fastq_file = open(fastq_files, 'r')
                        self.barcoded_fastq_informations(selected_barcode)

        return self.global_length_array, self.global_dico


    def get_fastq_without_barcode(self):
        '''
        Gets informations about the fastq sequence not barcoded
        :return: total nucleotides,
        sequence length contained in the fastq file,
        counting of the different nucleotide for the entire sample
        '''
        self.get_fastq_configuration()

        if os.path.isfile(self.fastq_directory):
            self.fastq_file = self.fastq_directory
            total_nucs_template, self.global_length_array, _, template_nucleotide_counter = self.fastq_metrics()

        elif self.fastq_file_extension == 'bz2':
            for bz2_fastq_file in glob.glob("{}/*.bz2".format(self.fastq_directory)):
                self.fastq_file = bz2_fastq_file
                total_nucs_template, self.global_length_array, _, template_nucleotide_counter = self.fastq_metrics()

        elif self.fastq_file_extension == 'gz':
            for bz2_fastq_file in glob.glob("{}/*.gz".format(self.fastq_directory)):
                self.fastq_file = bz2_fastq_file
                total_nucs_template, self.global_length_array, _, template_nucleotide_counter = self.fastq_metrics()

        elif self.fastq_file_extension == 'fq':
            for fastq_files in glob.glob("{}/*.fq".format(self.fastq_directory)):
                self.fastq_file = fastq_files
                total_nucs_template, self.global_length_array, _, template_nucleotide_counter = self.fastq_metrics()

        else:
            for fastq_files in glob.glob("{}/*.fastq".format(self.fastq_directory)):
                self.fastq_file = fastq_files
                total_nucs_template, self.global_length_array, _, template_nucleotide_counter = self.fastq_metrics()
        return total_nucs_template, self.global_length_array, _, template_nucleotide_counter


