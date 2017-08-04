import os
import glob
from matplotlib import pyplot as plt
from collections import Counter
import pandas as pd
import bz2
import seaborn as sns
import io
import parser
import numpy as np

class fastq():
    def __init__(self, pdf, result_directory, fastq_directory, dico_extension):
        self.global_dico = {}
        self.pdf = pdf
        self.global_length_array = []
        self.run_name, _, self.is_barcode, _ = parser.get_args()
        self.fastq_file = ''
        self.result_directory = result_directory
        # A voir pour le run_name si il faut l'indiquer dans le fichier config ou pas.
        self.fastq_directory = fastq_directory+self.run_name
        self.image_directory = self.result_directory+'images/'
        self.statistic_directory = self.result_directory+'statistics/'
        self.selection_global = []
        self.dico_extension = dico_extension

    def get_fastq_configuration(self):

        os.makedirs(self.image_directory)
        os.makedirs(self.statistic_directory)

    def fastq_metrics(self, file = ''):
        counter = 0
        barcode_length_array = []
        variable = ''
        if self.dico_extension['fastq_file_extension'] == 'bz2':
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
        """
        Get the fastq sequence
        """
        self.get_fastq_configuration()
        if glob.glob("{}/*.bz2".format(self.fastq_directory)):
            for bz2_fastq_file in glob.glob("{}/*.bz2".format(self.fastq_directory)):
                for selected_barcode in selection:
                    if selected_barcode in bz2_fastq_file:
                        self.fastq_file = bz2_fastq_file
                        self.barcoded_fastq_informations(selected_barcode)

        else:
            for fastq_files in glob.glob("{}/*.fastq".format(self.fastq_directory)):
                for selected_barcode in selection:
                    if selected_barcode in fastq_files:
                        self.fastq_file = open(fastq_files, 'r')
                        self.barcoded_fastq_informations(selected_barcode)

        return self.global_length_array, self.global_dico

    def get_fastq_without_barcode(self):
        """
        Gets the fastq sequence
        """
        self.get_fastq_configuration()

        if glob.glob("{}/*.bz2".format(self.fastq_directory)):
            for bz2_fastq_file in glob.glob("{}/*.bz2".format(self.fastq_directory)):
                self.fastq_file = bz2_fastq_file
                total_nucs_template, self.global_length_array, _, template_nucleotide_counter = self.fastq_metrics()
        else:
            for fastq_files in glob.glob("{}/*.fastq".format(self.fastq_directory)):
                self.fastq_file = open(fastq_files, 'r')
                total_nucs_template, self.global_length_array, _, template_nucleotide_counter = self.fastq_metrics()
        return total_nucs_template, self.global_length_array, _, template_nucleotide_counter


