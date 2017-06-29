import parser
import os
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import glob
from matplotlib import pyplot as plt
from collections import Counter
import pandas as pd
import bz2
import seaborn as sns
import numpy as np

class fastq():
    def __init__(self, pdf):
        self.pdf = pdf
        self.global_length_array = []
        self.run_name, self.selected_file_list, _, self.is_barcode = parser.get_args()
        self.path_dico = parser.file_path_initialization()
        self.fastq_file = 'fastq_sequence.txt'
        self.result_directory = self.path_dico['result_directory']
        self.fastq_directory = self.path_dico['fastq_directory']+self.run_name
        self.image_directory = self.result_directory+'images/'
        self.statistic_directory = self.result_directory+'statistics/'
        self.selection_global = []
        self.dico = {}

    def get_fastq_configuration(self):

     #   if os.path.isfile(self.fastq_file):
      #      open(self.fastq_file, 'w').close()

       # if not os.path.exists(self.image_directory):
        os.makedirs(self.image_directory)

        #if os.path.exists(self.image_directory):
        #    for file in os.listdir(self.image_directory):
        #        os.remove(self.image_directory+file)

       # if not os.path.exists(self.statistic_directory):
        os.makedirs(self.statistic_directory)

        #if os.path.exists(self.statistic_directory):
        #    for file in os.listdir(self.statistic_directory):
        #        os.remove(self.statistic_directory+file)

    def bz2_decompression(self, bz2_fastq_file):
        print(bz2_fastq_file)
        uncompressedData = bz2.BZ2File(bz2_fastq_file, 'rb').read()
        uncomp = uncompressedData.decode('utf-8')
        selected_file = open(self.fastq_file, 'w')
        selected_file.write(uncomp)
        selected_file.close()

    def fastq_metrics(self):
        total_nucs_template = 0
        barcode_length_array = []
        template_nucleotide_counter = Counter()
        with open(self.fastq_file) as in_handle:
            for title, seq, qual in FastqGeneralIterator(in_handle):
                if self.is_barcode:
                    barcode_length_array.append(len(seq))
                self.global_length_array.append(len(seq))

                for template_nucleotide in seq:
                    template_nucleotide_counter[template_nucleotide] += 1
                    total_nucs_template += 1
        return total_nucs_template, self.global_length_array, barcode_length_array, template_nucleotide_counter

    def barcoded_fastq_informations(self, selected_barcode):
        total_nucs_template, self.global_length_array, barcode_length_array, template_nucleotide_counter = self.fastq_metrics()
        completeName = os.path.join(self.statistic_directory, selected_barcode)
        barcode_file = open(completeName, 'w')
        series_read_size = pd.Series(barcode_length_array)
        selected_barcode_fastq_size_statistics = pd.Series.describe(series_read_size)

        for index, value in selected_barcode_fastq_size_statistics.iteritems():
            barcode_file.write("Read.fastq.length.{}={}\n".format(index, value))

        for nucleotide, count in template_nucleotide_counter.items():
            barcode_file.write("nucleotide.{}.template={}\n".format(nucleotide, count))
            if nucleotide == 'total':
                continue
            calcul = float(count) / float(total_nucs_template)
            barcode_file.write("nucleotide.{}.proportion={}\n".format(nucleotide, calcul))
           # plt.boxplot(barcode_length_array, showfliers=False)
           # plt.title('Read length boxplot for {}'.format(selected_barcode))
           # plt.xticks
           # plt.savefig(self.image_directory+'image_{}.png'.format(selected_barcode))
           # plt.close()
        barcode_file.close()
        self.dico[selected_barcode] = barcode_length_array
        #self.selection_global.append(barcode_length_array)

    def get_fastq_barcoded(self, selection):
        """
        Get the fastq sequence
        """
        self.get_fastq_configuration()
        if glob.glob("{}/*.bz2".format(self.fastq_directory)) != []:
            for bz2_fastq_file in glob.glob("{}/*.bz2".format(self.fastq_directory)):
                for selected_barcode in selection:
                    if selected_barcode in bz2_fastq_file:
                        self.bz2_decompression(bz2_fastq_file)
                        self.barcoded_fastq_informations(selected_barcode)
            mpl_fig = plt.figure()
            ax = mpl_fig.add_subplot(111)
            #sns.boxplot(self.selection_global, showfliers=False)
            #ax.boxplot(self.selection_global, showfliers=False)
            df = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in self.dico.items()]))
            sns.boxplot(data = df,showfliers=False)
            plt.xlabel('Barcodes')
            plt.ylabel('Read size(in pb)')
            plt.title('Read size for each barcode')
           # plt.xticks(list(range(1, len(selection)+1)),selection)
            plt.savefig(self.image_directory+'barcode_total.png')
            self.pdf.savefig()
            plt.close()


        else:
            for fastq_files in glob.glob("{}/*.fastq".format(self.fastq_directory)):
                for selected_barcode in selection:
                    if selected_barcode in fastq_files:
                        self.barcoded_fastq_informations(selected_barcode)

        return self.global_length_array

    def get_fastq_without_barcode(self):
        """
        Gets the fastq sequence
        """
        self.get_fastq_configuration()
        if not self.selected_file_list:
            if glob.glob("{}/*.bz2".format(self.fastq_directory)) != []:
                for bz2_fastq_file in glob.glob("{}/*.bz2".format(self.fastq_directory)):
                    self.bz2_decompression(bz2_fastq_file)
                    total_nucs_template, self.global_length_array, _, template_nucleotide_counter = self.fastq_metrics()

            else:
                for fastq_files in glob.glob("{}/*.fastq".format(self.fastq_directory)):
                    self.barcoded_fastq_informations(selected_barcode)
        else:
            for selected_file in self.selected_file_list:
                selected_file = self.fastq_directory + '/' + selected_file
                self.bz2_decompression(selected_file)
                total_nucs_template, self.global_length_array, _, template_nucleotide_counter = self.fastq_metrics()

        return template_nucleotide_counter, total_nucs_template, self.global_length_array

