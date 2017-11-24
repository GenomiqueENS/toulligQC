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

#Extraction of statistics from sequencing_summary.txt file

import pandas as pd
import sys
from toulligqc import graph_generator
import numpy as np
import re

class albacore_stats_extractor():
    '''
    Extraction of statistics from sequencing_summary.txt file and graph generation
    '''
    def __init__(self, config_dictionary):
        self.global_dictionnary = {}
        self.config_dictionary = config_dictionary
        self.albacore_log = pd.read_csv(config_dictionary['albacore_summary_source'], sep="\t")
        self.result_directory = config_dictionary['result_directory']
        self.channel = self.albacore_log['channel']
        self.sequence_length_template = self.albacore_log['sequence_length_template']
        self.null_event = self.albacore_log[self.albacore_log['num_events'] == 0]
        self.albacore_log = self.albacore_log.replace([np.inf, -np.inf], 0)
        self.albacore_log = self.albacore_log[self.albacore_log['num_events'] != 0]
        self.fast5_tot_number = len(self.albacore_log)
        self.is_barcode = config_dictionary['barcoding']

        if self.is_barcode == 'True':
            self.is_barcode = True
        elif self.is_barcode == 'False':
            self.is_barcode = False

        self.my_dpi = int(config_dictionary['dpi'])

        if self.is_barcode:

            self.barcode_selection = config_dictionary['barcode_selection']

            try:
                self.albacore_log.loc[~self.albacore_log['barcode_arrangement'].isin(
                    self.barcode_selection), 'barcode_arrangement'] = 'unclassified'

            except:
                sys.exit('No barcode found in sequencing summary file')

    def get_name(self):
        '''
        Get the name of the extractor.
        :return: the name of the extractor
        '''
        return 'Albacore statistics'

    def init(self):
        '''
        Initialisation
        :return:
        '''
        return

    def check_conf(self):
        '''Configuration checking'''
        return

    def extract(self, result_dict):
        if self.is_barcode:
            self.barcode_selection.append('unclassified')
            for index_barcode, barcode in enumerate(self.barcode_selection):
                barcode_selected_dataframe = self.albacore_log[self.albacore_log['barcode_arrangement'] == barcode]
                result_dict['mean_qscore_statistics_' + barcode] = \
                    barcode_selected_dataframe['mean_qscore_template'].describe()
                result_dict['sequence_length_statistics_' + barcode] = \
                    barcode_selected_dataframe['sequence_length_template'].describe()
        else:

            mean_qscore_template = self.albacore_log['mean_qscore_template']
            result_dict['mean_qscore_statistics'] = pd.DataFrame.describe(mean_qscore_template).drop("count")
            result_dict['sequence_length_statistics'] = self.albacore_log['sequence_length_template'].describe()

        result_dict['channel_occupancy_statistics'] = self._occupancy_channel()
        result_dict['sequence_length_template'] = self.sequence_length_template





    def graph_generation(self):
        '''
        Generation of the differents graphs containing in the graph_generator module
        :return: images array containing the title and the path toward the images
        '''
        images_directory = self.result_directory + '/images'
        images = []
        images.append(graph_generator.read_count_histogram(self.albacore_log, self.my_dpi, images_directory))
        images.append(graph_generator.read_length_histogram(self.albacore_log, self.my_dpi, images_directory))
        images.append(graph_generator.read_number_run(self.albacore_log, self.my_dpi, images_directory))
        images.append(graph_generator.read_quality_boxplot(self.albacore_log, self.my_dpi, images_directory))
        images.append(graph_generator.phred_score_frequency(self.albacore_log, self.my_dpi, images_directory))
        images.append(graph_generator.channel_count_histogram(self.albacore_log, self.my_dpi, images_directory))
        channel_count = self.channel
        total_number_reads_per_pore = pd.value_counts(channel_count)
        images.append(graph_generator.plot_performance(total_number_reads_per_pore, self.my_dpi,
                                                       images_directory))
        images.append(graph_generator.scatterplot(self.albacore_log, self.my_dpi, images_directory))
        if self.is_barcode:
            images.append(graph_generator.barcode_percentage_pie_chart(self.albacore_log, self.barcode_selection,
                                                                             self.my_dpi, images_directory))
            images.append(graph_generator.barcode_length_boxplot(self.albacore_log, self.barcode_selection,
                                                                       self.my_dpi, images_directory))
            images.append(graph_generator.barcoded_phred_score_frequency(self.albacore_log,
                                                                               self.barcode_selection, self.my_dpi,
                                                                               images_directory))
        return images

    def clean(self):
        '''
        Cleaning
        :return:
        '''
        return

    def _occupancy_channel(self):
        '''
        Statistics about the channels
        :return: channel_count_statistics containing statistics description about the channel occupancy
        '''
        channel_count = self.channel
        total_number_reads_per_channel = pd.value_counts(channel_count)
        channel_count_statistics = pd.DataFrame.describe(total_number_reads_per_channel)
        return channel_count_statistics


