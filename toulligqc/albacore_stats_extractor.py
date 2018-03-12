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
        self.result_directory = config_dictionary['result_directory']
        self.is_barcode = config_dictionary['barcoding']

        # panda's object for 1d_summary
        self.albacore_log_1d = pd.read_csv(config_dictionary['albacore_summary_source'], sep="\t")
        self.channel = self.albacore_log_1d['channel']
        self.passes_filtering_1d = self.albacore_log_1d['passes_filtering']
        self.sequence_length_template = self.albacore_log_1d['sequence_length_template']
        self.null_event_1d = self.albacore_log_1d[self.albacore_log_1d['num_events'] == 0]
        self.albacore_log_1d = self.albacore_log_1d.replace([np.inf, -np.inf], 0)
        self.albacore_log_1d = self.albacore_log_1d[self.albacore_log_1d['num_events'] != 0]
        self.fast5_tot_number_1d = len(self.albacore_log_1d)
        self.albacore_log_1d["Yield"]=sum(self.albacore_log_1d['sequence_length_template'])

        if self.is_barcode == 'True':
            self.is_barcode = True
        elif self.is_barcode == 'False':
            self.is_barcode = False

        self.my_dpi = int(config_dictionary['dpi'])

        if self.is_barcode:

            self.barcode_selection = config_dictionary['barcode_selection']
            try:
                self.albacore_log_1d.loc[~self.albacore_log_1d['barcode_arrangement'].isin(
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
                barcode_selected_dataframe = self.albacore_log_1d[self.albacore_log_1d['barcode_arrangement'] == barcode]
                result_dict['mean_qscore_statistics_' + barcode] = \
                    barcode_selected_dataframe['mean_qscore_template'].describe()
                result_dict['sequence_length_statistics_' + barcode] = \
                    barcode_selected_dataframe['sequence_length_template'].describe()
        else:

            mean_qscore_template = self.albacore_log_1d['mean_qscore_template']
            result_dict['mean_qscore_statistics'] = pd.DataFrame.describe(mean_qscore_template).drop("count")
            result_dict['sequence_length_statistics'] = self.albacore_log_1d['sequence_length_template'].describe()

        result_dict['channel_occupancy_statistics'] = self._occupancy_channel()
        result_dict['sequence_length_template'] = self.sequence_length_template
        result_dict['yield']=sum(self.sequence_length_template)

    def graph_generation(self):
        '''
        Generation of the differents graphs containing in the graph_generator module
        :return: images array containing the title and the path toward the images
        '''

        #  Yield plot of basecalled reads (1D in orange). The basecalled reads are filtered with a 7.5 quality score threshold in pass (1D pass in green) or fail (1D fail in red) categories.
        # Boxplot of 1D reads (in orange) quality.  The basecalled reads are filtered with a 7.5 quality score threshold in pass (1D pass in green) or fail (1D fail in red) categories.
        #  Phred score frequency of 1D reads.
        # The basecalled reads are filtered with a 7.5 quality score threshold in pass (1D pass in green) or fail (1D fail in red) categories.
        #  Number of reads sequenced per pore channel.
        #  The Mean Phred score varies according to the read length. The basecalled reads are filtered with a 7.5 quality score threshold in pass (1D pass in green) or fail (1D fail in red) categories.
        #  1D read pass distribution per barcode.
        #  1D read fail distribution per barcode.
        #  Read length boxplot per barcode of pass (in green) and fail (in red) 1D reads.

        images_directory = self.result_directory + '/images'
        images = []
        images.append(graph_generator.read_count_histogram(self.albacore_log_1d, 'Read count histogram', self.my_dpi,images_directory, "Number of reads produced before (Fast 5 in blue) and after (1D in orange) basecalling. The basecalled reads are filtered with a 7.5 quality score threshold in pass (1D pass in green) or fail (1D fail in red) categories."))

        #images.append(graph_generator.read_length_histogram(self.albacore_log_1d, 'Read size histogram', self.my_dpi,images_directory))
        images.append(graph_generator.read_length_multihistogram(self.albacore_log_1d, 'Read length histogram', self.my_dpi,images_directory,"Size distribution of basecalled reads (1D in orange). The basecalled reads are filtered with a 7.5 quality score threshold in pass (1D pass in green) or fail (1D fail in red) categories."))

        images.append(graph_generator.allread_number_run(self.albacore_log_1d, 'Yield plot of 1D read type', self.my_dpi,images_directory,"Yield plot of basecalled reads (1D in orange). The basecalled reads are filtered with a 7.5 quality score threshold in pass (1D pass in green) or fail (1D fail in red) categories."))

        images.append(graph_generator.read_quality_multiboxplot(self.albacore_log_1d, "Read type quality boxplot", self.my_dpi,images_directory,"Boxplot of 1D reads (in orange) quality.  The basecalled reads are filtered with a 7.5 quality score threshold in pass (1D pass in green) or fail (1D fail in red) categories."))

        images.append(graph_generator.phred_score_frequency(self.albacore_log_1d, 'Mean Phred score frequency of 1D reads', self.my_dpi, images_directory,"Phred score frequency of 1D reads."))
        images.append(graph_generator.allphred_score_frequency(self.albacore_log_1d, 'Mean Phred score frequency of all 1D read type', self.my_dpi,images_directory,"The basecalled reads are filtered with a 7.5 quality score threshold in pass (1D pass in green) or fail (1D fail in red) categories."))

        #images.append(graph_generator.channel_count_histogram(self.albacore_log_1d, 'Channel occupancy', self.my_dpi, images_directory))
        channel_count = self.channel
        total_number_reads_per_pore = pd.value_counts(channel_count)
        images.append(graph_generator.plot_performance(total_number_reads_per_pore, 'Channel occupancy of the flowcell', self.my_dpi,images_directory,"Number of reads sequenced per pore channel."))

        images.append(graph_generator.all_scatterplot(self.albacore_log_1d, 'Mean Phred score function of 1D read length',self.my_dpi, images_directory,"The Mean Phred score varies according to the read length. The basecalled reads are filtered with a 7.5 quality score threshold in pass (1D pass in green) or fail (1D fail in red) categories."))

        if self.is_barcode:
            images.append(graph_generator.barcode_percentage_pie_chart_pass(self.albacore_log_1d,'1D pass reads percentage of different barcodes', self.barcode_selection,
                                                                             self.my_dpi, images_directory,"1D pass read distribution per barcode."))
            images.append(graph_generator.barcode_percentage_pie_chart_fail(self.albacore_log_1d,'1D fail reads percentage of different barcodes', self.barcode_selection,
                                                                             self.my_dpi, images_directory,"1D fail read distribution per barcode."))
            images.append(graph_generator.barcode_length_boxplot(self.albacore_log_1d,'1D reads size distribution for each barcode',  self.barcode_selection,
                                                                       self.my_dpi, images_directory,"Read length boxplot per barcode of pass (in green) and fail (in red) 1D reads."))
            images.append(graph_generator.barcoded_phred_score_frequency(self.albacore_log_1d,'1D reads Mean Phred score distribution for each barcode',
                                                                               self.barcode_selection, self.my_dpi,
                                                                        images_directory,"Read Mean Phred score boxplot per barcode of pass (in green) and fail (in red) 1D reads."))
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


