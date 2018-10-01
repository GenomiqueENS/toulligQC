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
        self.get_report_data_file_id()

        # panda's object for 1d_summary
        self.albacore_log_1d = pd.read_csv(config_dictionary['albacore_summary_source'], sep="\t")
        self.channel = self.albacore_log_1d['channel']
        self.passes_filtering_1d = self.albacore_log_1d['passes_filtering']
        self.sequence_length_template = self.albacore_log_1d['sequence_length_template']
        self.null_event_1d = self.albacore_log_1d[self.albacore_log_1d['num_events'] == 0]
        self.albacore_log_1d = self.albacore_log_1d.replace([np.inf, -np.inf], 0)
        self.albacore_log_1d = self.albacore_log_1d[self.albacore_log_1d['num_events'] != 0]
        self.fast5_tot_number_1d = len(self.albacore_log_1d)

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

    def get_report_data_file_id(self):
        '''
        Get the report.data id of the extractor.
        :return: the report.data id
        '''
        return 'albacore.stats.1d.extractor'

    def init(self):
        '''
        Initialisation
        :return:
        '''
        return

    def check_conf(self):
        '''Configuration checking'''
        return

    def add_key_to_result_dict(self, key):
        return self.get_report_data_file_id() + '.' + key

    def add_value_to_unwritten_key(self,result_dict, value):
        return result_dict['unwritten.keys'].append(value)

    def extract(self, result_dict):


        # read count
        result_dict[self.add_key_to_result_dict("fastq.entries")] = len(self.albacore_log_1d['num_events'])
        result_dict[self.add_key_to_result_dict("fast5.template.basecalled")] = len(self.albacore_log_1d[self.albacore_log_1d["num_called_template"] != 0])
        result_dict[self.add_key_to_result_dict("read.pass.count")] = len(self.albacore_log_1d[self.albacore_log_1d['passes_filtering'] == True])
        result_dict[self.add_key_to_result_dict("read.fail.count")] = len(self.albacore_log_1d[self.albacore_log_1d['passes_filtering'] == False])

        # read length information
        result_dict[self.add_key_to_result_dict("sequence.length.template")] = self.albacore_log_1d.sequence_length_template[self.albacore_log_1d['num_called_template'] != 0]
        result_dict[self.add_key_to_result_dict("passes.filtering")] = self.albacore_log_1d['passes_filtering']
        result_dict[self.add_key_to_result_dict("read.pass")] = self.albacore_log_1d.sequence_length_template.loc[True == self.albacore_log_1d['passes_filtering']]
        result_dict[self.add_key_to_result_dict("read.fail")] = self.albacore_log_1d.sequence_length_template.loc[False == self.albacore_log_1d['passes_filtering']]

        #yield information
        result_dict[self.add_key_to_result_dict("yield")] = sum(self.albacore_log_1d['sequence_length_template'])
        result_dict[self.add_key_to_result_dict("start.time.sorted")] = sorted(sorted(self.albacore_log_1d['start_time'] / 3600))
        result_dict[self.add_key_to_result_dict("read.pass.sorted")] = sorted(self.albacore_log_1d.start_time.loc[True == self.albacore_log_1d['passes_filtering']]/3600)
        result_dict[self.add_key_to_result_dict("read.fail.sorted")] = sorted(self.albacore_log_1d.start_time.loc[False == self.albacore_log_1d['passes_filtering']]/3600)
        result_dict[self.add_key_to_result_dict("run.time")] = (max(result_dict[self.add_key_to_result_dict("start.time.sorted")]))

        #qscore information
        result_dict[self.add_key_to_result_dict("mean.qscore")] = self.albacore_log_1d.loc[:,"mean_qscore_template"]
        result_dict[self.add_key_to_result_dict("qscore.read.pass")] = self.albacore_log_1d.mean_qscore_template.loc[True == self.albacore_log_1d['passes_filtering']]
        result_dict[self.add_key_to_result_dict("qscore.read.fail")] = self.albacore_log_1d.mean_qscore_template.loc[False == self.albacore_log_1d['passes_filtering']]

        result_dict[self.add_key_to_result_dict('channel.occupancy.statistics')] = self._occupancy_channel()
        channel_occupancy_statistics = result_dict[self.add_key_to_result_dict('channel.occupancy.statistics')]
        for index, value in channel_occupancy_statistics.iteritems():
            result_dict[self.add_key_to_result_dict('channel.occupancy.statistics') + index] = value

        if self.is_barcode:
            self.barcode_selection.append('unclassified')
            pattern = '(\d{2})'
            length = {}
            length['passes_filtering'] = result_dict[self.add_key_to_result_dict("passes.filtering")]
            phred = {}
            phred['passes_filtering'] = result_dict[self.add_key_to_result_dict("passes.filtering")]
            for index_barcode, barcode in enumerate(self.barcode_selection):
                barcode_selected_dataframe = self.albacore_log_1d[self.albacore_log_1d['barcode_arrangement'] == barcode]

                for index,value in barcode_selected_dataframe['mean_qscore_template'].describe().iteritems():
                    result_dict[self.add_key_to_result_dict('') + barcode + '.qscore.' + index] = value

                for index,value in barcode_selected_dataframe['sequence_length_template'].describe().iteritems():
                    result_dict[self.add_key_to_result_dict('') + barcode + '.length.' + index] = value

                match = re.search(pattern, barcode)
                if match:
                    length[match.group(0)] = barcode_selected_dataframe['sequence_length_template']
                    phred[match.group(0)] = barcode_selected_dataframe['mean_qscore_template']
                else:
                    length['Unclassified'] = barcode_selected_dataframe['sequence_length_template']
                    phred['Unclassified'] = barcode_selected_dataframe['mean_qscore_template']


            result_dict[self.add_key_to_result_dict('barcode_selection_sequence_length_dataframe')] = pd.DataFrame(
                dict([(k, pd.Series(v)) for k, v in length.items()]))
            result_dict[self.add_key_to_result_dict('barcode_selection_sequence_length_melted_dataframe')] = pd.melt(
                result_dict[self.add_key_to_result_dict('barcode_selection_sequence_length_dataframe')],
                id_vars=['passes_filtering'],
                var_name="barcodes", value_name="length")

            result_dict[self.add_key_to_result_dict('barcode_selection_sequence_phred_dataframe')] = pd.DataFrame(
                dict([(k, pd.Series(v)) for k, v in phred.items()]))
            result_dict[self.add_key_to_result_dict('barcode_selection_sequence_phred_melted_dataframe')] = pd.melt(
                result_dict[self.add_key_to_result_dict('barcode_selection_sequence_phred_dataframe')],
                id_vars=['passes_filtering'],
                var_name="barcodes", value_name="qscore")

            result_dict[self.add_key_to_result_dict("barcode.arrangement")] = self.albacore_log_1d["barcode_arrangement"]
            result_dict[self.add_key_to_result_dict("read.pass.barcode")] = self.albacore_log_1d.barcode_arrangement.loc[True == self.albacore_log_1d['passes_filtering']]
            result_dict[self.add_key_to_result_dict("read.fail.barcode")] = self.albacore_log_1d.barcode_arrangement.loc[False == self.albacore_log_1d['passes_filtering']]
            length.clear()
            phred.clear()

        else:

            result_dict[self.add_key_to_result_dict('mean.qscore.statistics')] = pd.DataFrame.describe(self.albacore_log_1d['mean_qscore_template']).drop("count")
            for index,value in result_dict[self.add_key_to_result_dict('mean.qscore.statistics')].iteritems():
                result_dict[self.add_key_to_result_dict('mean.qscore.statistics') + index] = value
            result_dict[self.add_key_to_result_dict('sequence.length.statistics')] = self.albacore_log_1d['sequence_length_template'].describe()
            for index,value in result_dict[self.add_key_to_result_dict('sequence.length.statistics')].iteritems():
                result_dict[self.add_key_to_result_dict('sequence.length.statistics') + index] = value



    def graph_generation(self,result_dict):
        '''
        Generation of the differents graphs containing in the graph_generator module
        :return: images array containing the title and the path toward the images
        '''
        images_directory = self.result_directory + '/images'
        images = []
        images.append(graph_generator.read_count_histogram(result_dict, 'Read count histogram', self.my_dpi,images_directory, "Number of reads produced before (Fast 5 in blue) and after (1D in orange) basecalling. The basecalled reads are filtered with a 7.5 quality score threshold in pass (1D pass in green) or fail (1D fail in red) categories."))
        images.append(graph_generator.read_length_multihistogram(result_dict, 'Read length histogram', self.my_dpi,images_directory,"Size distribution of basecalled reads (1D in orange). The basecalled reads are filtered with a 7.5 quality score threshold in pass (1D pass in green) or fail (1D fail in red) categories."))

        images.append(graph_generator.allread_number_run(result_dict, 'Yield plot of 1D read type', self.my_dpi,images_directory,"Yield plot of basecalled reads (1D in orange). The basecalled reads are filtered with a 7.5 quality score threshold in pass (1D pass in green) or fail (1D fail in red) categories."))

        images.append(graph_generator.read_quality_multiboxplot(result_dict, "Read type quality boxplot", self.my_dpi,images_directory,"Boxplot of 1D reads (in orange) quality.  The basecalled reads are filtered with a 7.5 quality score threshold in pass (1D pass in green) or fail (1D fail in red) categories."))

        images.append(graph_generator.phred_score_frequency(result_dict, 'Mean Phred score frequency of 1D reads', self.my_dpi, images_directory,"Phred score frequency of 1D reads."))
        images.append(graph_generator.allphred_score_frequency(result_dict, 'Mean Phred score frequency of all 1D read type', self.my_dpi,images_directory,"The basecalled reads are filtered with a 7.5 quality score threshold in pass (1D pass in green) or fail (1D fail in red) categories."))

        #images.append(graph_generator.channel_count_histogram(self.albacore_log_1d, 'Channel occupancy', self.my_dpi, images_directory))
        channel_count = self.channel
        total_number_reads_per_pore = pd.value_counts(channel_count)
        images.append(graph_generator.plot_performance(total_number_reads_per_pore, 'Channel occupancy of the flowcell', self.my_dpi,images_directory,"Number of reads sequenced per pore channel."))

        images.append(graph_generator.all_scatterplot(result_dict, 'Mean Phred score function of 1D read length',self.my_dpi, images_directory,"The Mean Phred score varies according to the read length. The basecalled reads are filtered with a 7.5 quality score threshold in pass (1D pass in green) or fail (1D fail in red) categories."))

        if self.is_barcode:
            images.append(graph_generator.barcode_percentage_pie_chart_pass(result_dict,'1D pass reads percentage of different barcodes', self.barcode_selection,
                                                                             self.my_dpi, images_directory,"1D pass read distribution per barcode."))
            images.append(graph_generator.barcode_percentage_pie_chart_fail(result_dict,'1D fail reads percentage of different barcodes', self.barcode_selection,
                                                                             self.my_dpi, images_directory,"1D fail read distribution per barcode."))
            images.append(graph_generator.barcode_length_boxplot(result_dict,'1D reads size distribution for each barcode',
                                                                       self.my_dpi, images_directory,"Read length boxplot per barcode of pass (in green) and fail (in red) 1D reads."))
            images.append(graph_generator.barcoded_phred_score_frequency(result_dict,'1D reads Mean Phred score distribution for each barcode',
                                                                               self.my_dpi,
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


