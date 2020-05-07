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

# Extraction of statistics from sequencing_summary.txt file (1D chemistry)

import pandas as pd
import sys
import graph_generator
import numpy as np
import re
import os.path


class SequencingSummaryExtractor:

    """
    Extraction of data from sequencing_summary.txt and optional barcoding files.
    The data is extracted from dataframes and placed in the result_dict in the form of key-value pairs
    """

    def __init__(self, config_dictionary):
        """
        Constructor that initialize the values of the config_dictionary and check in the case of 1 argument in 
        sequencing_summary_source if the path points to a file, the others cases are managed in check_conf 
        and _load_sequencing_summary_data methods
        :param config_dictionary: dictionary containing all files or directories paths for sequencing_summary.txt and barcoding files
        """

        self.config_dictionary = config_dictionary
        self.sequencing_summary_source = self.config_dictionary['sequencing_summary_source']
        self.result_directory = config_dictionary['result_directory']
        self.sequencing_summary_files = self.sequencing_summary_source.split('\t')

        if len(self.sequencing_summary_files) == 1:
            if os.path.isfile(self.sequencing_summary_source):
                self.sequencing_summary_files = [self.sequencing_summary_source]
            elif os.path.isdir(self.sequencing_summary_source):
                raise ValueError("The sequencing summary file must be a file path not a directory path")

        if config_dictionary['barcoding'] == 'True':
            self.is_barcode = True
        else:
            self.is_barcode = False

        self.my_dpi = int(self.config_dictionary['dpi'])


    def check_conf2(self):
        """
        Check if the sequencing summary source contains files
        If true, check for sequencing_summary_file within the sequencing_summary_source
        """

        if len(self.sequencing_summary_files) == 0:
            return False, "No sequencing summary file has been defined"

        found = False
        while not found:
            for f in self.sequencing_summary_files:
                try:
                    if self._is_sequencing_summary_file(f):
                        found = True
                except IsADirectoryError as exception:
                        return False, "There is a path that points to a directory : " + f
        if not found:
            return False, "There is no sequencing summary file in sequencing_summary_source"
        else:
            return found, ""
                        

    def check_conf(self):
        """
        Check if the sequencing summary source contains files
        If true, check for sequencing_summary_file within the sequencing_summary_source
        """

        if len(self.sequencing_summary_files) == 0:
            return False, "No sequencing summary file has been defined"

        for f in self.sequencing_summary_files:
            if not os.path.isfile(f):
                return False, "The sequencing summary file is not a file: " + f

            if self._is_sequencing_summary_file(f):
                return True, ""

        return False, "There is no sequencing summary file in sequencing_summary_source"


    def init(self):
        """
        Creation of the dataframe containing all info from sequencing_summary.txt
        :return: pd.Dataframe object
        """
        #TODO: garder le nom "dataframe_1d" ??
        self.dataframe_1d = self._load_sequencing_summary_data()
        self.channel = self.dataframe_1d['channel']
        self.passes_filtering_1d = self.dataframe_1d['passes_filtering']
        self.sequence_length_template = self.dataframe_1d['sequence_length_template']
        self.null_event_1d = self.dataframe_1d[self.dataframe_1d['num_events'] == 0]
        self.dataframe_1d = self.dataframe_1d[self.dataframe_1d['num_events'] != 0]

        if self.is_barcode:
            self.barcode_selection = self.config_dictionary['barcode_selection']

    @staticmethod
    def get_name():
        """
        Get the name of the extractor.
        :return: the name of the extractor
        """
        return 'Basecaller sequencing summary'


    @staticmethod
    def get_report_data_file_id():
        """
        Get the report.data id of the extractor.
        :return: a string with the name of the ID of the extractor
        """
        #changer la valeur du report data file ??
        return 'basecaller.sequencing.summary.1d.extractor'


    def add_key_to_result_dict(self, key):
        """
        Add a key to the result_dict dictionary with the module name as a prefix
        :param key: key suffix
        :return: result_dict entry (string)
        """
        return self.get_report_data_file_id() + '.' + key

    @staticmethod
    def describe_dict(result_dict, prefix):
        """
        Gives basic statistic like mean, min, max, median and percentiles (without the count value) filled in the result_dict dictionary
        :param result_dict: result_dict dictionary
        :param prefix: Key prefix
        :return: result_dict dictionary filled with statistics values
        """
        dictionary_statistics = pd.Series.describe(result_dict[prefix]).drop("count")
        for key in dict(dictionary_statistics):
            result_dict[prefix + '.' + key] = dictionary_statistics[key]


    def barcode_frequency(self, result_dict, entry, prefix=''):
        """
        Count reads by values of barcode_selection, computes sum of counts by barcode_selection, and sum of unclassified counts.
        Regroup all non used barcodes in index "other"
        Compute all frequency values for each number of barcoded reads
        :param entry: entry about barcoded counts
        :param prefix: key prefix
        :return: Series with all barcodes (used, non used, and unclassified) frequencies
        """
        # Regroup all barcoded read in Series
        all_barcode_count = result_dict[self.add_key_to_result_dict(entry)].value_counts()
        # Sort by list of used barcodes
        count_sorted = all_barcode_count.sort_index()[self.barcode_selection]
        # Compute sum of all used barcodes without barcode 'unclassified'
        result_dict[self.add_key_to_result_dict(prefix + 'barcoded.count')] = sum(count_sorted.drop("unclassified"))
        # Compute all reads of barcodes that are not in the barcode_selection list
        result_dict[self.add_key_to_result_dict(prefix + "non.used.barcodes.count")] = (sum(all_barcode_count)-sum(count_sorted))
        # Create Series for all non-used barcode counts and rename index array with "other"
        other_all_barcode_count = pd.Series(result_dict[self.add_key_to_result_dict(prefix + "non.used.barcodes.count")], index=['other'])
        # Append Series of non-used barcodes to the Series with the counts of barcode_selection
        count_sorted = count_sorted.append(other_all_barcode_count).sort_index()

        # Compute frequency for all counts by barcodes
        for key in count_sorted.to_dict():
            result_dict[self.add_key_to_result_dict(prefix) + key + ".frequency"] = \
            count_sorted[key]*100/sum(count_sorted)
                
        return count_sorted


    def extract_old(self, result_dict):
        """
        Get Phred score (Qscore) and Length details per read
        :param result_dict:
        :return:
        """
        # Basecaller analysis
        if 'sequencing.telemetry.extractor.software.analysis' not in result_dict:
            result_dict['sequencing.telemetry.extractor.software.analysis'] = '1d_basecalling'

        # Read count
        result_dict[self.add_key_to_result_dict("fastq.entries")] = len(self.dataframe_1d['num_events'])

        result_dict[self.add_key_to_result_dict("read.count")] = \
            len(self.dataframe_1d[self.dataframe_1d["num_events_template"] != 0])

        result_dict[self.add_key_to_result_dict("read.with.length.equal.zero.count")] = \
            len(self.dataframe_1d[self.dataframe_1d['sequence_length_template'] == 0])

        # 1D pass information
        result_dict[self.add_key_to_result_dict("read.pass.count")] = \
            len(self.dataframe_1d.loc[self.dataframe_1d['passes_filtering'] == bool(True)])
        result_dict[self.add_key_to_result_dict("read.pass.length")] = \
            self.dataframe_1d.sequence_length_template.loc[self.dataframe_1d['passes_filtering'] == bool(True)]
        result_dict[self.add_key_to_result_dict("read.pass.sorted")] = \
            sorted(self.dataframe_1d.start_time.loc[self.dataframe_1d['passes_filtering'] == bool(True)] / 3600)
        result_dict[self.add_key_to_result_dict("read.pass.qscore")] = \
            self.dataframe_1d.mean_qscore_template.loc[self.dataframe_1d['passes_filtering'] == bool(True)]

        # 1D fail information
        result_dict[self.add_key_to_result_dict("read.fail.count")] = \
            len(self.dataframe_1d.loc[self.dataframe_1d['passes_filtering'] == bool(False)])
        result_dict[self.add_key_to_result_dict("read.fail.length")] = \
            self.dataframe_1d.sequence_length_template.loc[self.dataframe_1d['passes_filtering'] == bool(False)]
        result_dict[self.add_key_to_result_dict("read.fail.sorted")] = \
            sorted(self.dataframe_1d.start_time.loc[self.dataframe_1d['passes_filtering'] == bool(False)] / 3600)
        result_dict[self.add_key_to_result_dict("read.fail.qscore")] = \
            self.dataframe_1d.mean_qscore_template.loc[self.dataframe_1d['passes_filtering'] == bool(False)]

        # Read proportion
        result_dict[self.add_key_to_result_dict("fastq.entries.ratio")] = \
            result_dict[self.add_key_to_result_dict('fastq.entries')] / \
            result_dict[self.add_key_to_result_dict('fastq.entries')]

        result_dict[self.add_key_to_result_dict("read.count.ratio")] = \
            result_dict[self.add_key_to_result_dict("read.count")] / \
            result_dict[self.add_key_to_result_dict("read.count")]

        result_dict[self.add_key_to_result_dict("read.with.length.equal.zero.ratio")] = \
            result_dict[self.add_key_to_result_dict("read.with.length.equal.zero.count")] / \
            result_dict[self.add_key_to_result_dict("read.count")]

        result_dict[self.add_key_to_result_dict("read.pass.ratio")] = \
            result_dict[self.add_key_to_result_dict("read.pass.count")] / \
            result_dict[self.add_key_to_result_dict("read.count")]

        result_dict[self.add_key_to_result_dict("read.fail.ratio")] = \
            result_dict[self.add_key_to_result_dict("read.fail.count")] / \
            result_dict[self.add_key_to_result_dict("read.count")]

        result_dict[self.add_key_to_result_dict("fastq.entries.frequency")] = \
            result_dict[self.add_key_to_result_dict('fastq.entries')] / \
            result_dict[self.add_key_to_result_dict('fastq.entries')]*100

        result_dict[self.add_key_to_result_dict("read.count.frequency")] = \
            result_dict[self.add_key_to_result_dict("read.count")] / \
            result_dict[self.add_key_to_result_dict("read.count")]*100

        result_dict[self.add_key_to_result_dict("read.with.length.equal.zero.frequency")] = \
            result_dict[self.add_key_to_result_dict("read.with.length.equal.zero.count")] / \
            result_dict[self.add_key_to_result_dict("read.count")]*100

        result_dict[self.add_key_to_result_dict("read.pass.frequency")] = \
            result_dict[self.add_key_to_result_dict("read.pass.count")] / \
            result_dict[self.add_key_to_result_dict("read.count")]*100

        result_dict[self.add_key_to_result_dict("read.fail.frequency")] = \
            result_dict[self.add_key_to_result_dict("read.fail.count")] / \
            result_dict[self.add_key_to_result_dict("read.count")]*100

        # Read length information
        result_dict[self.add_key_to_result_dict("sequence.length")] = \
            self.dataframe_1d.sequence_length_template[self.dataframe_1d["num_events_template"] != 0]

        result_dict[self.add_key_to_result_dict("passes.filtering")] = \
            self.dataframe_1d['passes_filtering']

        # Yield
        result_dict[self.add_key_to_result_dict("yield")] = sum(self.dataframe_1d['sequence_length_template'])

        result_dict[self.add_key_to_result_dict("start.time.sorted")] = \
            sorted(sorted(self.dataframe_1d['start_time'] / 3600))

        result_dict[self.add_key_to_result_dict("run.time")] = \
            (max(result_dict[self.add_key_to_result_dict("start.time.sorted")]))

        # Qscore information
        result_dict[self.add_key_to_result_dict("mean.qscore")] = self.dataframe_1d.loc[:, "mean_qscore_template"]

        # Channel occupancy information
        result_dict[self.add_key_to_result_dict('channel.occupancy.statistics')] = self._occupancy_channel()
        channel_occupancy_statistics = result_dict[self.add_key_to_result_dict('channel.occupancy.statistics')]

        for index, value in channel_occupancy_statistics.iteritems():
            result_dict[self.add_key_to_result_dict('channel.occupancy.statistics.') + index] = value

        # Length's statistic information provided in the result_dict
        result_dict[self.add_key_to_result_dict('all.read.length')] = \
            self.dataframe_1d['sequence_length_template'].describe()

        for index, value in result_dict[self.add_key_to_result_dict('all.read.length')].iteritems():
            result_dict[self.add_key_to_result_dict('all.read.length.') + index] = value
        self.describe_dict(result_dict, self.add_key_to_result_dict("read.pass.length"))
        self.describe_dict(result_dict, self.add_key_to_result_dict("read.fail.length"))

        # Qscore's statistic information provided in the result_dict
        result_dict[self.add_key_to_result_dict('all.read.qscore')] = \
            pd.DataFrame.describe(self.dataframe_1d['mean_qscore_template']).drop("count")

        for index, value in result_dict[self.add_key_to_result_dict('all.read.qscore')].iteritems():
            result_dict[self.add_key_to_result_dict('all.read.qscore.') + index] = value
        self.describe_dict(result_dict, self.add_key_to_result_dict("read.pass.qscore"))
        self.describe_dict(result_dict, self.add_key_to_result_dict("read.fail.qscore"))

        # In case of barcoded samples
        if self.is_barcode:

            self.barcode_selection.append('unclassified')

            result_dict[self.add_key_to_result_dict("barcode.arrangement")] = \
                self.dataframe_1d["barcode_arrangement"]
            result_dict[self.add_key_to_result_dict("read.pass.barcode")] = \
                self.dataframe_1d.barcode_arrangement.loc[self.dataframe_1d['passes_filtering'] == bool(True)]
            result_dict[self.add_key_to_result_dict("read.fail.barcode")] = \
                self.dataframe_1d.barcode_arrangement.loc[self.dataframe_1d['passes_filtering'] == bool(False)]

            # Get barcodes frequency by read type
            result_dict[self.add_key_to_result_dict("all.read.barcoded")] = self.barcode_frequency(result_dict, "barcode.arrangement", 'all.read.')
            result_dict[self.add_key_to_result_dict("read.pass.barcoded")] = self.barcode_frequency(result_dict, "read.pass.barcode", 'read.pass.')
            result_dict[self.add_key_to_result_dict("read.fail.barcoded")] = self.barcode_frequency(result_dict, "read.fail.barcode", 'read.fail.')

            result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.frequency"] = result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count"]/result_dict["basecaller.sequencing.summary.1d.extractor.read.count"] *100
            result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.frequency"] = result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.count"]/result_dict["basecaller.sequencing.summary.1d.extractor.read.count"] *100

            pattern = '(\d{2})'

            self.dataframe_1d.loc[~self.dataframe_1d['barcode_arrangement'].isin(
                self.barcode_selection), 'barcode_arrangement'] = 'other'

            result_dict[self.add_key_to_result_dict("passes.filtering")] = \
                self.dataframe_1d['passes_filtering']


            length = {'passes_filtering': result_dict[self.add_key_to_result_dict("passes.filtering")]}
            phred = {'passes_filtering': result_dict[self.add_key_to_result_dict("passes.filtering")]}


            self.barcode_selection.append('other')

            for index_barcode, barcode in enumerate(self.barcode_selection):

                barcode_selected_dataframe = \
                    self.dataframe_1d[self.dataframe_1d['barcode_arrangement'] == barcode]

                barcode_selected_read_pass_dataframe = \
                    barcode_selected_dataframe.loc[self.dataframe_1d['passes_filtering'] == bool(True)]

                barcode_selected_read_fail_dataframe = \
                    barcode_selected_dataframe.loc[self.dataframe_1d['passes_filtering'] == bool(False)]

                match = re.search(pattern, barcode)

                if match:
                    length[match.group(0)] = barcode_selected_dataframe['sequence_length_template']
                    phred[match.group(0)] = barcode_selected_dataframe['mean_qscore_template']

                    for index, value in barcode_selected_dataframe['sequence_length_template']\
                            .describe().iteritems():
                        result_dict[self.add_key_to_result_dict('all.read.') + barcode + '.length.' + index] = value

                    for index, value in barcode_selected_read_pass_dataframe['sequence_length_template']\
                            .describe().iteritems():
                        result_dict[self.add_key_to_result_dict('read.pass.') + barcode + '.length.' + index] = value

                    for index, value in barcode_selected_read_fail_dataframe['sequence_length_template']\
                            .describe().iteritems():
                        result_dict[self.add_key_to_result_dict('read.fail.') + barcode + '.length.' + index] = value

                    for index, value in barcode_selected_dataframe['mean_qscore_template']\
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('all.read.') + barcode + '.qscore.' + index] = value

                    for index, value in barcode_selected_read_pass_dataframe['mean_qscore_template']\
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('read.pass.') + barcode + '.qscore.' + index] = value

                    for index, value in barcode_selected_read_fail_dataframe['mean_qscore_template']\
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('read.fail.') + barcode + '.qscore.' + index] = value
                if barcode == 'unclassified':
                    length['Unclassified'] = barcode_selected_dataframe['sequence_length_template']
                    phred['Unclassified'] = barcode_selected_dataframe['mean_qscore_template']

                    for index, value in barcode_selected_dataframe['sequence_length_template'].describe().iteritems():
                        result_dict[self.add_key_to_result_dict('all.read.unclassified.length.') + index] = value

                    for index, value in barcode_selected_read_pass_dataframe['sequence_length_template']\
                            .describe().iteritems():
                        result_dict[self.add_key_to_result_dict('read.pass.unclassified.length.') + index] = value

                    for index, value in barcode_selected_read_fail_dataframe['sequence_length_template']\
                            .describe().iteritems():
                        result_dict[self.add_key_to_result_dict('read.fail.unclassified.length.') + index] = value

                    for index, value in barcode_selected_dataframe['mean_qscore_template']\
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('all.read.unclassified.qscore.') + index] = value

                    for index, value in barcode_selected_read_pass_dataframe['mean_qscore_template']\
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('read.pass.unclassified.qscore.') + index] = value

                    for index, value in barcode_selected_read_fail_dataframe['mean_qscore_template']\
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('read.fail.unclassified.qscore.') + index] = value

                if barcode == 'other':
                    length['Other Barcodes'] = barcode_selected_dataframe['sequence_length_template']
                    phred['Other Barcodes'] = barcode_selected_dataframe['mean_qscore_template']

                    for index, value in barcode_selected_dataframe['sequence_length_template'].describe().iteritems():
                        result_dict[self.add_key_to_result_dict('all.read.with.other.barcodes.length.') + index] = value

                    for index, value in barcode_selected_read_pass_dataframe['sequence_length_template'] \
                            .describe().iteritems():
                        result_dict[self.add_key_to_result_dict('read.pass.with.other.barcodes.length.') + index] = value

                    for index, value in barcode_selected_read_fail_dataframe['sequence_length_template'] \
                            .describe().iteritems():
                        result_dict[self.add_key_to_result_dict('read.fail.with.other.barcodes.length.') + index] = value

                    for index, value in barcode_selected_dataframe['mean_qscore_template'] \
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('all.read.with.other.barcodes.qscore.') + index] = value

                    for index, value in barcode_selected_read_pass_dataframe['mean_qscore_template'] \
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('read.pass.with.other.barcodes.qscore.') + index] = value

                    for index, value in barcode_selected_read_fail_dataframe['mean_qscore_template'] \
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('read.fail.with.other.barcodes.qscore.') + index] = value

            # Provide statistic per barcode in the result_dict dictionary

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

            length.clear()
            phred.clear()


    def _count_elements(self, dataframe, column):
        """
        Returns the number of values in the dataframe's column
        """
        return len(dataframe[column])
    
        
    def _count_not_zero_elements(self, dataframe, column):
        """
        Returns the number of values different than zero in the dataframe's column
        """
        return len(dataframe[dataframe[column] != 0])
    
    
    def _count_zero_elements(self, dataframe, column):
        """
        Returns the number of values equals to zero in the dataframe's column
        """
        return len(dataframe[dataframe[column] == 0])
    
    
    def _count_boolean_elements(self, dataframe, column, boolean_value):
        """
        Returns the number of values of a column filtered by a boolean
        """
        return len(dataframe.loc[dataframe[column] == bool(boolean_value)])
    
    
    def _series_cols_boolean_elements(self, dataframe, column1, column2, boolean_value):
        """
        Returns the number of values of different columns filtered by a boolean
        """
        return dataframe[column1].loc[dataframe[column2] == bool(boolean_value)]
    
    
    def _count_sorted_boolean_elements(self, dataframe, column1, column2, boolean_value):
        """
        Returns a sorted list of values of different columns filtered by a boolean
        """
        return sorted(dataframe[column1].loc[dataframe[column2] == bool(boolean_value)])
    
    
    def _count_sorted_boolean_elements_divided(self, dataframe, column1, column2, boolean_value, denominator):
        """
        Returns a sorted list of values of different columns filtered by a boolean and divided by the denominator
        """
        return sorted(dataframe[column1].loc[dataframe[column2] == bool(boolean_value)] / denominator)
        

    def _set_result_to_dict(self, result_dict, key, function, *args):
        """
        Add a new item in result_dict : {key : value returned by function}
        """
        result_dict[self.add_key_to_result_dict(key)] = function
    
    
    def extract(self, result_dict):
        """
        Get Phred score (Qscore) and Length details per read
        :param result_dict:
        :return: None, info is saved in result_dict
        """
        # Basecaller analysis

        if 'sequencing.telemetry.extractor.software.analysis' not in result_dict:
            result_dict['sequencing.telemetry.extractor.software.analysis'] = '1d_basecalling'

        # Fastq entries
        self._set_result_to_dict(result_dict, "fastq.entries", self._count_elements(self.dataframe_1d, 'num_events'))
        # Read count
        self._set_result_to_dict(result_dict, "read.count", self._count_not_zero_elements(self.dataframe_1d, "num_events_template"))
        # Read count with length equals zero
        self._set_result_to_dict(result_dict, "read.with.length.equal.zero.count", self._count_zero_elements(self.dataframe_1d, 'sequence_length_template'))
        
        # 1D pass information : count, length, qscore values and sorted Series
        self._set_result_to_dict(result_dict, "read.pass.count", self._count_boolean_elements(self.dataframe_1d, 'passes_filtering', True))
        self._set_result_to_dict(result_dict, "read.pass.length", self._series_cols_boolean_elements(self.dataframe_1d, 'sequence_length_template', 'passes_filtering', True))
        self._set_result_to_dict(result_dict, "read.pass.qscore", self._series_cols_boolean_elements(self.dataframe_1d, 'mean_qscore_template', 'passes_filtering', True))
        self._set_result_to_dict(result_dict, "read.pass.sorted", self._count_sorted_boolean_elements_divided(self.dataframe_1d, 'start_time', 'passes_filtering', True, 3600))
        
        # 1D fail information : count, length, qscore values and sorted Series
        self._set_result_to_dict(result_dict, "read.fail.count", self._count_boolean_elements(self.dataframe_1d, 'passes_filtering', False))
        self._set_result_to_dict(result_dict, "read.fail.length", self._series_cols_boolean_elements(self.dataframe_1d, 'sequence_length_template', 'passes_filtering', False))
        self._set_result_to_dict(result_dict, "read.fail.qscore", self._series_cols_boolean_elements(self.dataframe_1d, 'mean_qscore_template', 'passes_filtering', False))
        self._set_result_to_dict(result_dict, "read.fail.sorted", self._count_sorted_boolean_elements_divided(self.dataframe_1d, 'start_time', 'passes_filtering', False, 3600))
        
        total_reads = result_dict.get("basecaller.sequencing.summary.1d.extractor.read.count")
        
        # Ratios 
        result_dict[self.add_key_to_result_dict("read.with.length.equal.zero.ratio")] = \
        result_dict.get("basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.count") / total_reads

        result_dict[self.add_key_to_result_dict("read.pass.ratio")] = \
        result_dict.get("basecaller.sequencing.summary.1d.extractor.read.pass.count") / total_reads

        result_dict[self.add_key_to_result_dict("read.fail.ratio")] = \
        result_dict.get("basecaller.sequencing.summary.1d.extractor.read.fail.count") / total_reads

        # Frequencies
        # delete those 2 lines after redoing read count histogram (values always equal to 100)
        result_dict[self.add_key_to_result_dict("fastq.entries.frequency")] = 100
        result_dict[self.add_key_to_result_dict("read.count.frequency")] = 100
        
        result_dict[self.add_key_to_result_dict("read.with.length.equal.zero.frequency")] = \
        (result_dict.get("basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.count") / total_reads) * 100

        result_dict[self.add_key_to_result_dict("read.pass.frequency")] = \
        (result_dict.get("basecaller.sequencing.summary.1d.extractor.read.pass.count") / total_reads) * 100
        
        result_dict[self.add_key_to_result_dict("read.fail.frequency")] = \
        (result_dict.get("basecaller.sequencing.summary.1d.extractor.read.fail.count") / total_reads) * 100

        # Read length information
        result_dict[self.add_key_to_result_dict("sequence.length")] = \
            self.dataframe_1d.sequence_length_template[self.dataframe_1d["num_events_template"] != 0]

        result_dict[self.add_key_to_result_dict("passes.filtering")] = self.dataframe_1d['passes_filtering']
        
        # Yield
        result_dict[self.add_key_to_result_dict("yield")] = sum(self.dataframe_1d['sequence_length_template'])

        start_time_sorted = sorted(self.dataframe_1d['start_time'] / 3600)
        self._set_result_to_dict(result_dict, "start.time.sorted", sorted(start_time_sorted))

        result_dict[self.add_key_to_result_dict("run.time")] = max(result_dict.get("basecaller.sequencing.summary.1d.extractor.start.time.sorted"))
        
        # Retrieve Qscore column information and save it in mean.qscore entry
        result_dict[self.add_key_to_result_dict("mean.qscore")] = self.dataframe_1d['mean_qscore_template']
        
        # Get channel occupancy statistics and store each value into result_dict
        channel_statistics = self._occupancy_channel()
        for index, value in channel_statistics.items():
            result_dict[self.add_key_to_result_dict('channel.occupancy.statistics.') + index] = value

        # Get statistics about all reads and store each value into result_dict
        sequence_length_statistics = self.dataframe_1d['sequence_length_template'].describe()
        
        for index, value in sequence_length_statistics.items():
            result_dict[self.add_key_to_result_dict('all.read.length.') + index] = value
        
        # Add statistics (without count) about read pass/fail length in the result_dict
        self.describe_dict(result_dict, "basecaller.sequencing.summary.1d.extractor.read.pass.length")
        self.describe_dict(result_dict, "basecaller.sequencing.summary.1d.extractor.read.fail.length")
        
        # Get Qscore statistics without count value and store them into result_dict
        qscore_statistics = self.dataframe_1d['mean_qscore_template'].describe().drop("count")
        
        for index, value in qscore_statistics.items():
            result_dict[self.add_key_to_result_dict('all.read.qscore.') + index] = value
        
        if self._is_barcode_file:
            self.extract_barcode_info(result_dict)

          
    def extract_barcode_info(self, result_dict):
        """
        Gather barcode information for graphs
        """
        # Add values unclassified and other to barcode list
        self.barcode_selection.append('unclassified')
        
        # Create keys barcode_arrangement, and read.pass/fail.barcode with all values of column barcode_arrangement when reads are passed/failed
        result_dict[self.add_key_to_result_dict("barcode.arrangement")] = self.dataframe_1d["barcode_arrangement"]
        self._set_result_to_dict(result_dict, "read.pass.barcode", self._series_cols_boolean_elements(self.dataframe_1d, 'barcode_arrangement', 'passes_filtering', True))
        self._set_result_to_dict(result_dict, "read.fail.barcode", self._series_cols_boolean_elements(self.dataframe_1d, 'barcode_arrangement', 'passes_filtering', False))
        
        # Get barcodes frequency by read type
        self._set_result_to_dict(result_dict,"all.read.barcoded", self.barcode_frequency(result_dict, "barcode.arrangement", 'all.read.'))
        self._set_result_to_dict(result_dict,"read.pass.barcoded", self.barcode_frequency(result_dict, "read.pass.barcode", 'read.pass.'))
        self._set_result_to_dict(result_dict,"read.fail.barcoded", self.barcode_frequency(result_dict, "read.fail.barcode", 'read.fail.'))

        read_pass_barcoded_count = result_dict.get("basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count")
        read_fail_barcoded_count = result_dict.get("basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.count")
        
        # Add key "read.pass.barcoded.frequency"
        total_reads = result_dict.get("basecaller.sequencing.summary.1d.extractor.read.count")
        result_dict[self.add_key_to_result_dict("read.pass.barcoded.frequency")] = (read_pass_barcoded_count / total_reads) * 100
        
        # Add key "read.fail.barcoded.frequency"
        result_dict[self.add_key_to_result_dict("read.fail.barcoded.frequency")] = (read_fail_barcoded_count / total_reads) * 100

                # Replaces all rows with unused barcodes (ie not in barcode_selection) in column barcode_arrangement with the 'other' value
        self.dataframe_1d.loc[~self.dataframe_1d['barcode_arrangement'].isin(
            self.barcode_selection), 'barcode_arrangement'] = 'other'

        pattern = '(\d{2})'
        self.barcode_selection.append('other')
        for index_barcode, barcode in enumerate(self.barcode_selection):
                        barcode_selected_dataframe           = self.dataframe_1d[self.dataframe_1d['barcode_arrangement'] == barcode]
                        barcode_selected_read_pass_dataframe = barcode_selected_dataframe.loc[self.dataframe_1d['passes_filtering'] == bool(True)]
                        barcode_selected_read_fail_dataframe = barcode_selected_dataframe.loc[self.dataframe_1d['passes_filtering'] == bool(False)]
                        match = re.search(pattern, barcode) # search for number of barcode used
        if match:
            barcode_name = match.group(0)
        else:
            barcode_name = barcode
            self._barcode_stats(result_dict, self.get_report_data_file_id(), barcode_selected_dataframe, barcode_selected_read_pass_dataframe, barcode_selected_read_fail_dataframe, barcode_name)
        
        #Refactor
        #index = pd.Int64Index([0, self.dataframe_1d.shape[0]])
        filtered_length_df = self.dataframe_1d.filter(items=['passes_filtering', 'sequence_length_template', 'barcode_arrangement'])
    
        barcode_selection_sequence_length_dataframe = filtered_length_df.set_index([pd.RangeIndex(start=0, stop=self.dataframe_1d.shape[0]), 'passes_filtering'], drop=True).pivot(columns="barcode_arrangement")
        barcode_selection_sequence_length_dataframe.columns.droplevel(level=0)
        col_names = [values.replace('barcode', '') for values in self.barcode_selection]
        barcode_selection_sequence_length_dataframe.columns.set_levels(col_names,
                                    level=1, inplace=True)

        barcode_selection_sequence_length_dataframe.columns = barcode_selection_sequence_length_dataframe.columns.droplevel(0)
        barcode_selection_sequence_length_dataframe.reset_index(level='passes_filtering', inplace=True)
        
        result_dict[self.add_key_to_result_dict('barcode_selection_sequence_length_dataframe')] = barcode_selection_sequence_length_dataframe
        
        result_dict[self.add_key_to_result_dict('barcode_selection_sequence_length_melted_dataframe')] = pd.melt(
            barcode_selection_sequence_length_dataframe,
            id_vars=['passes_filtering'],
            var_name="barcodes", value_name="length")

        filtered_qscore_df = self.dataframe_1d.filter(items=['passes_filtering', 'mean_qscore_template', 'barcode_arrangement'])
    
        barcode_selection_sequence_phred_dataframe = filtered_qscore_df.set_index([pd.RangeIndex(start=0, stop=self.dataframe_1d.shape[0]), 'passes_filtering'], drop=True).pivot(columns="barcode_arrangement")
        barcode_selection_sequence_phred_dataframe.columns.droplevel(level=0)
        col_names = [values.replace('barcode', '') for values in self.barcode_selection]
        barcode_selection_sequence_phred_dataframe.columns.set_levels(col_names,
                                    level=1, inplace=True)

        barcode_selection_sequence_phred_dataframe.columns = barcode_selection_sequence_phred_dataframe.columns.droplevel(0)
        barcode_selection_sequence_phred_dataframe.reset_index(level='passes_filtering', inplace=True)
        
        result_dict[self.add_key_to_result_dict('barcode_selection_sequence_phred_dataframe')] = barcode_selection_sequence_phred_dataframe
        
        result_dict[self.add_key_to_result_dict('barcode_selection_sequence_phred_melted_dataframe')] = pd.melt(
            barcode_selection_sequence_phred_dataframe,
            id_vars=['passes_filtering'],
            var_name="barcodes", value_name="qscore")


    def _barcode_stats(self, result_dict, prefix, barcode_selected_dataframe, barcode_selected_read_pass_dataframe, barcode_selected_read_fail_dataframe, barcode_name):
        
        df_dict = {'all.read.' : barcode_selected_dataframe,
                    'read.pass.' : barcode_selected_read_pass_dataframe,
                     'read.fail.' : barcode_selected_read_fail_dataframe }
            
        for df_name, df in df_dict.items(): #df_dict.items = all.read/read.pass/read.fail
            for stats_index, stats_value in df['sequence_length_template'].describe().items():
                result_dict[prefix + '.' + df_name + barcode_name + '.length.' + stats_index] = stats_value
            for stats_index, stats_value in df['mean_qscore_template'].describe().drop('count').items():
                result_dict[prefix + '.' + df_name + barcode_name + '.qscore.' + stats_index] = stats_value

        
    def graph_generation(self, result_dict):
        """
        Generation of the different graphs containing in the graph_generator module
        :return: images array containing the title and the path toward the images
        """
        images_directory = self.result_directory + '/images'
        images = list([graph_generator.read_count_histogram(result_dict, 'Read count histogram',
                                                            self.my_dpi, images_directory,
                                                            "Number of reads produced before (Fast 5 in blue) "
                                                            "and after (1D in orange) basecalling. "
                                                            "The basecalled reads are filtered with a 7.5 quality "
                                                            "score threshold in pass (1D pass in green) "
                                                            "or fail (1D fail in red) categories.")])
        images.append(graph_generator.read_length_multihistogram(result_dict, 'Read length histogram',
                                                                 self.my_dpi, images_directory,
                                                                 "Size distribution of basecalled reads (1D in orange)."
                                                                 "The basecalled reads are filtered with a 7.5 quality "
                                                                 "score threshold in pass (1D pass in green) "
                                                                 "or fail (1D fail in red) categories."))
        images.append(graph_generator.allread_number_run(result_dict, 'Yield plot of 1D read type',
                                                         self.my_dpi, images_directory,
                                                         "Yield plot of basecalled reads (1D in orange)."
                                                         " The basecalled reads are filtered with a 7.5 quality "
                                                         "score threshold in pass (1D pass in green) "
                                                         "or fail (1D fail in red) categories."))
        images.append(graph_generator.read_quality_multiboxplot(result_dict, "Read type quality boxplot",
                                                                self.my_dpi, images_directory,
                                                                "Boxplot of 1D reads (in orange) quality."
                                                                "The basecalled reads are filtered with a 7.5 quality "
                                                                "score threshold in pass (1D pass in green) "
                                                                "or fail (1D fail in red) categories."))
        images.append(graph_generator.allphred_score_frequency(result_dict,
                                                               'Mean Phred score frequency of all 1D read type',
                                                               self.my_dpi, images_directory,
                                                               "The basecalled reads are filtered with a 7.5 quality "
                                                               "score threshold in pass (1D pass in green) "
                                                               "or fail (1D fail in red) categories."))
        channel_count = self.channel
        total_number_reads_per_pore = pd.value_counts(channel_count)
        images.append(graph_generator.plot_performance(total_number_reads_per_pore, 'Channel occupancy of the flowcell',
                                                       self.my_dpi, images_directory,
                                                       "Number of reads sequenced per pore channel."))

        images.append(graph_generator.all_scatterplot(result_dict, 'Mean Phred score function of 1D read length',
                                                      self.my_dpi, images_directory,
                                                      "The Mean Phred score varies according to the read length."
                                                      "The basecalled reads are filtered with a 7.5 quality "
                                                      "score threshold in pass (1D pass in green) "
                                                      "or fail (1D fail in red) categories."))
        if self.is_barcode:
            images.append(graph_generator.barcode_percentage_pie_chart_pass(result_dict,
                                                                            '1D pass reads percentage of different '
                                                                            'barcodes', self.barcode_selection,
                                                                            self.my_dpi, images_directory,
                                                                            "1D pass read distribution per barcode."))

            images.append(graph_generator.barcode_percentage_pie_chart_fail(result_dict,
                                                                            '1D fail reads percentage of different '
                                                                            'barcodes', self.barcode_selection,
                                                                            self.my_dpi, images_directory,
                                                                            "1D fail read distribution per barcode."))

            images.append(graph_generator.barcode_length_boxplot(result_dict,
                                                                 '1D reads size distribution for each barcode',
                                                                 self.my_dpi, images_directory,
                                                                 "Read length boxplot per barcode of pass (in green) "
                                                                 "and fail (in red) 1D reads."))

            images.append(graph_generator.barcoded_phred_score_frequency(result_dict,
                                                                         '1D reads Mean Phred score distribution '
                                                                         'for each barcode',
                                                                         self.my_dpi, images_directory,
                                                                         "Read Mean Phred score boxplot per barcode of "
                                                                         "pass (in green) and fail (in red) 1D reads."))
        return images


    def clean(self, result_dict):
        """
        Removing dictionary entries that will not be kept in the report.data file
        :return:
        """
        #TODO: supprimer les clefs non utilisées et enlever les commentaires
        keys = ["sequence.length", "passes.filtering", "read.pass.length", "read.fail.length",
                "start.time.sorted", "read.pass.sorted", "read.fail.sorted",
                "mean.qscore", "read.pass.qscore", "read.fail.qscore",
                'channel.occupancy.statistics', #TODO: delete key 'channel.occupancy.statistics' if not useful after refactoring
                'all.read.qscore', 'all.read.length', #TODO: delete keys 'all.read.qscore' & 'all.read.length' if not useful after refactoring
                "barcode.arrangement", "read.pass.barcode", "read.fail.barcode",
                'barcode_selection_sequence_length_dataframe', #utilisé dans graph_generator.py
                'barcode_selection_sequence_length_melted_dataframe', #utilisé dans graph_generator.py
                'barcode_selection_sequence_phred_dataframe', #utilisé dans graph_generator.py
                'barcode_selection_sequence_phred_melted_dataframe', #utilisé dans graph_generator.py
                "all.read.barcoded", #TODO: delete key after refactoring
                "read.pass.barcoded", #utilisé dans graph_generator.py
                "read.fail.barcoded" #utilisé dans graph_generator.py
        ]

        key_list = []
        for key in keys:
            key_list.append(self.add_key_to_result_dict(key))
        result_dict['unwritten.keys'].extend(key_list)


    def _occupancy_channel(self):
        """
          Statistics about the channels of the flowcell
        :return: pd.Series object containing statistics about the channel occupancy
        """
          # Returns a pd.Series containing counts of unique values, by descending order
        total_reads_per_channel = pd.value_counts(self.channel)
        return pd.DataFrame.describe(total_reads_per_channel)

    
    def _load_sequencing_summary_data(self):
        """
        Load sequencing summary data frame.
        :return: a Pandas DataFrame object
        """
        # Initialization
        files = self.sequencing_summary_files

        summary_dataframe = None
        barcode_dataframe = None

        sequencing_summary_columns = ['read_id', 'channel', 'start_time', 'duration',
        'num_events', 'num_events_template', 'passes_filtering', 'sequence_length_template', 'mean_qscore_template']

        sequencing_summary_datatypes = {
        'read_id' : object,
        'channel': np.int16,
        'start_time': np.float,
        'duration': np.float,
        'num_events': np.int16,
        'num_events_template' : np.int16,
        'passes_filtering': np.bool,
        'sequence_length_template': np.int16,
        'mean_qscore_template': np.float}

        barcoding_summary_columns = ['read_id', 'barcode_arrangement']

        barcoding_summary_datatypes = {
            'read_id' : object,
            'barcode_arrangement' : object
            }

        # If only one file and it's a sequencing_summary.txt
        if (len(files) == 1 and self._is_sequencing_summary_file(files[0])):
            return pd.read_csv(files[0], sep="\t", usecols=sequencing_summary_columns, dtype=sequencing_summary_datatypes)
        
        # If multiple files, check if there's a barcoding one and a sequencing one :
        for f in files:
            # check if barcoding_summary is in files
            if self._is_barcode_file(f):
                dataframe = pd.read_csv(f, sep="\t", usecols=barcoding_summary_columns, dtype=barcoding_summary_datatypes)
                if barcode_dataframe is None:
                    barcode_dataframe = dataframe
                # if a barcoding file has already been read, append the 2 dataframes
                else:
                    barcode_dataframe = barcode_dataframe.append(dataframe, ignore_index=True)

            # when barcoding_summary is not here, check if the files contain a sequencing_summary
            else:
                if self._is_sequencing_summary_file(f):
                    dataframe = pd.read_csv(f, sep="\t", usecols=sequencing_summary_columns, dtype=sequencing_summary_datatypes)
                    if summary_dataframe is None:
                        summary_dataframe = dataframe
                    else:
                        summary_dataframe = summary_dataframe.append(dataframe, ignore_index=True)
                        
        #TODO: Condition jamais atteinte dans les tests unitaires, car si aucun fichier, les méthodes _is_barcoding/sequencing_summary renvoient FileNotFoundError: [Errno 2] No such file or directory: ''
        if (summary_dataframe is None and barcode_dataframe is None):
            raise ValueError("Sequencing summary file not found nor barcoding summary file(s)")

        elif barcode_dataframe is None:
            return summary_dataframe
        else:
            dataframes_merged = pd.merge(summary_dataframe, barcode_dataframe, on='read_id', how = 'left')
            # delete column read_id after merging
            del dataframes_merged['read_id']

            return dataframes_merged


    def _is_barcode_file(self, filename):
        """
        Check if input is a barcoding summary file i.e. has the column barcode_arrangement
        :param filename: path of the file to test
        :return: True if the filename is a barcoding summary file
        """
        with open(filename, 'r') as f:
            header = f.readline()
            return "barcode_arrangement" in header


    def _is_sequencing_summary_file(self, filename):
        """
        Check if input is a sequencing summary file i.e. first word is "filename"
        :param filename: path of the file to test
        :return: True if the file is indeed a sequencing summary file
        """
        with open(filename, 'r') as f:
            header = f.readline()
        return header.startswith("filename")
