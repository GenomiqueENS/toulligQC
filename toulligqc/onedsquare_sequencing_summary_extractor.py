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
# First author: Bérengère Laffay
# Maintainer: Bérengère Laffay
# Since version 0.6

# Extraction of statistics from 1dsqr_sequencing_summary.txt file (1Dsquare chemistry)

import pandas as pd
import sys
from toulligqc import graph_generator
import numpy as np
import re
import os.path


class OneDSquareSequencingSummaryExtractor:
    """
    Extraction of statistics from 1dsqr_sequencing_summary.txt file and graph generation
    """
    def __init__(self, config_dictionary):

        self.config_dictionary = config_dictionary
        self.sequencing_summary_source = self.config_dictionary['sequencing_summary_source']
        self.sequencing_1dsqr_summary_source = self.config_dictionary['sequencing_summary_1dsqr_source']
        self.result_directory = config_dictionary['result_directory']
        self.my_dpi = int(self.config_dictionary['dpi'])

        # Attribute initialized in the init() method
        self.dataframe_1d = None
        self.channel = None
        self.dataframe_1dsqr = None
        self.barcode_selection = None

        if config_dictionary['barcoding'] == 'True':
            self.is_barcode = True
        else:
            self.is_barcode = False

        self.sequencing_summary_files = self.sequencing_summary_source.split('\t')

        if len(self.sequencing_summary_files) == 1:
            if os.path.isdir(self.sequencing_summary_source):
                self.sequencing_summary_files = [self.sequencing_summary_source + "/sequencing_summary.txt"]
            else:
                self.sequencing_summary_files = [self.sequencing_summary_source]

        self.sequencing_1dsqr_summary_files = self.sequencing_1dsqr_summary_source.split('\t')

        if len(self.sequencing_1dsqr_summary_files) == 1:
            if os.path.isdir(self.sequencing_1dsqr_summary_source):
                self.sequencing_1dsqr_summary_files = [self.sequencing_1dsqr_summary_source + "/sequencing_1dsq_summary.txt"]
            else:
                self.sequencing_1dsqr_summary_files = [self.sequencing_1dsqr_summary_source]

    def check_conf(self):
        """Configuration checking"""

        # Check if a sequencing summary file has been defined
        if len(self.sequencing_summary_files) == 0:
            return False, "No sequencing summary file defined"

        # Check if sequencing summary files exists
        for f in self.sequencing_summary_files:
            if not os.path.isfile(f):
                return False, "Sequencing summary file does not exists: " + f

        # Check if a 1d2 sequencing summary file has been defined
        if len(self.sequencing_1dsqr_summary_files) == 0:
            return False, "No 1dsqr sequencing summary file defined"

        # Check if 1d2 sequencing summary files exists
        for f in self.sequencing_1dsqr_summary_files:
            if not os.path.isfile(f):
                return False, "Sequencing summary file does not exists: " + f

        return True, ""

    def init(self):
        """
        Initialisation
        :return:
        """

        # Panda's object for 1d_summary
        self.dataframe_1d = pd.read_csv(self.sequencing_summary_files[0], sep="\t")
        self.channel = self.dataframe_1d['channel']
        self.dataframe_1d = self.dataframe_1d.replace([np.inf, -np.inf], 0)
        self.dataframe_1d = self.dataframe_1d[self.dataframe_1d['num_events'] != 0]
        self.dataframe_1d["Yield"] = sum(self.dataframe_1d['sequence_length_template'])

        # Panda's object for 1dsqr_summary
        self.dataframe_1dsqr = self._load_sequencing_summary_data()

        if self.is_barcode:

            self.barcode_selection = self.config_dictionary['barcode_selection']

            # try:
            #     self.dataframe_1dsqr.loc[~self.dataframe_1dsqr['barcode_arrangement'].isin(
            #         self.barcode_selection), 'barcode_arrangement'] = 'unclassified'
            # except ValueError:
            #     sys.exit('No barcode found in sequencing summary file')

    @staticmethod
    def get_name():
        """
        Get the name of the extractor.
        :return: the name of the extractor
        """
        return 'Basecaller 1d square sequencing summary'

    @staticmethod
    def get_report_data_file_id():
        """
        Get the report.data id of the extractor.
        :return: the report.data id
        """
        return 'basecaller.sequencing.summary.1dsqr.extractor'

    def add_key_to_result_dict(self, key):
        """
        :param key:
        :return:
        """
        return '{0}.{1}'.format(self.get_report_data_file_id(), key)

    @staticmethod
    def describe_dict(result_dict, attribute):
        """
        :param result_dict:
        :param attribute:
        :return:
        """
        dictionary = pd.Series.describe(result_dict[attribute])
        for key in dict(dictionary):
            result_dict[attribute + '.' + key] = dictionary[key]

    def barcode_frequency(self, result_dict, entry, prefix=''):
        """
        Adds barcode frequency to the result_dict dictionary
        :param result_dict:
        :param entry:
        :param prefix: key prefix
        :return: result_dict dictionary filled and barcodes frequency in pandas series type
        """
        barcode_count = result_dict[self.add_key_to_result_dict(entry)].value_counts()
        count_sorted = barcode_count.sort_index()[self.barcode_selection]
        result_dict[self.add_key_to_result_dict(prefix + 'barcoded.count')] = sum(count_sorted.drop("unclassified"))
        result_dict[self.add_key_to_result_dict(prefix + "with.other.barcodes.count")] = (sum(barcode_count)-sum(count_sorted))
        other_barcode_count = pd.Series(result_dict[self.add_key_to_result_dict(prefix + "with.other.barcodes.count")], index=['other'])
        count_sorted = count_sorted.append(other_barcode_count).sort_index()

        for key in dict(count_sorted):
            result_dict[self.add_key_to_result_dict(prefix) + key + ".frequency"] = \
                count_sorted[key]*100/sum(count_sorted)
        return count_sorted

    def extract(self, result_dict):
        """
        :param result_dict:
        :return:
        """
        #
        # Extract from 1D summary source
        #

        # Read count
        result_dict['basecaller.sequencing.summary.1d.extractor.fastq.entries'] = \
            len(self.dataframe_1d['num_events'])

        # 1D pass information
        result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"] = \
            len(self.dataframe_1d.loc[self.dataframe_1d['passes_filtering'] == bool(True)])
        result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.length"] = \
            self.dataframe_1d.sequence_length_template.loc[self.dataframe_1d['passes_filtering'] == bool(True)]
        result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.sorted"] = \
            sorted(self.dataframe_1d.start_time.loc[self.dataframe_1d['passes_filtering'] == bool(True)] / 3600)
        result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.qscore"] = \
            self.dataframe_1d.mean_qscore_template.loc[self.dataframe_1d['passes_filtering'] == bool(True)]

        # 1D fail information
        result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"] = \
            len(self.dataframe_1d.loc[self.dataframe_1d['passes_filtering'] == bool(False)])
        result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.length"] = \
            self.dataframe_1d.sequence_length_template.loc[self.dataframe_1d['passes_filtering'] == bool(False)]
        result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.sorted"] = \
            sorted(self.dataframe_1d.start_time.loc[self.dataframe_1d['passes_filtering'] == bool(False)] / 3600)
        result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.qscore"] = \
            self.dataframe_1d.mean_qscore_template.loc[self.dataframe_1d['passes_filtering'] == bool(False)]

        # The field  "num_called_template" has been renamed "num_events_template" in Guppy
        if "num_called_template" in self.dataframe_1d:
            num_called_template_field = "num_called_template"
        else:
            num_called_template_field = "num_events_template"

        result_dict['basecaller.sequencing.summary.1d.extractor.read.count'] = \
            len(self.dataframe_1d[self.dataframe_1d[num_called_template_field] != 0])
        result_dict["basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.count"] = \
            len(self.dataframe_1d[self.dataframe_1d['sequence_length_template'] == 0])

        # Read proportion
        result_dict["basecaller.sequencing.summary.1d.extractor.fastq.entries.ratio"] = \
            result_dict['basecaller.sequencing.summary.1d.extractor.fastq.entries'] / \
            result_dict['basecaller.sequencing.summary.1d.extractor.fastq.entries']

        result_dict["basecaller.sequencing.summary.1d.extractor.read.count.ratio"] = \
            result_dict["basecaller.sequencing.summary.1d.extractor.read.count"] / \
            result_dict["basecaller.sequencing.summary.1d.extractor.read.count"]

        result_dict["basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.ratio"] = \
            result_dict["basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.count"] / \
            result_dict["basecaller.sequencing.summary.1d.extractor.read.count"]

        result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.ratio"] = \
            result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"] / \
            result_dict["basecaller.sequencing.summary.1d.extractor.read.count"]

        result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.ratio"] = \
            result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"] / \
            result_dict["basecaller.sequencing.summary.1d.extractor.read.count"]

        result_dict["basecaller.sequencing.summary.1d.extractor.fastq.entries.frequency"] = \
            result_dict['basecaller.sequencing.summary.1d.extractor.fastq.entries'] / \
            result_dict['basecaller.sequencing.summary.1d.extractor.fastq.entries']*100

        result_dict["basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.frequency"] = \
            result_dict["basecaller.sequencing.summary.1d.extractor.read.with.length.equal.zero.count"] / \
            result_dict["basecaller.sequencing.summary.1d.extractor.read.count"]*100

        result_dict["basecaller.sequencing.summary.1d.extractor.read.count.frequency"] = \
            result_dict["basecaller.sequencing.summary.1d.extractor.read.count"] / \
            result_dict["basecaller.sequencing.summary.1d.extractor.read.count"]*100

        result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.frequency"] = \
            result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"] / \
            result_dict["basecaller.sequencing.summary.1d.extractor.read.count"]*100

        result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.frequency"] = \
            result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"] / \
            result_dict["basecaller.sequencing.summary.1d.extractor.read.count"]*100

        # Read length information
        result_dict["basecaller.sequencing.summary.1d.extractor.sequence.length"] = \
            self.dataframe_1d.sequence_length_template[self.dataframe_1d[num_called_template_field] != 0]

        result_dict["basecaller.sequencing.summary.1d.extractor.passes.filtering"] = \
            self.dataframe_1d['passes_filtering']

        # Yield
        result_dict["basecaller.sequencing.summary.1d.extractor.yield"] = \
            sum(self.dataframe_1d['sequence_length_template'] / 1000000000)

        result_dict["basecaller.sequencing.summary.1d.extractor.start.time.sorted"] = \
            sorted(sorted(self.dataframe_1d['start_time'] / 3600))

        result_dict["basecaller.sequencing.summary.1d.extractor.run.time"] = \
            int(max(result_dict["basecaller.sequencing.summary.1d.extractor.start.time.sorted"]))

        # Qscore information
        result_dict["basecaller.sequencing.summary.1d.extractor.mean.qscore"] = self.dataframe_1d.loc[:, "mean_qscore_template"]

        # Channel occupancy information
        result_dict["basecaller.sequencing.summary.1d.extractor.channel.occupancy.statistics"] = self._occupancy_channel()
        channel_occupancy_statistics = result_dict['basecaller.sequencing.summary.1d.extractor.channel.occupancy.statistics']

        for index, value in channel_occupancy_statistics.iteritems():
            result_dict['basecaller.sequencing.summary.1d.extractor.channel.occupancy.statistics.' + index] = value

        # Length's statistic information provided in the result_dict
        result_dict["basecaller.sequencing.summary.1d.extractor.all.read.length"] = \
            self.dataframe_1d['sequence_length_template'].describe()

        for index, value in result_dict["basecaller.sequencing.summary.1d.extractor.all.read.length"].iteritems():
            result_dict["basecaller.sequencing.summary.1d.extractor.all.read.length." + index] = value
        self.describe_dict(result_dict, "basecaller.sequencing.summary.1d.extractor.read.pass.length")
        self.describe_dict(result_dict, "basecaller.sequencing.summary.1d.extractor.read.fail.length")

        # Qscore's statistic information provided in the result_dict
        result_dict['basecaller.sequencing.summary.1d.extractor.all.read.qscore'] = \
            pd.DataFrame.describe(self.dataframe_1d['mean_qscore_template']).drop("count")

        for index, value in result_dict["basecaller.sequencing.summary.1d.extractor.all.read.qscore"].iteritems():
            result_dict["basecaller.sequencing.summary.1d.extractor.all.read.qscore." + index] = value

        self.describe_dict(result_dict, "basecaller.sequencing.summary.1d.extractor.read.pass.qscore")
        self.describe_dict(result_dict, "basecaller.sequencing.summary.1d.extractor.read.fail.qscore")

        #
        # Extract from 1dsqr sequencing summary
        #

        # The field  "sequence_length_2d" has been renamed "sequence_length" in Guppy
        if "sequence_length_2d" in self.dataframe_1dsqr.columns:
            sequence_length_field = "sequence_length_2d"
        else:
            sequence_length_field = "sequence_length"

        # The field  "mean_qscore_2d" has been renamed "mean_qscore" in Guppy
        if "mean_qscore_2d" in self.dataframe_1dsqr.columns:
            mean_qscore_field = "mean_qscore_2d"
        else:
            mean_qscore_field = "mean_qscore"

        # Read count
        result_dict[self.add_key_to_result_dict('read.count')] = \
            len(self.dataframe_1dsqr['passes_filtering'])

        # 1Dsquare pass information
        result_dict[self.add_key_to_result_dict('read.pass.count')] = \
            len(self.dataframe_1dsqr.loc[self.dataframe_1dsqr['passes_filtering'] == bool(True)])
        result_dict[self.add_key_to_result_dict('read.pass.length')] = \
            self.dataframe_1dsqr[sequence_length_field].loc[self.dataframe_1dsqr['passes_filtering'] == bool(True)]
        result_dict[self.add_key_to_result_dict('read.pass.qscore')] = \
            self.dataframe_1dsqr[mean_qscore_field].loc[self.dataframe_1dsqr['passes_filtering'] == bool(True)]

        # 1Dsquare fail information
        result_dict[self.add_key_to_result_dict('read.fail.count')] = \
            len(self.dataframe_1dsqr.loc[self.dataframe_1dsqr['passes_filtering'] == bool(False)])
        result_dict[self.add_key_to_result_dict('read.fail.length')] = \
            self.dataframe_1dsqr[sequence_length_field].loc[self.dataframe_1dsqr['passes_filtering'] == bool(False)]
        result_dict[self.add_key_to_result_dict('read.fail.qscore')] = \
            self.dataframe_1dsqr[mean_qscore_field].loc[self.dataframe_1dsqr['passes_filtering'] == bool(False)]

        # Read count proportion
        result_dict[self.add_key_to_result_dict("read.count.ratio")] = \
            result_dict[self.add_key_to_result_dict('read.count')] / \
            result_dict[self.add_key_to_result_dict('read.count')]

        result_dict[self.add_key_to_result_dict("read.pass.ratio")] = \
            result_dict[self.add_key_to_result_dict('read.pass.count')] / \
            result_dict[self.add_key_to_result_dict('read.count')]

        result_dict[self.add_key_to_result_dict("read.fail.ratio")] = \
            result_dict[self.add_key_to_result_dict('read.fail.count')] / \
            result_dict[self.add_key_to_result_dict('read.count')]

        result_dict[self.add_key_to_result_dict("read.count.frequency")] = \
            result_dict[self.add_key_to_result_dict('read.count')] / \
            result_dict[self.add_key_to_result_dict('read.count')]*100

        result_dict[self.add_key_to_result_dict("read.pass.frequency")] = \
            result_dict[self.add_key_to_result_dict('read.pass.count')] / \
            result_dict[self.add_key_to_result_dict('read.count')]*100

        result_dict[self.add_key_to_result_dict("read.fail.frequency")] = \
            result_dict[self.add_key_to_result_dict('read.fail.count')] / \
            result_dict[self.add_key_to_result_dict('read.count')]*100

        # Read length information
        result_dict[self.add_key_to_result_dict('sequence.length')] = \
            self.dataframe_1dsqr.loc[:, sequence_length_field]

        result_dict[self.add_key_to_result_dict("passes.filtering")] = self.dataframe_1dsqr['passes_filtering']

        # Qscore information
        result_dict[self.add_key_to_result_dict('mean.qscore')] = self.dataframe_1dsqr.loc[:, mean_qscore_field]

        # Length's statistic information provided in the result_dict
        result_dict[self.add_key_to_result_dict('all.read.length')] = \
            self.dataframe_1dsqr[sequence_length_field].describe()

        for index, value in result_dict[self.add_key_to_result_dict('all.read.length')].iteritems():
            result_dict[self.add_key_to_result_dict('all.read.length.') + index] = value
        self.describe_dict(result_dict, self.add_key_to_result_dict("read.pass.length"))
        self.describe_dict(result_dict, self.add_key_to_result_dict("read.fail.length"))

        # Qscore's statistic information provided in the result_dict
        result_dict[self.add_key_to_result_dict('all.read.qscore')] = \
            pd.DataFrame.describe(self.dataframe_1dsqr[mean_qscore_field]).drop("count")

        for index, value in result_dict[self.add_key_to_result_dict('all.read.qscore')].iteritems():
            result_dict[self.add_key_to_result_dict('all.read.qscore.') + index] = value
        self.describe_dict(result_dict, self.add_key_to_result_dict("read.pass.qscore"))
        self.describe_dict(result_dict, self.add_key_to_result_dict("read.fail.qscore"))

        # In case of barcoded samples

        if self.is_barcode:
            self.barcode_selection.append('unclassified')

            result_dict[self.add_key_to_result_dict("barcode.arrangement")] = \
                self.dataframe_1dsqr["barcode_arrangement"]

            result_dict[self.add_key_to_result_dict("read.pass.barcode")] = \
                self.dataframe_1dsqr.barcode_arrangement.loc[
                    self.dataframe_1dsqr['passes_filtering'] == bool(True)]

            result_dict[self.add_key_to_result_dict("read.fail.barcode")] = \
                self.dataframe_1dsqr.barcode_arrangement.loc[
                    self.dataframe_1dsqr['passes_filtering'] == bool(False)]

            # Get barcodes frequency by read type
            result_dict[self.add_key_to_result_dict("all.read.barcoded")] = self.barcode_frequency(result_dict, "barcode.arrangement", 'all.read.')
            result_dict[self.add_key_to_result_dict("read.pass.barcoded")] = self.barcode_frequency(result_dict, "read.pass.barcode", 'read.pass.')
            result_dict[self.add_key_to_result_dict("read.fail.barcoded")] = self.barcode_frequency(result_dict, "read.fail.barcode", 'read.fail.')

            result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.barcoded.frequency"] = result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.pass.barcoded.count"]/result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.count"] *100
            result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.barcoded.frequency"] = result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.fail.barcoded.count"]/result_dict["basecaller.sequencing.summary.1dsqr.extractor.read.count"] *100

            self.dataframe_1dsqr.loc[~self.dataframe_1dsqr['barcode_arrangement'].isin(
                    self.barcode_selection), 'barcode_arrangement'] = 'other'

            result_dict[self.add_key_to_result_dict("passes.filtering")] = self.dataframe_1dsqr['passes_filtering']

            pattern = '(\d{2})'
            length = {'passes_filtering': result_dict[self.add_key_to_result_dict("passes.filtering")]}
            phred = {'passes_filtering': result_dict[self.add_key_to_result_dict("passes.filtering")]}

            self.barcode_selection.append('other')

            for index_barcode, barcode in enumerate(self.barcode_selection):
                barcode_selected_dataframe = \
                    self.dataframe_1dsqr[self.dataframe_1dsqr['barcode_arrangement'] == barcode]

                barcode_selected_read_pass_dataframe = \
                    barcode_selected_dataframe.loc[barcode_selected_dataframe['passes_filtering'] == bool(True)]

                barcode_selected_read_fail_dataframe = \
                    barcode_selected_dataframe.loc[barcode_selected_dataframe['passes_filtering'] == bool(False)]

                match = re.search(pattern, barcode)
                if match:
                    length[match.group(0)] = barcode_selected_dataframe[sequence_length_field]
                    phred[match.group(0)] = barcode_selected_dataframe[mean_qscore_field]

                    for index, value in barcode_selected_dataframe[sequence_length_field]\
                            .describe().iteritems():
                        result_dict[self.add_key_to_result_dict('all.read.') + barcode + '.length.' + index] = value

                    for index, value in barcode_selected_read_pass_dataframe[sequence_length_field]\
                            .describe().iteritems():
                        result_dict[self.add_key_to_result_dict('read.pass.') + barcode + '.length.' + index] = value

                    for index, value in barcode_selected_read_fail_dataframe[sequence_length_field]\
                            .describe().iteritems():
                        result_dict[self.add_key_to_result_dict('read.fail.') + barcode + '.length.' + index] = value

                    for index, value in barcode_selected_dataframe[mean_qscore_field]\
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('all.read.') + barcode + '.qscore.' + index] = value

                    for index, value in barcode_selected_read_pass_dataframe[mean_qscore_field]\
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('read.pass.') + barcode + '.qscore.' + index] = value

                    for index, value in barcode_selected_read_fail_dataframe[mean_qscore_field]\
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('read.fail.') + barcode + '.qscore.' + index] = value

                if barcode == 'unclassified':
                    length['Unclassified'] = barcode_selected_dataframe[sequence_length_field]
                    phred['Unclassified'] = barcode_selected_dataframe[mean_qscore_field]

                    for index, value in barcode_selected_dataframe[sequence_length_field]\
                            .describe().iteritems():
                        result_dict[self.add_key_to_result_dict('all.read.unclassified.length.') + index] = value

                    for index, value in barcode_selected_read_pass_dataframe[sequence_length_field]\
                            .describe().iteritems():
                        result_dict[self.add_key_to_result_dict('read.pass.unclassified.length.') + index] = value

                    for index, value in barcode_selected_read_fail_dataframe[sequence_length_field]\
                            .describe().iteritems():
                        result_dict[self.add_key_to_result_dict('read.fail.unclassified.length.') + index] = value

                    for index, value in barcode_selected_dataframe[mean_qscore_field]\
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('all.read.unclassified.qscore.') + index] = value

                    for index, value in barcode_selected_read_pass_dataframe[mean_qscore_field]\
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('read.pass.unclassified.qscore.') + index] = value

                    for index, value in barcode_selected_read_fail_dataframe[mean_qscore_field]\
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('read.fail.unclassified.qscore.') + index] = value

                if barcode == 'other':

                    length['Other Barcodes'] = barcode_selected_dataframe[sequence_length_field]
                    phred['Other Barcodes'] = barcode_selected_dataframe[mean_qscore_field]

                    for index, value in barcode_selected_dataframe[sequence_length_field] \
                            .describe().iteritems():
                        result_dict[self.add_key_to_result_dict('all.read.with.other.barcodes.length.') + index] = value

                    for index, value in barcode_selected_read_pass_dataframe[sequence_length_field] \
                            .describe().iteritems():
                        result_dict[self.add_key_to_result_dict('read.pass.with.other.barcodes.length.') + index] = value

                    for index, value in barcode_selected_read_fail_dataframe[sequence_length_field] \
                            .describe().iteritems():
                        result_dict[self.add_key_to_result_dict('read.fail.with.other.barcodes.length.') + index] = value

                    for index, value in barcode_selected_dataframe[mean_qscore_field] \
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('all.read.with.other.barcodes.qscore.') + index] = value

                    for index, value in barcode_selected_read_pass_dataframe[mean_qscore_field] \
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('read.pass.with.other.barcodes.qscore.') + index] = value

                    for index, value in barcode_selected_read_fail_dataframe[mean_qscore_field] \
                            .describe().drop('count').iteritems():
                        result_dict[self.add_key_to_result_dict('read.fail.with.other.barcodes.qscore.') + index] = value

            # Provide statistic per barcode in the result_dict dictionary

            result_dict[self.add_key_to_result_dict('barcode_selection_sequence_length_dataframe')] = \
                pd.DataFrame(dict([(k, pd.Series(v)) for k, v in length.items()]))

            result_dict[self.add_key_to_result_dict('barcode_selection_sequence_length_melted_dataframe')] = \
                pd.melt(result_dict[self.add_key_to_result_dict('barcode_selection_sequence_length_dataframe')],
                        id_vars=['passes_filtering'], var_name="barcodes", value_name="length")

            result_dict[self.add_key_to_result_dict('barcode_selection_sequence_phred_dataframe')] = \
                pd.DataFrame(dict([(k, pd.Series(v)) for k, v in phred.items()]))

            result_dict[self.add_key_to_result_dict('barcode_selection_sequence_phred_melted_dataframe')] = \
                pd.melt(result_dict[self.add_key_to_result_dict('barcode_selection_sequence_phred_dataframe')],
                        id_vars=['passes_filtering'], var_name="barcodes", value_name="qscore")

            length.clear()
            phred.clear()

    def graph_generation(self, result_dict):
        """
        Generation of the differents graphs containing in the graph_generator module
        :return: images array containing the title and the path toward the images
        """
        images_directory = self.result_directory + '/images/'
        images = list([graph_generator.read_count_histogram(result_dict, '1D read count histogram',
                                                            self.my_dpi, images_directory,
                                                            "Number of reads produced before (Fast 5 in blue) and after"
                                                            " (1D in orange) basecalling. The basecalled reads are "
                                                            "filtered with a 7.5 quality score threshold in pass "
                                                            "(1D pass in green) or fail (1D fail in red) categories.")])

        images.append(graph_generator.dsqr_read_count_histogram(result_dict, "1Dsquare read count histogram",
                                                                self.my_dpi, images_directory,
                                                                "Number of reads produced basecalled (1D in orange) and"
                                                                " 1Dsquare reads (in gold). The 1Dsquare reads are "
                                                                "filtered with a 7.5 quality score threshold in pass "
                                                                "(1Dsquare pass in green) or fail "
                                                                "(1Dsquare fail in red) categories."))

        images.append(graph_generator.read_length_multihistogram(result_dict, '1D read size histogram',
                                                                 self.my_dpi, images_directory,
                                                                 "Size distribution of basecalled reads (1D in orange)."
                                                                 " The basecalled reads are filtered with a 7.5 quality"
                                                                 " score threshold in pass (1D pass in green) or fail "
                                                                 "(1D fail in red) categories."))

        images.append(graph_generator.dsqr_read_length_multihistogram(result_dict, '1Dsquare read size histogram',
                                                                      self.my_dpi, images_directory,
                                                                      "Size distribution of basecalled reads "
                                                                      "(1D in orange) and 1Dsquare reads (in gold). "
                                                                      "The 1Dsquare reads are filtered with a 7.5 "
                                                                      "quality score threshold in pass "
                                                                      "(1Dsquare pass in green) or fail "
                                                                      "(1Dsquare fail in red) categories."))

        images.append(graph_generator.allread_number_run(result_dict, 'Yield plot of 1D read type',
                                                         self.my_dpi, images_directory,
                                                         "Yield plot of basecalled reads (1D in orange). "
                                                         "The basecalled reads are filtered with a 7.5 quality score "
                                                         "threshold in pass (1D pass in green) or fail "
                                                         "(1D fail in red) categories."))

        images.append(graph_generator.read_quality_multiboxplot(result_dict, "1D reads quality boxplot",
                                                                self.my_dpi, images_directory,
                                                                "Boxplot of 1D reads (in orange) quality. "
                                                                "The basecalled reads are filtered with a 7.5 quality "
                                                                "score threshold in pass (1D pass in green) or fail "
                                                                "(1D fail in red) categories."))

        images.append(graph_generator.dsqr_read_quality_multiboxplot(result_dict, "1Dsquare reads quality boxplot",
                                                                     self.my_dpi, images_directory,
                                                                     "Boxplot of 1D (in orange) and 1Dsquare (in gold) "
                                                                     "reads quality. The 1Dsquare reads are filtered "
                                                                     "with a 7.5 quality score threshold in pass "
                                                                     "(1Dsquare pass in green) or fail (1Dsquare "
                                                                     "fail in red) categories."))

        images.append(graph_generator.allphred_score_frequency(result_dict, 'Mean Phred score frequency of '
                                                                            '1D read type',
                                                               self.my_dpi, images_directory,
                                                               "The basecalled reads are filtered with a 7.5 quality "
                                                               "score threshold in pass (1D pass in green) or fail "
                                                               "(1D fail in red) categories."))

        images.append(graph_generator.dsqr_allphred_score_frequency(result_dict, "Mean Phred score frequency of "
                                                                                 "1Dsquare read type",
                                                                    self.my_dpi, images_directory,
                                                                    "The 1Dsquare reads are filtered with a 7.5 "
                                                                    "quality score threshold in pass (1Dsquare pass "
                                                                    "in green) or fail (1Dsquare fail in red) "
                                                                    "categories."))

        channel_count = self.channel
        total_number_reads_per_pore = pd.value_counts(channel_count)
        images.append(graph_generator.plot_performance(total_number_reads_per_pore, 'Channel occupancy of the flowcell',
                                                       self.my_dpi, images_directory,
                                                       "Number of reads sequenced per pore channel."))

        images.append(graph_generator.all_scatterplot(result_dict, 'Mean Phred score function of 1D read length',
                                                      self.my_dpi, images_directory,
                                                      "The Mean Phred score varies according to the read length. "
                                                      "The basecalled reads are filtered with a 7.5 quality score "
                                                      "threshold in pass (1D pass in green) or fail "
                                                      "(1D fail in red) categories."))

        images.append(graph_generator.scatterplot_1dsqr(result_dict, "Mean Phred score function of 1Dsquare "
                                                                     "read length",
                                                        self.my_dpi, images_directory,
                                                        "The Mean Phred score varies according to the read length. "
                                                        "The 1Dsquare reads are filtered with a 7.5 quality score "
                                                        "threshold in pass (1Dsquare pass in green) or fail "
                                                        "(1Dsquare fail in red) categories."))

        if self.is_barcode:
            images.append(graph_generator.barcode_percentage_pie_chart_1dsqr_pass(result_dict,
                                                                                  "1Dsquare pass reads percentage of "
                                                                                  "different barcodes",
                                                                                  self.barcode_selection, self.my_dpi,
                                                                                  images_directory,
                                                                                  "1Dsquare pass reads distribution "
                                                                                  "per barcode."))

            images.append(graph_generator.barcode_percentage_pie_chart_1dsqr_fail(result_dict,
                                                                                  "1Dsquare fail reads percentage of "
                                                                                  "different barcodes",
                                                                                  self.barcode_selection, self.my_dpi,
                                                                                  images_directory,
                                                                                  "1Dsquare fail reads distribution "
                                                                                  "per barcode."))

            images.append(graph_generator.barcode_length_boxplot_1dsqr(result_dict,
                                                                       "1Dsquare read size distribution for "
                                                                       "each barcode",
                                                                       self.my_dpi, images_directory,
                                                                       "Read length boxplot per barcode of pass "
                                                                       "(in green) and fail (in red) 1Dsquare reads."))

            images.append(graph_generator.barcoded_phred_score_frequency_1dsqr(result_dict,
                                                                               "1Dsquare read phred score distribution "
                                                                               "for each barcode",
                                                                               self.my_dpi, images_directory,
                                                                               "Read Mean Phred score boxplot per "
                                                                               "barcode of pass (in green) and fail "
                                                                               "(in red) 1Dsquare reads."))
        return images

    def clean(self, result_dict):
        """
        Removing dictionary entries that will not be kept in the report.data file
        :param result_dict:
        :return:
        """

        keys = ['sequence.length', 'passes.filtering', 'read.pass.length', 'read.fail.length',
                'mean.qscore', 'read.pass.qscore', 'read.fail.qscore',
                'all.read.qscore', 'all.read.length',
                "barcode.arrangement", "read.pass.barcode", "read.fail.barcode",
                'barcode_selection_sequence_length_dataframe', 'barcode_selection_sequence_length_melted_dataframe',
                'barcode_selection_sequence_phred_dataframe', 'barcode_selection_sequence_phred_melted_dataframe',
                "all.read.barcoded",
                "read.pass.barcoded",
                "read.fail.barcoded"]

        key_list = ["basecaller.sequencing.summary.1d.extractor.sequence.length", "basecaller.sequencing.summary.1d.extractor.passes.filtering",
                    "basecaller.sequencing.summary.1d.extractor.read.pass.length", "basecaller.sequencing.summary.1d.extractor.read.fail.length",
                    "basecaller.sequencing.summary.1d.extractor.start.time.sorted",
                    "basecaller.sequencing.summary.1d.extractor.read.pass.sorted", "basecaller.sequencing.summary.1d.extractor.read.fail.sorted",
                    "basecaller.sequencing.summary.1d.extractor.mean.qscore", "basecaller.sequencing.summary.1d.extractor.read.pass.qscore",
                    "basecaller.sequencing.summary.1d.extractor.read.fail.qscore",
                    "basecaller.sequencing.summary.1d.extractor.channel.occupancy.statistics",
                    "basecaller.sequencing.summary.1d.extractor.all.read.qscore", "basecaller.sequencing.summary.1d.extractor.all.read.length"
                    ]

        for key in keys:
            key_list.append(self.add_key_to_result_dict(key))
        result_dict['unwritten.keys'].extend(key_list)

    def _occupancy_channel(self):
        """
        Statistics about the channels
        :return: channel_count_statistics containing statistics description about the channel occupancy
        """
        channel_count = self.channel
        total_number_reads_per_channel = pd.value_counts(channel_count)
        channel_count_statistics = pd.DataFrame.describe(total_number_reads_per_channel)
        return channel_count_statistics

    def _load_sequencing_summary_data(self):
        """
        Load sequencing summary data frame.
        :return: a Pandas DataFrame object
        """
        files = self.sequencing_1dsqr_summary_files

        if len(files) == 1:
            return pd.read_csv(files[0], sep="\t")

        summary_df = None
        barcode_df = None

        for f in files:
            if self._is_barcode_file(f):
                df = pd.read_csv(f, sep="\t")

                if barcode_df is None:
                    barcode_df = df
                else:
                    barcode_df = barcode_df.append(df, ignore_index=True)

            else:
                df = pd.read_csv(f, sep="\t")

                if summary_df is None:
                    summary_df = df
                else:
                    summary_df = summary_df.append(df, ignore_index=True)

        if summary_df is None:
            sys.exit("Only barcode sequencing summary found")

        if barcode_df is None:
            return summary_df

        result = summary_df.merge(barcode_df, left_on="read_id1", right_on="read_id", how="left")

        return result

    def _is_barcode_file(self, filename):
        """
        Chech if a sequencing summary file is a barcode summary file.
        :param filename: path of the file to test
        :return: True if the sequencing summary file is a barcode summary file
        """

        with open(filename) as f:
            line = f.readline()
            if not line.startswith('filename'):
                return True
            return False
