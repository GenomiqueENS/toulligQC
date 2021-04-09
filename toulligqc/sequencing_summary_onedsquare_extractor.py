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

import sys

import numpy as np
import pandas as pd

from toulligqc import plotly_graph_generator as pgg
from toulligqc import plotly_graph_onedsquare_generator as pgg2
from toulligqc.sequencing_summary_common import check_result_values
from toulligqc.sequencing_summary_common import count_boolean_elements
from toulligqc.sequencing_summary_common import describe_dict
from toulligqc.sequencing_summary_common import extract_barcode_info
from toulligqc.sequencing_summary_common import get_result_value
from toulligqc.sequencing_summary_common import series_cols_boolean_elements
from toulligqc.sequencing_summary_common import set_result_value
from toulligqc.sequencing_summary_extractor import SequencingSummaryExtractor as SSE


class OneDSquareSequencingSummaryExtractor(SSE):
    """
    Extraction of statistics from 1dsqr_sequencing_summary.txt file and graph generation
    """

    def __init__(self, config_dictionary):
        """
        Constructor that initialize the values of the config_dictionary and check in the case of 1 argument in
        sequencing_summary_source and seqencing_summary_1dsqr_source if the path points to a file,
        the others cases are managed in check_conf and _load_sequencing_summary_1dsqr_data
        :param config_dictionary: dictionary containing all files or directories paths for sequencing_summary, sequencing_1dsq_summary.txt and barcoding files
        """
        super().__init__(config_dictionary)
        self.sse = SSE(config_dictionary)
        self.sequencing_summary_1dsqr_source = self.config_dictionary[
            'sequencing_summary_1dsqr_source']
        self.sequencing_summary_1dsqr_files = self.sequencing_summary_1dsqr_source.split(
            '\t')

        # overiding attribute .is_barcode
        self.is_barcode = False
        if config_dictionary['barcoding'] == 'True':
            for f in self.sequencing_summary_1dsqr_files:
                if self._is_barcode_file(f) or self._is_sequencing_summary_with_barcodes(
                        f) or self._is_sequencing_summary_1dsqr_with_barcodes(
                    f):
                    self.is_barcode = True

    def check_conf(self):
        """
        Check if the sequencing summary 1dsqr source contains a 1dsqr sequencing summary file
        :return: boolean and a string for error message
        """
        # Check the presence of sequencing_summary.txt
        if self.sse.check_conf()[0] == False:
            return False, "No sequencing summary file has been found"

        # Check the presence of sequencing_1dsq_summary.txt
        found = False
        while not found:
            for f in self.sequencing_summary_1dsqr_files:
                try:
                    if self._is_sequencing_summary_1dsqr_file(
                            f
                    ) or self._is_sequencing_summary_1dsqr_with_barcodes(f):
                        found = True
                except FileNotFoundError:
                    return False, "No such file or directory " + f
            break
        if not found:
            return False, "No 1D squared sequencing summary file has been found"

        return True, ""

    def init(self):
        """
        Initialisation
        :return:
        """
        self.sse.init()

        # Load dataframe_1d and remove duplicate columns that are also present in dataframe_1dsqr
        self.dataframe_1d = self.sse.dataframe_1d

        # Copy dataframe to avoid changing original df when dropping columns
        dataframe_1d_copy = self.dataframe_1d.copy(deep=True)
        dataframe_1d_copy.drop(columns=["sequence_length", "mean_qscore", "passes_filtering"], inplace=True)

        # Load dataframe_1dsqr df from 1D² files
        self.dataframe_1dsqr = self._load_sequencing_summary_1dsqr_data()

        # Create duration column in dataframe_1dsqr
        self.dataframe_1dsqr['duration'] = self.dataframe_1dsqr['trimmed_duration1'] + self.dataframe_1dsqr[
            'trimmed_duration2']  # duration of the 2 strands sequenced
        self.dataframe_1dsqr.drop(columns=['trimmed_duration1', 'trimmed_duration2'], inplace=True)

        # dataframe_dicts
        self.dataframe_dict_1dsqr = {}
        self.dataframe_dict = self.sse.dataframe_dict

        if self.is_barcode:
            self.barcode_selection = self.config_dictionary[
                'barcode_selection']

        if self.dataframe_1d.empty or self.dataframe_1dsqr.empty:
            raise pd.errors.EmptyDataError("Dataframe is empty")

        # Dictionary for storing all pd.Series and pd.Dataframe entries
        self.dataframe_dict = {}

        if self.is_barcode:
            self.barcode_selection = self.config_dictionary[
                'barcode_selection']

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

    def extract(self, result_dict):
        """
        :param result_dict:
        :return:
        """
        #
        # Extract from 1D summary source
        #

        # If key already in result_dict, replace 1D by 1D² Basecaller analysis
        if 'sequencing.telemetry.extractor.software.analysis' in result_dict:
            result_dict['sequencing.telemetry.extractor.software.analysis'] = '1dsqr_basecalling'

        # Call to extract parent method to get all keys, values from 1D extractor
        self.sse.extract(result_dict)
        self.dataframe_dict = self.sse.dataframe_dict

        #
        # Extract info from 1D²
        #

        self._fill_series_dict(self.dataframe_dict_1dsqr, self.dataframe_1dsqr)

        # Read count
        set_result_value(self, result_dict, "read.count", len(self.dataframe_1dsqr))

        # 1D² pass information : count, length and qscore values
        set_result_value(self, result_dict, "read.pass.count",
                         count_boolean_elements(self.dataframe_1dsqr, 'passes_filtering', True))

        # 1D² fail information : count, length and qscore values
        set_result_value(self, result_dict, "read.fail.count",
                         count_boolean_elements(self.dataframe_1dsqr, 'passes_filtering', False))

        # Ratios & frequencies
        set_result_value(self, result_dict, "read.count.frequency", 100)
        total_reads = get_result_value(self, result_dict, "read.count")
        set_result_value(self, result_dict, "read.pass.ratio",
                         (get_result_value(self, result_dict, "read.pass.count") / total_reads))
        set_result_value(self, result_dict, "read.fail.ratio",
                         (get_result_value(self, result_dict, "read.fail.count") / total_reads))

        read_pass_frequency = (get_result_value(self, result_dict, "read.pass.count") / total_reads) * 100
        set_result_value(self, result_dict, "read.pass.frequency", read_pass_frequency)

        read_fail_frequency = (get_result_value(self, result_dict, "read.fail.count") / total_reads) * 100
        set_result_value(self, result_dict, "read.fail.frequency", read_fail_frequency)

        # Get statistics about all reads length and store each value into result_dict
        sequence_length_statistics = self.dataframe_1dsqr['sequence_length'].describe()

        for index, value in sequence_length_statistics.items():
            set_result_value(self,
                             result_dict, "all.read.length." + index, value)

        # Add statistics (without count) about read pass/fail length in the result_dict
        describe_dict(self, result_dict, self.dataframe_dict_1dsqr["pass.reads.sequence.length"],
                      "pass.reads.sequence.length")
        describe_dict(self, result_dict, self.dataframe_dict_1dsqr["fail.reads.sequence.length"],
                      "fail.reads.sequence.length")

        # Get Qscore statistics without count value and store them into result_dict
        qscore_statistics = self.dataframe_1dsqr['mean_qscore'].describe().drop(
            "count")

        for index, value in qscore_statistics.items():
            set_result_value(self,
                             result_dict, "all.reads.mean.qscore." + index, value)

        # Add statistics (without count) about read pass/fail qscore in the result_dict
        describe_dict(self, result_dict, self.dataframe_dict_1dsqr["pass.reads.mean.qscore"], "pass.reads.mean.qscore")
        describe_dict(self, result_dict, self.dataframe_dict_1dsqr["fail.reads.mean.qscore"], "fail.reads.mean.qscore")

        if self.is_barcode:
            extract_barcode_info(self,
                                 result_dict,
                                 self.barcode_selection,
                                 self.dataframe_dict_1dsqr,
                                 self.dataframe_1dsqr)

    def _fill_series_dict(self, df_dict, df):

        for read_type in ['pass', 'fail']:
            read_type_bool = True if read_type == 'pass' else False

            self.dataframe_dict_1dsqr[read_type + '.reads.sequence.length'] = \
                series_cols_boolean_elements(self.dataframe_1dsqr,
                                             'sequence_length',
                                             'passes_filtering',
                                             read_type_bool)

            self.dataframe_dict_1dsqr[read_type + '.reads.mean.qscore'] = \
                series_cols_boolean_elements(self.dataframe_1dsqr,
                                             'mean_qscore',
                                             'passes_filtering',
                                             read_type_bool)

        # Read length & passes_filtering & qscore information
        self.dataframe_dict_1dsqr["all.reads.sequence.length"] = self.dataframe_1dsqr["sequence_length"]
        self.dataframe_dict_1dsqr["passes.filtering"] = self.dataframe_1dsqr["passes_filtering"]
        self.dataframe_dict_1dsqr["all.reads.mean.qscore"] = self.dataframe_1dsqr["mean_qscore"]

        self.dataframe_dict_1dsqr["all.reads.start.time1"] = self.dataframe_1dsqr['start_time1']
        self.dataframe_dict_1dsqr["all.reads.duration"] = self.dataframe_1dsqr['duration']

    def graph_generation(self, result_dict):
        """
        Generation of the differents graphs containing in the plotly_graph_generator modules
        :return: images array containing the title and the path toward the images
        """

        images = list([pgg.read_count_histogram(result_dict, self.images_directory)])
        images.append(pgg2.dsqr_read_count_histogram(result_dict, self.images_directory))
        images.append(pgg.read_length_scatterplot(self.dataframe_dict, self.images_directory))
        images.append(pgg2.dsqr_read_length_scatterplot(self.dataframe_dict_1dsqr, self.images_directory))
        images.append(pgg.yield_plot(self.dataframe_1dsqr, self.images_directory, oneDsquare=True))
        images.append(pgg.read_quality_multiboxplot(self.dataframe_dict, self.images_directory, ))
        images.append(
            pgg2.dsqr_read_quality_multiboxplot(result_dict, self.dataframe_dict_1dsqr, self.images_directory))
        images.append(pgg.allphred_score_frequency(self.dataframe_dict, self.images_directory))
        images.append(pgg2.dsqr_allphred_score_frequency(result_dict, self.dataframe_dict_1dsqr,
                                                         self.images_directory))
        images.append(pgg.all_scatterplot(self.dataframe_dict, self.images_directory))
        images.append(pgg2.scatterplot_1dsqr(self.dataframe_dict_1dsqr, self.images_directory))
        images.append(pgg.plot_performance(self.dataframe_dict, self.images_directory))
        images.append(
            pgg2.sequence_length_over_time_dsqr(self.dataframe_dict_1dsqr, self.images_directory))
        images.append(pgg2.phred_score_over_time_dsqr(result_dict, self.dataframe_dict_1dsqr, self.images_directory))
        images.append(pgg2.speed_over_time_dsqr(self.dataframe_dict_1dsqr, self.images_directory))

        if self.is_barcode:
            images.append(pgg2.barcode_percentage_pie_chart_1dsqr_pass(self.dataframe_dict_1dsqr,
                                                                       self.barcode_selection,
                                                                       self.images_directory))

            images.append(pgg2.barcode_percentage_pie_chart_1dsqr_fail(self.dataframe_dict_1dsqr,
                                                                       self.barcode_selection,
                                                                       self.images_directory))

            images.append(pgg2.barcode_length_boxplot_1dsqr(self.dataframe_dict_1dsqr,
                                                            self.images_directory))

            images.append(pgg2.barcoded_phred_score_frequency_1dsqr(self.dataframe_dict_1dsqr,
                                                                    self.images_directory))
        return images

    def clean(self, result_dict):
        """
        Removing dictionary entries that will not be kept in the report.data file
        :param result_dict:
        :return:
        """

        # Clean SSE
        self.sse.clean(result_dict)

        # Check values in result_dict (avoid Series and Dataframes)
        check_result_values(self, result_dict)

        # Clear dictionary for Series and Dataframe
        self.dataframe_dict_1dsqr.clear()

        # Clear DataFrame
        self.dataframe_1dsqr.iloc[0:0]

    def _occupancy_channel(self):
        """
        Statistics about the channels
        :return: channel_count_statistics containing statistics description about the channel occupancy
        """
        channel_count = self.channel_df
        total_number_reads_per_channel = pd.value_counts(channel_count)
        channel_count_statistics = pd.DataFrame.describe(total_number_reads_per_channel)
        return channel_count_statistics

    def _load_sequencing_summary_1dsqr_data(self):
        """
        Load sequencing summary dataframe with or without barcodes
        :return: a Pandas Dataframe object
        """
        # Initialization
        files = self.sequencing_summary_1dsqr_files

        summary_dataframe = None
        barcode_dataframe = None

        sequencing_summary_columns = [
            'passes_filtering',
            'sequence_length', 'mean_qscore',
            'start_time1',
            'trimmed_duration1', 'trimmed_duration2'
        ]

        sequencing_summary_datatypes = {
            'passes_filtering': np.bool,
            'sequence_length': np.uint32,
            'mean_qscore': np.float,
            'start_time1': np.float,
            'trimmed_duration1': np.float,
            'trimmed_duration2': np.float,
        }

        # If barcoding files are provided, merging of dataframes must be done on read_id column
        barcoding_summary_columns = ['read_id', 'barcode_arrangement']

        barcoding_summary_datatypes = {
            'read_id': object,
            'barcode_arrangement': object
        }

        try:
            # If 1 file and it's a 1dsqr_sequencing_summary.txt
            if len(files) == 1 and self._is_sequencing_summary_1dsqr_file(
                    files[0]):
                return pd.read_csv(files[0],
                                   sep="\t",
                                   usecols=sequencing_summary_columns,
                                   dtype=sequencing_summary_datatypes)

            # If 1 file and it's a 1_dsqr_sequencing_summary.txt with barcode info, load column barcode_arrangement
            elif len(
                    files
            ) == 1 and self._is_sequencing_summary_1dsqr_with_barcodes(
                files[0]):
                sequencing_summary_columns.append('barcode_arrangement')
                sequencing_summary_datatypes.update(
                    {'barcode_arrangement': object})

                return pd.read_csv(files[0],
                                   sep="\t",
                                   usecols=sequencing_summary_columns,
                                   dtype=sequencing_summary_datatypes)

            # If multiple files, check if there's a barcoding one and a sequencing one :
            for f in files:

                # check for presence of barcoding files
                if self._is_barcode_file(f):
                    dataframe = pd.read_csv(f,
                                            sep="\t",
                                            usecols=barcoding_summary_columns,
                                            dtype=barcoding_summary_datatypes)
                    if barcode_dataframe is None:
                        barcode_dataframe = dataframe
                    # if a barcoding file has already been read, append the 2 dataframes
                    else:
                        barcode_dataframe = barcode_dataframe.append(
                            dataframe, ignore_index=True)

                # check for presence of sequencing_summary file, if True add column read_id for merging with barcode dataframe
                else:
                    if self._is_sequencing_summary_1dsqr_file(f):
                        sequencing_summary_columns.append('read_id1')
                        sequencing_summary_datatypes.update(
                            {'read_id1': object})

                        dataframe = pd.read_csv(
                            f,
                            sep="\t",
                            usecols=sequencing_summary_columns,
                            dtype=sequencing_summary_datatypes)
                        if summary_dataframe is None:
                            summary_dataframe = dataframe
                        else:
                            summary_dataframe = summary_dataframe.append(
                                dataframe, ignore_index=True)

            if barcode_dataframe is None:
                # If no barcodes in files, no merged dataframes on column 'read_id'
                return summary_dataframe.drop(columns=['read_id1'])
            else:
                summary_dataframe.rename(columns={"read_id1": "read_id"},
                                         inplace=True)
                dataframes_merged = pd.merge(summary_dataframe,
                                             barcode_dataframe,
                                             on='read_id',
                                             how='left')

                missing_barcodes_count = dataframes_merged['barcode_arrangement'].isna().sum()
                if missing_barcodes_count > 0:
                    sys.stderr.write('Warning: {} barcodes values are missing in sequencing summary file(s).'
                                     ' They will be marked as "unclassified".\n'.format(missing_barcodes_count))

                # Replace missing barcodes values by 'unclassified'
                dataframes_merged['barcode_arrangement'] = dataframes_merged['barcode_arrangement'].fillna(
                    'unclassified')

                # delete column read_id after merging
                del dataframes_merged['read_id']

                return dataframes_merged

        except IOError:
            raise FileNotFoundError("Sequencing summary file not found")

    @staticmethod
    def _is_barcode_file(filename):
        """
        Check if input is a barcoding summary file i.e. has the column barcode_arrangement
        :param filename: path of the file to test
        :return: True if the filename is a barcoding summary file
        """
        try:
            with open(filename, 'r') as f:
                header = f.readline()
            return header.startswith(
                'read_id') and 'barcode_arrangement' in header
        except FileNotFoundError:
            "No barcoding file was found"

    # Static methods for checking 1dsqr files

    @staticmethod
    def _is_sequencing_summary_1dsqr_file(filename):
        """
        Check if input is a sequencing summary file i.e. first word is "filename" and does not have column 'barcode_arrangement'
        :param filename: path of the file to test
        :return: True if the file is indeed a sequencing summary file
        """
        try:
            with open(filename, 'r') as f:
                header = f.readline()
            return header.startswith(
                'filename1') and not 'barcode_arrangement' in header
        except IOError:
            raise FileNotFoundError

    @staticmethod
    def _is_sequencing_summary_1dsqr_with_barcodes(filename):
        pass
