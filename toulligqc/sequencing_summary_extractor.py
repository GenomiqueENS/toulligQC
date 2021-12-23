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
# Maintainer: Karine Dias
# Since version 0.1

# Extraction of statistics from sequencing_summary.txt file (1D chemistry)

import sys
import time

import numpy as np
import pandas as pd

from toulligqc import plotly_graph_generator as pgg
from toulligqc.sequencing_summary_common import check_result_values
from toulligqc.sequencing_summary_common import count_boolean_elements
from toulligqc.sequencing_summary_common import describe_dict
from toulligqc.sequencing_summary_common import extract_barcode_info
from toulligqc.sequencing_summary_common import get_result_value
from toulligqc.sequencing_summary_common import series_cols_boolean_elements
from toulligqc.sequencing_summary_common import set_result_value
from toulligqc.sequencing_summary_common import log_task
from toulligqc.sequencing_summary_common import add_image_to_result
from toulligqc.sequencing_summary_common import read_first_line_file


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
        self.sequencing_summary_source = config_dictionary['sequencing_summary_source']
        self.images_directory = config_dictionary['images_directory']
        self.sequencing_summary_files = self.sequencing_summary_source.split('\t')
        if 'quiet' not in config_dictionary or config_dictionary['quiet'].lower() != 'true':
            self.quiet = False
        else:
            self.quiet = True

        self.is_barcode = False
        if config_dictionary['barcoding'] == 'True':
            for f in self.sequencing_summary_files:
                if self._is_barcode_file(f) or self._is_sequencing_summary_with_barcodes(f):
                    self.is_barcode = True

    def check_conf(self):
        """
        Check if the sequencing summary source contains a sequencing summary file
        :return: boolean and a string for error message
        """

        if not self.sequencing_summary_files[0]:
            return False, "No file has been defined"

        found = False
        while not found:
            for f in self.sequencing_summary_files:
                try:
                    if self._is_sequencing_summary_file(f) or self._is_sequencing_summary_with_barcodes(f):
                        found = True
                except FileNotFoundError:
                    return False, "No such file or directory " + f
            break

        if not found:
            return False, "No sequencing summary file has been found"
        return True, ""

    def init(self):
        """
        Creation of the dataframe containing all info from sequencing_summary.txt
        :return: Panda's Dataframe object
        """

        start_time = time.time()

        self.dataframe_1d = self._load_sequencing_summary_data()
        if self.dataframe_1d.empty:
            raise pd.errors.EmptyDataError("Dataframe is empty")

        # Rename 'sequence_length_template' and 'mean_qscore_template'
        self.dataframe_1d.rename(columns={'sequence_length_template': 'sequence_length',
                                          'mean_qscore_template': 'mean_qscore'}, inplace=True)

        # Add missing categories
        if 'barcode_arrangement' in self.dataframe_1d.columns:
            self.dataframe_1d['barcode_arrangement'].cat.add_categories([0, 'other barcodes', 'passes_filtering'],
                                                                        inplace=True)

        # Replace all NaN values by 0 to avoid data manipulation errors when columns are not the same length
        self.dataframe_1d = self.dataframe_1d.fillna(0)

        # Dictionary for storing all pd.Series and pd.Dataframe entries
        self.dataframe_dict = {}

        if self.is_barcode:
            self.barcode_selection = self.config_dictionary['barcode_selection']

        log_task(self.quiet,
                 'Load sequencing summary file ({:,.2f} MB used)'.format(self.dataframe_1d.memory_usage(deep=True).sum()/1024/1024),
                 start_time,
                 time.time())

    @staticmethod
    def get_name() -> str:
        """
        Get the name of the extractor.
        :return: the name of the extractor
        """
        return 'Basecaller sequencing summary'

    @staticmethod
    def get_report_data_file_id() -> str:
        """
        Get the report.data id of the extractor.
        :return: a string with the name of the ID of the extractor
        """
        return 'basecaller.sequencing.summary.1d.extractor'

    def extract(self, result_dict):
        """
        Get Phred score (Qscore) and Length details (frequencies, ratios, yield and statistics) per type read (pass or fail)
        :param result_dict:
        """

        start_time = time.time()

        # Basecaller analysis
        if 'sequencing.telemetry.extractor.software.analysis' not in result_dict:
            result_dict['sequencing.telemetry.extractor.software.analysis'] = '1d_basecalling'

        self._fill_series_dict(self.dataframe_dict, self.dataframe_1d)

        # Read count
        set_result_value(self, result_dict, "read.count", len(self.dataframe_1d))

        # 1D pass information : count, length, qscore values and sorted Series
        set_result_value(self, result_dict, "read.pass.count",
                         count_boolean_elements(self.dataframe_1d, 'passes_filtering', True))

        # 1D fail information : count, length, qscore values and sorted Series
        set_result_value(self, result_dict, "read.fail.count",
                         count_boolean_elements(self.dataframe_1d, 'passes_filtering', False))

        total_reads = get_result_value(self, result_dict, "read.count")

        # Ratios
        set_result_value(self, result_dict, "read.pass.ratio",
                         (get_result_value(self, result_dict, "read.pass.count") / total_reads))
        set_result_value(self, result_dict, "read.fail.ratio",
                         (get_result_value(self, result_dict, "read.fail.count") / total_reads))

        # Frequencies
        set_result_value(self, result_dict, "read.count.frequency", 100)

        read_pass_frequency = (get_result_value(self,
                                                result_dict, "read.pass.count") / total_reads) * 100
        set_result_value(self,
                         result_dict, "read.pass.frequency", read_pass_frequency)

        read_fail_frequency = (get_result_value(self,
                                                result_dict, "read.fail.count") / total_reads) * 100
        set_result_value(self,
                         result_dict, "read.fail.frequency", read_fail_frequency)

        # Yield, n50, run time
        set_result_value(self, result_dict, "yield", sum(self.dataframe_dict["all.reads.sequence.length"]))

        set_result_value(self, result_dict, "n50", self._compute_NXX(50))
        set_result_value(self, result_dict, "l50", self._compute_LXX(50))

        set_result_value(self, result_dict, "run.time", max(self.dataframe_1d['start_time']))

        # Get channel occupancy statistics and store each value into result_dict
        for index, value in self._occupancy_channel().items():
            set_result_value(self,
                             result_dict, "channel.occupancy.statistics." + index, value)

        # Get statistics about all reads length and store each value into result_dict
        sequence_length_statistics = self.dataframe_dict["all.reads.sequence.length"].describe()

        for index, value in sequence_length_statistics.items():
            set_result_value(self,
                             result_dict, "all.read.length." + index, value)

        # Add statistics (without count) about read pass/fail length in the result_dict
        describe_dict(self, result_dict, self.dataframe_dict["pass.reads.sequence.length"],
                      "pass.reads.sequence.length")
        describe_dict(self, result_dict, self.dataframe_dict["fail.reads.sequence.length"],
                      "fail.reads.sequence.length")

        # Get Qscore statistics without count value and store them into result_dict
        qscore_statistics = self.dataframe_1d['mean_qscore'].describe().drop(
            "count")

        for index, value in qscore_statistics.items():
            set_result_value(self,
                             result_dict, "all.read.qscore." + index, value)

        # Add statistics (without count) about read pass/fail qscore in the result_dict
        describe_dict(self, result_dict, self.dataframe_dict["pass.reads.mean.qscore"], "pass.reads.mean.qscore")
        describe_dict(self, result_dict, self.dataframe_dict["fail.reads.mean.qscore"], "fail.reads.mean.qscore")

        if self.is_barcode:
            extract_barcode_info(self, result_dict,
                                 self.barcode_selection,
                                 self.dataframe_dict,
                                 self.dataframe_1d)

        log_task(self.quiet, 'Extract info from sequencing summary file', start_time, time.time())

    def _fill_series_dict(self, df_dict, df):

        for read_type in ['pass', 'fail']:
            read_type_bool = True if read_type == 'pass' else False

            # Read length series
            df_dict[read_type + '.reads.sequence.length'] = series_cols_boolean_elements(df,
                                                                                         'sequence_length',
                                                                                         'passes_filtering',
                                                                                         read_type_bool)

            # Read qscore series
            df_dict[read_type + '.reads.mean.qscore'] = series_cols_boolean_elements(df,
                                                                                     'mean_qscore',
                                                                                     'passes_filtering',
                                                                                     read_type_bool)

        # Read length series
        df_dict["all.reads.sequence.length"] = df['sequence_length']

        # Mean QScore
        df_dict["all.reads.mean.qscore"] = df['mean_qscore']

        # Time series
        df_dict["all.reads.start.time"] = df['start_time']

        # Duration series
        df_dict["all.reads.duration"] = df['duration']

    def graph_generation(self, result_dict):
        """
        Generation of the different graphs containing in the plotly_graph_generator module
        :return: images array containing the title and the path toward the images
        """
        images = list()

        add_image_to_result(self.quiet, images, time.time(), pgg.read_count_histogram(result_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.read_length_scatterplot(self.dataframe_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.yield_plot(self.dataframe_1d, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.read_quality_multiboxplot(self.dataframe_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.allphred_score_frequency(self.dataframe_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.plot_performance(self.dataframe_1d, self.images_directory))

        add_image_to_result(self.quiet, images, time.time(), pgg.all_scatterplot(self.dataframe_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.sequence_length_over_time(self.dataframe_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.phred_score_over_time(self.dataframe_dict, result_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.speed_over_time(self.dataframe_dict, self.images_directory))

        if self.is_barcode:
            add_image_to_result(self.quiet, images, time.time(), pgg.barcode_percentage_pie_chart_pass(self.dataframe_dict,
                                                                                                       self.barcode_selection,
                                                                                                       self.images_directory))

            add_image_to_result(self.quiet, images, time.time(), pgg.barcode_percentage_pie_chart_fail(self.dataframe_dict,
                                                                                                       self.barcode_selection,
                                                                                                       self.images_directory))

            add_image_to_result(self.quiet, images, time.time(), pgg.barcode_length_boxplot(self.dataframe_dict,
                                                                                            self.images_directory))

            add_image_to_result(self.quiet, images, time.time(), pgg.barcoded_phred_score_frequency(self.dataframe_dict,
                                                                                                    self.images_directory))
        return images

    def clean(self, result_dict):
        """
        Removing dictionary entries that will not be kept in the report.data file
        :return:
        """

        # Check values in result_dict (avoid Series and Dataframes)
        check_result_values(self, result_dict)

        # Clear dictionary for Series and Dataframe
        self.dataframe_dict.clear()

        # Clear DataFrame
        self.dataframe_1d.iloc[0:0]

    def _occupancy_channel(self):
        """
        Statistics about the channels of the flowcell
        :return: pd.Series object containing statistics about the channel occupancy without count value
        """
        total_reads_per_channel = pd.value_counts(self.dataframe_1d["channel"])
        return pd.DataFrame.describe(total_reads_per_channel)

    def _load_sequencing_summary_data(self):
        """
        Load sequencing summary dataframe with or without barcodes
        :return: a Pandas Dataframe object
        """
        # Initialization
        files = self.sequencing_summary_files

        summary_dataframe = None
        barcode_dataframe = None

        sequencing_summary_columns = ['channel', 'start_time',
                                      'passes_filtering',
                                      'sequence_length_template',
                                      'mean_qscore_template',
                                      'duration']

        sequencing_summary_datatypes = {
            'channel': np.int16,
            'start_time': np.float64,
            'passes_filtering': np.bool,
            'sequence_length_template': np.uint32,
            'mean_qscore_template': np.float32,
            'duration': np.float32}

        # If barcoding files are provided, merging of dataframes must be done on read_id column
        barcoding_summary_columns = ['read_id', 'barcode_arrangement']

        barcoding_summary_datatypes = {
            'read_id': object,
            'barcode_arrangement': 'category'
        }

        try:
            # If 1 file and it's a sequencing_summary.txt
            if len(files) == 1 and self._is_sequencing_summary_file(files[0]):
                return pd.read_csv(files[0], sep="\t", usecols=sequencing_summary_columns,
                                   dtype=sequencing_summary_datatypes)

            # If 1 file and it's a sequencing_summary.txt with barcode info, load column barcode_arrangement
            elif len(files) == 1 and self._is_sequencing_summary_with_barcodes(files[0]):
                sequencing_summary_columns.append('barcode_arrangement')
                sequencing_summary_datatypes.update(
                    {'barcode_arrangement': 'category'})

                return pd.read_csv(files[0], sep="\t", usecols=sequencing_summary_columns,
                                   dtype=sequencing_summary_datatypes)

            # If multiple files, check if there's a barcoding one and a sequencing one :
            for f in files:

                # check for presence of barcoding files
                if self._is_barcode_file(f):
                    dataframe = pd.read_csv(
                        f, sep="\t", usecols=barcoding_summary_columns, dtype=barcoding_summary_datatypes)
                    if barcode_dataframe is None:
                        barcode_dataframe = dataframe
                    # if a barcoding file has already been read, append the 2 dataframes
                    else:
                        barcode_dataframe = barcode_dataframe.append(
                            dataframe, ignore_index=True)

                # check for presence of sequencing_summary file, if True add column read_id for merging with barcode dataframe
                else:
                    if self._is_sequencing_summary_file(f):
                        sequencing_summary_columns.append('read_id')
                        sequencing_summary_datatypes.update(
                            {'read_id': object})

                        dataframe = pd.read_csv(f, sep="\t", usecols=sequencing_summary_columns,
                                                dtype=sequencing_summary_datatypes)
                        if summary_dataframe is None:
                            summary_dataframe = dataframe
                        else:
                            summary_dataframe = summary_dataframe.append(
                                dataframe, ignore_index=True)

            if barcode_dataframe is None:
                # If no barcodes in files, no merged dataframes on column 'read_id'
                return summary_dataframe.drop(columns=['read_id'])
            else:
                dataframes_merged = pd.merge(
                    summary_dataframe, barcode_dataframe, on='read_id', how='left')

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

    def _compute_NXX(self, x):
        """Compute NXX value of total sequence length"""
        data = self.dataframe_dict["all.reads.sequence.length"].dropna().values
        data.sort()
        half_sum = data.sum() * x / 100
        cum_sum = 0
        for v in data:
            cum_sum += v
            if cum_sum >= half_sum:
                return int(v)

    def _compute_LXX(self, x):
        """Compute LXX value of total sequence length"""
        data = self.dataframe_dict["all.reads.sequence.length"].dropna().values
        data.sort()
        half_sum = data.sum() * x / 100
        cum_sum = 0
        count = 0
        for v in data:
            cum_sum += v
            count += 1
            if cum_sum >= half_sum:
                return count

    @staticmethod
    def _is_barcode_file(filename):
        """
        Check if input is a barcoding summary file i.e. has the column barcode_arrangement
        :param filename: path of the file to test
        :return: True if the filename is a barcoding summary file
        """
        header = read_first_line_file(filename)
        return header.startswith('read_id') and 'barcode_arrangement' in header

    @staticmethod
    def _is_sequencing_summary_file(filename):
        """
        Check if input is a sequencing summary file i.e. first word is "filename" and does not have column 'barcode_arrangement'
        :param filename: path of the file to test
        :return: True if the file is indeed a sequencing summary file
        """
        header = read_first_line_file(filename)
        return header.startswith('filename') and not 'barcode_arrangement' in header

    @staticmethod
    def _is_sequencing_summary_with_barcodes(filename):
        """
        Check if the sequencing summary has also barcode information :
        - check for presence of columns "filename" and "barcode_arrangement"
        :param filename: path of the file to test
        :return: True if the filename is a sequencing summary file with barcodes
        """
        header = read_first_line_file(filename)
        return header.startswith('filename') and 'barcode_arrangement' in header



