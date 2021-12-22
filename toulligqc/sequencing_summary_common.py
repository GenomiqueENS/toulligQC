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
# First author: Laurent Jourdren
# Maintainer: Laurent Jourdren
# Since version 2.0

# This module contains common methods for sequencing summary modules.

import sys
import gzip
import bz2
import time
import pandas as pd
from toulligqc import common


def set_result_value(extractor, result_dict, key: str, value):
    """
    Set a key, value pair to the result_dict
    :param result_dict:
    :param key: string entry to add to result_dict
    :param value: int, float, list, pd.Series or pd.Dataframe value of the corresponding key
    """

    _check_result_key_value(key, value)
    result_dict[extractor.get_report_data_file_id() + '.' + key] = value


def get_result_value(extractor, result_dict, key: str):
    """
    :param result_dict:
    :param key: string entry to add to result_dict
    Returns the value associated with the result_dict key
    """
    if not (extractor.get_report_data_file_id() + '.' + key) in result_dict.keys():
        raise KeyError("Key {key} not found").__format__(key)
    return result_dict.get(extractor.get_report_data_file_id() + '.' + key)


def check_result_values(extractor, result_dict):
    prefix = extractor.get_report_data_file_id() + '.'
    for key, value in result_dict.items():
        if key.startswith(prefix):
            _check_result_key_value(key, value)


def _check_result_key_value(key, value):
    if not isinstance(key, str):
        raise TypeError("Invalid type for key: {}".format(type(key)))

    if not isinstance(value, int) and not isinstance(value, float) and not isinstance(value, str):
        raise TypeError("Invalid type for the value of the key {}: {} ".format(key, type(value)))


def describe_dict(extractor, result_dict: dict, function, entry: str):
    """
    Set statistics for a key like mean, min, max, median and percentiles (without the count value) filled in the _set_result_value dictionary
    :param result_dict:
    :param function: function returning the values to describe
    :param entry: entry to put in result_dict completed with the statistics
    """
    stats = pd.Series.describe(function).drop("count")
    for key, value in stats.iteritems():
        set_result_value(extractor, result_dict, entry + '.' + key, value)


def count_boolean_elements(dataframe, column_name, boolean: bool) -> int:
    """
    Returns the number of values of a column filtered by a boolean
    :param colum_name: name of the dafatrame column
    :boolean: bool to filter
    """
    return len(dataframe.loc[dataframe[column_name] == bool(boolean)])


def series_cols_boolean_elements(dataframe, column_name1: str, column_name2: str, boolean: bool) -> pd.Series:
    """
    Returns a Panda's Series object with the number of values of different columns filtered by a boolean
    :param dataframe: dataframe_1d
    :param column_name1: 1st column to filter
    :param column_name2: 2nd column to filter
    :param boolean: access columns of dataframe by boolean array
    """
    return dataframe[column_name1].loc[dataframe[column_name2] == bool(boolean)]


def sorted_series_boolean_elements_divided(dataframe, column_name1: str, column_name2: str, boolean: bool,
                                           denominator: int):
    """
    Returns a sorted series of values of different columns filtered by a boolean and divided by the denominator
    :param dataframe: dataframe_1d
    :param column_name1: 1st column to filter
    :param column_name2: 2nd column to filter
    :param boolean: access columns of dataframe by boolean array
    :param denominator: number to divide by
    """
    return (dataframe[column_name1].loc[dataframe[column_name2] == bool(boolean)] / denominator).sort_values()


def extract_barcode_info(extractor, result_dict, barcode_selection, dataframe_dict, df):
    """
    :param result_dict:
    Gather all barcode info for graphs : reads pass/fail and frequency per barcodes
    """
    # Add values unclassified and other to barcode list
    if "unclassified" not in barcode_selection:
        barcode_selection.append("unclassified")

    # Create keys barcode.arrangement, and read.pass/fail.barcode in dataframe_dict with all values of
    # column barcode_arrangement when reads are passed/failed
    dataframe_dict["barcode.arrangement"] = df["barcode_arrangement"]

    # Print warning message if a barcode is unknown
    barcodes_found = set(dataframe_dict["barcode.arrangement"].unique())
    for element in barcode_selection:
        if element not in barcodes_found and element != 'other barcodes':
            sys.stderr.write("Warning: The barcode {} doesn't exist in input data\n".format(element))

    # Get barcodes frequency by read type
    series_read_pass_barcode = series_cols_boolean_elements(df, "barcode_arrangement",
                                                            "passes_filtering", True)

    dataframe_dict["read.pass.barcoded"] = _barcode_frequency(extractor, barcode_selection, result_dict,
                                                              "read.pass.barcoded",
                                                              series_read_pass_barcode)

    series_read_fail_barcode = series_cols_boolean_elements(df, "barcode_arrangement",
                                                            "passes_filtering", False)

    dataframe_dict["read.fail.barcoded"] = _barcode_frequency(extractor, barcode_selection, result_dict,
                                                              "read.fail.barcoded",
                                                              series_read_fail_barcode)

    read_pass_barcoded_count = get_result_value(extractor, result_dict, 'read.pass.barcoded.count')
    read_fail_barcoded_count = get_result_value(extractor, result_dict, 'read.fail.barcoded.count')

    # Add key "read.pass.barcoded.frequency"
    total_reads = get_result_value(extractor, result_dict, "read.count")
    set_result_value(extractor, result_dict, "read.pass.barcoded.frequency",
                     (read_pass_barcoded_count / total_reads) * 100)

    # Add key "read.fail.barcoded.frequency"
    set_result_value(extractor, result_dict, "read.fail.barcoded.frequency",
                     (read_fail_barcoded_count / total_reads) * 100)

    # Replaces all rows with unused barcodes (ie not in barcode_selection) in column barcode_arrangement with the 'other' value
    df.loc[~df['barcode_arrangement'].isin(
        barcode_selection), 'barcode_arrangement'] = 'other barcodes'

    if 'other barcodes' not in barcode_selection:
        barcode_selection.append('other barcodes')

    # Create dataframes filtered by barcodes and read quality
    for index_barcode, barcode in enumerate(barcode_selection):
        barcode_all_reads_df = df[df['barcode_arrangement'] == barcode]
        barcode_pass_reads_df = barcode_all_reads_df.loc[df['passes_filtering'] == bool(True)]
        barcode_fail_reads_df = barcode_all_reads_df.loc[df['passes_filtering'] == bool(False)]

        # Add all barcode statistics to result_dict based on values of selected dataframes
        _barcode_stats(extractor,
                       result_dict,
                       barcode_all_reads_df,
                       barcode_pass_reads_df,
                       barcode_fail_reads_df,
                       barcode)

    # Add filtered dataframes (all info by barcode and by length or qscore) to dataframe_dict
    _barcode_selection_dataframe(dataframe_dict, df, "sequence_length",
                                 "barcode_selection_sequence_length_dataframe",
                                 "length")
    _barcode_selection_dataframe(dataframe_dict, df, "mean_qscore",
                                 "barcode_selection_sequence_phred_dataframe",
                                 "qscore")


def _barcode_selection_dataframe(dataframe_dict, df, df_column_name: str, df_key_name: str,
                                 melted_column_name: str):
    """
    Create custom dataframes by grouping all reads per barcodes and per read type (pass/fail) for read length or phred score info
    Reshape the dataframes from wide to long format to display barcode, read type and read length or phred score per read
    These dataframes are used for sequence length and qscore boxplots
    :param key: string name to put in dataframe_dict
    :param df_column_name: name of the dataframe_1d column used for the new barcode_selection_dataframes
    :param melted_column_name: value (qscore or length) to use for renaming column of melted dataframe
    """
    # Count total number of rows
    nrows = df.shape[0]
    # Create a new dataframe with 3 columns : 'passes_filtering', 'barcode_arrangement' and the column name parameter
    filtered_df = df.filter(
        items=['passes_filtering', df_column_name, 'barcode_arrangement'])

    # Reshape dataframe with new MultiIndex : numbered index of df length + passes filtering index and then shape data by barcode
    barcode_selection_dataframe = filtered_df.set_index([pd.RangeIndex(start=0, stop=nrows), 'passes_filtering'],
                                                        drop=True).pivot(columns="barcode_arrangement")

    # Remove the column parameter index
    barcode_selection_dataframe.columns.droplevel(level=0)

    # Remove sequence_length Multindex to only have barcode_arrangement column labels
    barcode_selection_dataframe.columns = barcode_selection_dataframe.columns.droplevel(
        level=0)

    # Reset index to have all labels in the same level
    barcode_selection_dataframe.reset_index(
        level='passes_filtering', inplace=True)

    # Add final dataframe to dataframe_dict
    dataframe_dict[df_key_name] = barcode_selection_dataframe


def _barcode_stats(extractor, result_dict, barcode_selected_dataframe, barcode_selected_read_pass_dataframe,
                   barcode_selected_read_fail_dataframe, barcode_name):
    """
    :param result_dict:
    :param prefix: report.data id of the extractor (string)
    :param barcode_selected: barcode filtered dataframes
    Put statistics (with describe method) about barcode length and qscore in result_dict for each selected dataframe : all.read/read.pass and read.fail
    N.b. does not include count statistic for qscore
    """
    df_dict = {'all.read.': barcode_selected_dataframe,
               'read.pass.': barcode_selected_read_pass_dataframe,
               'read.fail.': barcode_selected_read_fail_dataframe}

    for df_name, df in df_dict.items():  # df_dict.items = all.read/read.pass/read.fail
        for stats_index, stats_value in df['sequence_length'].describe().items():
            key_to_result_dict = df_name + barcode_name.replace(' ', '.') + '.length.' + stats_index
            set_result_value(extractor,
                             result_dict, key_to_result_dict, stats_value)

        for stats_index, stats_value in df['mean_qscore'].describe().drop('count').items():
            key_to_result_dict = df_name + barcode_name + '.qscore.' + stats_index
            set_result_value(extractor,
                             result_dict, key_to_result_dict, stats_value)


def _barcode_frequency(extractor, barcode_selection, result_dict, entry: str, df_filtered) -> pd.Series:
    """
    Count reads by values of barcode_selection, computes sum of counts by barcode_selection, and sum of unclassified counts.
    Regroup all non used barcodes in index "other"
    Compute all frequency values for each number of barcoded reads
    :param result_dict: result dictionary with statistics
    :param entry: entry about barcoded counts
    :param prefix: key prefix
    :return: Series with all barcodes (used, non used, and unclassified) frequencies
    """
    # Regroup all barcoded read in Series
    all_barcode_count = df_filtered.value_counts()

    # Retain only existing barcodes from barcode_selection list
    barcodes_found = set(df_filtered.unique())
    barcode_selection_existing = [x for x in barcode_selection if x in barcodes_found]

    # Sort by list of barcode_selection
    count_sorted = all_barcode_count.sort_index()[barcode_selection_existing]
    # Replace all NaN values to zero
    count_sorted.fillna(0, downcast='int16', inplace=True)

    # Compute sum of all used barcodes without barcode 'unclassified'
    set_result_value(extractor, result_dict, entry + '.count', sum(count_sorted.drop("unclassified")))

    # Replace entry name ie read.pass/fail.barcode with read.pass/fail.non.used.barcodes.count
    non_used_barcodes_count_key = entry.replace(".barcoded", ".non.used.barcodes.count")

    # Compute all reads of barcodes that are not in the barcode_selection list
    other_barcode_count = sum(all_barcode_count) - sum(count_sorted)
    set_result_value(extractor, result_dict, non_used_barcodes_count_key, other_barcode_count)

    # Create Series for all non-used barcode counts and rename index array with "other"
    other_all_barcode_count = pd.Series(other_barcode_count, index=['other barcodes'])

    # Append Series of non-used barcode counts to the Series of barcode_selection counts
    count_sorted = count_sorted.append(other_all_barcode_count).sort_index()

    # Compute frequency for all barcode counts and save into dataframe_dict
    for barcode in count_sorted.to_dict():
        frequency_value = count_sorted[barcode] * 100 / sum(count_sorted)
        set_result_value(extractor, result_dict, entry.replace(".barcoded", ".") + barcode + ".frequency",
                         frequency_value)

    return count_sorted


def log_task(quiet, msg, start_time, end_time):
    if not quiet:
        delta = end_time - start_time
        print('  - {:} in {:}'.format(msg, common.format_duration(delta)))


def add_image_to_result(quiet, image_list, start_time, image):
    end_time = time.time()
    log_task(quiet, 'Creation of image "{0}"'.format(image[0]), start_time, end_time)


def read_first_line_file(filename):
    """
    Load the first line of a file.
    :param filename: name of the file to load.
    :return: the first line of the file
    """
    try:
        if filename.endswith('.gz'):
            with gzip.open(filename, 'rt') as f:
                return f.readline()
        elif filename.endswith('.bz2'):
            with bz2.open(filename, 'rt') as f:
                return f.readline()
        else:
            with open(filename, 'r') as f:
                return f.readline()
    except IOError:
        raise FileNotFoundError
