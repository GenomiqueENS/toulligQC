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

import pandas as pd


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
