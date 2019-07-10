# -*- coding: utf-8 -*-
#
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

# Creation of a text file containing statistics retrieved from the result_dict dictionary.
# The information by modules is retained in a key-value form.
# The prefix of a key being the report data file id of the module.

import os


def add_values_to_unwritten_key(result_dict, values):
    """
    :param result_dict:
    :param values: must be a list
    :return:
    """
    return result_dict['unwritten.keys'].extend(values)


def statistics_generator(config_dictionary, result_dict):
    """
    Create a log file where different information and statistics about the minion run are printed
    :param result_dict:
    :param config_dictionary:
    """
    result_directory = config_dictionary['result_directory']

    complete_name = os.path.join(result_directory, "report.data")

    add_values_to_unwritten_key(result_dict, ["sequence_length_template", "passes_filtering",
                                              "read_pass", "read_fail", "start_time_sorted",
                                              "read_pass_sorted", "read_fail_sorted"])

    with open(complete_name, 'w') as file_data:
        for key, value in result_dict.items():
            if key not in result_dict['unwritten.keys']:
                file_data.write("{0}={1}\n".format(key, value))
