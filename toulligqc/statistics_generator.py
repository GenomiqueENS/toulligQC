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
#

#Creation of a text file containing statistics retrieved from FAST5, FASTQ files and sequencing_summary.txt file

import os
import numpy as np


def statistics_generator(config_dictionary, result_dict):
    '''
    Create a log file where different informations and statistics about the minion run are printed
    :param result_dict:
    :param config_dictionary:
    '''
    result_directory = config_dictionary['result_directory']
    barcode_selection = config_dictionary['barcode_selection']

    minknown_version = result_dict['minknow_version']
    flow_cell_id = result_dict['flow_cell_id']
    hostname = result_dict['hostname']
    minion_run_id = result_dict['minion_run_id']
    protocol_run_id = result_dict['protocol_run_id']


    completeName = os.path.join(result_directory+'statistics/', "run_statistics_file.txt")


    with open(completeName, 'w') as file_data:
        if config_dictionary['barcoding'] == True:
            for barcode in barcode_selection:
                if barcode == 'unclassified':
                    pass
                else:
                    fastq_length_statistics = result_dict['sequence_length_statistics_' + barcode]
                    mean_qscore_statistics = result_dict['mean_qscore_statistics_' + barcode]
                    barcode_total_nucleotide = result_dict['total_nucleotide_' + barcode]

                    for index, value in fastq_length_statistics.iteritems():
                        file_data.write(
                            "Read.fastq.length.{}.{}={}\n".format(index, barcode, np.round(value, decimals=2)))

                    nucleotide_counter = result_dict['nucleotide_count_' + barcode]
                    for nucleotide, count in nucleotide_counter.items():
                        file_data.write("nucleotide.{}.{}.template={}\n".format(nucleotide, barcode, np.round(count, decimals=2)))
                        calcul = float(count) / float(barcode_total_nucleotide)
                        calcul *= 100
                        file_data.write("nucleotide.{}.{}.proportion={}\n".format(nucleotide, barcode, np.round(calcul, decimals=2)))
                    file_data.write("barcode_total_nucleotide={}".format(barcode_total_nucleotide))

                    for key, mean_qscore_stat in mean_qscore_statistics.iteritems():
                        file_data.write("meanq_score.{}.{}={}".format(key, barcode, mean_qscore_stat))

        else:
            fastq_length_statistics = result_dict['sequence_length_statistics']
            nucleotide_counter = result_dict['nucleotide_count']
            mean_qscore_statistics = result_dict['mean_qscore_statistics']
            total_nucleotide = result_dict['total_nucleotide']


            for index, value in fastq_length_statistics.iteritems():
                file_data.write("Read.fastq.length.{}={}\n".format(index, np.round(value, decimals=2)))

            for nucleotide, count in nucleotide_counter.items():
                file_data.write(
                    "nucleotide.{}.template={}\n".format(nucleotide, np.round(count, decimals=2)))
                calcul = float(count) / float(total_nucleotide)
                calcul *= 100
                file_data.write("nucleotide.{}.proportion={}\n".format(nucleotide, np.round(calcul, decimals=2)))
            for key, mean_qscore_stat in mean_qscore_statistics.iteritems():
                file_data.write("meanq_score.{} ={}\n".format(key, mean_qscore_stat))

            file_data.write("barcode_total_nucleotide={}".format(total_nucleotide))

        channel_occupancy_statistics = result_dict['channel_occupancy_statistics']
        for index, value in channel_occupancy_statistics.iteritems():
            file_data.write("channel.occupancy.{}={}\n".format(index, value))

        file_data.write("Number.of.reads={}\n".format(len(result_dict['sequence_length_template'])))
        file_data.write("flowcell.serial.number={}\n".format(flow_cell_id))
        file_data.write("minknown.version={}\n".format(minknown_version))
        file_data.write("hostname={}\n".format(hostname))
        file_data.write("minion.serial.number={}\n".format(protocol_run_id))
        file_data.write(("run.id={}\n".format(minion_run_id)))

def save_result_file(config_dictionary, result_dict):
    result_directory = config_dictionary['result_directory']
    completeName = os.path.join(result_directory + 'statistics/', "save_result_statistics.txt")

    with open(completeName, 'w') as file_data:
        for key,value in result_dict.items():
            file_data.write('{}={}'.format(key, value))
