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

import os
import numpy as np


def log_file1D(config_dictionary, result_dict):
    '''
    Create a log file like where different information about the minion run are printed
    :param fast5_data: tuple containing informations from a raw FAST5 file
    :param basecall_stat: basecalling_stat_plotter instance
    :param result_directory: result directory
    '''
    result_directory = config_dictionary['result_directory']
    barcode_selection = config_dictionary['barcode_selection']

    minknown_version = result_dict['minknow_version']
    flowcell_id = result_dict['flowcell_id']
    hostname = result_dict['hostname']
    minion_run_id = result_dict['minion_run_id']
    protocol_run_id = result_dict['protocol_run_id']


    completeName = os.path.join(result_directory+'statistics/', "run_statistics_file.txt")


    with open(completeName, 'w') as file_data:
        if config_dictionary['barcoding']:
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
                        calcul = calcul * 100
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
                calcul = calcul * 100
                file_data.write("nucleotide.{}.proportion={}\n".format(nucleotide, np.round(calcul, decimals=2)))
            for key, mean_qscore_stat in mean_qscore_statistics.iteritems():
                file_data.write("meanq_score.{} ={}\n".format(key, mean_qscore_stat))

            file_data.write("barcode_total_nucleotide={}".format(total_nucleotide))

        channel_occupancy_statistics = result_dict['channel_occupancy_statistics']
        for index, value in channel_occupancy_statistics.iteritems():
            file_data.write("channel.occupancy.{}={}\n".format(index, value))

        file_data.write("Number.of.reads={}\n".format(len(result_dict['sequence_length_template'])))
        file_data.write("flowcell.serial.number={}\n".format(flowcell_id))
        file_data.write("minknown.version={}\n".format(minknown_version))
        file_data.write("hostname={}\n".format(hostname))
        file_data.write("minion.serial.number={}\n".format(protocol_run_id))
        file_data.write(("run.id={}\n".format(minion_run_id)))

# def statistics_dataframe(self, result_dict):
#
#     """
#     Returns the statistics retrieved from the statistics files in the statistics directory for each barcode as a dataframe to make
#     the reading easier.
#     """
#
#     nucleotide = ['A', 'T', 'C', 'G']
#     sequence_length_matrix = [[0 for column in range(len(self.barcode_selection) + 1)] for line in range(8)]
#     channel_occupancy_matrix = [[0 for column in range(len(self.barcode_selection) + 1)] for line in range(8)]
#     mean_qscore_matrix = [[0 for column in range(len(self.barcode_selection) + 1)] for line in range(8)]
#     nucleotide_count_matrix = [[0 for column in range(len(self.barcode_selection) + 1)] for line in range(4)]
#     nucleotide_proportion_matrix = [[0 for column in range(len(self.barcode_selection) + 1)] for line in range(4)]
#
#
#     column = 0
#     for index_barcode, barcode in enumerate(self.barcode_selection):
#        # barcode_selected_dataframe = self.albacore_log[self.albacore_log['barcode_arrangement'] == barcode]
#         #channel_occupancy_statistics = barcode_selected_dataframe['channel'].describe()
#         #mean_qscore_statistics = barcode_selected_dataframe['mean_qscore_template'].describe()
#         #sequence_length_statistics = barcode_selected_dataframe['sequence_length_template'].describe()
#         if barcode != 'unclassified':
#             sorted_list = sorted(result_dict['nucleotide_count_'+barcode].items())
#         if column == 0:
#             for i in range(8):
#                 mean_qscore_matrix[i][column] = 'phred_score_' + result_dict['mean_qscore_statistics_' + barcode].keys()[i]
#                 #channel_occupancy_matrix[i][0] = 'channel_occupancy_' + result_dict['mean_qscore_statistics_' + barcode].keys()[i]
#                 sequence_length_matrix[i][0] = 'sequence_length_' + sequence_length_statistics.keys()[i]
#
#             for line, nucleotide_number_list in enumerate(sorted_list):
#                 if nucleotide_number_list[0] in nucleotide and barcode != 'unclassified':
#                     nucleotide_count_matrix[line][column] = 'nucleotide_count_' + nucleotide_number_list[0]
#                     nucleotide_proportion_matrix[line][column] = 'nucleotide_proportion_'+ nucleotide_number_list[0]
#                 else:
#                     continue
#
#             column += 1
#
#         for line, metric in enumerate(mean_qscore_statistics):
#             mean_qscore_matrix[line][column] = round(metric, 3)
#
#         for line, metric in enumerate(channel_occupancy_statistics):
#             channel_occupancy_matrix[line][column] = round(metric, 3)
#
#         for line, metric in enumerate(sequence_length_statistics):
#             sequence_length_matrix[line][column] = round(metric, 3)
#
#         for line, nucleotide_number_list in enumerate(sorted_list):
#             if nucleotide_number_list[0] in nucleotide and barcode != 'unclassified':
#                 nucleotide_count_matrix[line][column] = float(nucleotide_number_list[1])
#                 calcul = nucleotide_number_list[1] / result_dict['total_nucleotide_'+barcode]
#                 nucleotide_proportion_matrix[line][column] = round(calcul, 3)
#             else:
#                 continue
#         column += 1
#
#     barcode_selection_matrix = list(self.barcode_selection)
#     barcode_selection_matrix.insert(0, '')
#     result_dict["minknown.version"] = self.minknown_version
#     result_dict["hostname"] = self.hostname
#     result_dict["minion.serial.number"] = self.numMinion
#     result_dict["run.id"] = self.run_id
#     general_information_list = [['minknown', result_dict["minknown.version"]],
#                                 ['hostname', result_dict["hostname"]],
#                                 ['minion.serial.number', result_dict["minion.serial.number"]],
#                                 ['run.id', result_dict["run.id"]]]
#
#     with open(self.result_directory + 'dataframe.csv', 'a') as csv_file:
#         writer = csv.writer(csv_file, delimiter='\t')
#
#         writer.writerow(['Number of reads:', len(self.sequence_length_template)])
#         writer.writerow('')
#
#         writer.writerow(barcode_selection_matrix)
#         for element in channel_occupancy_matrix:
#             writer.writerow(element)
#
#         writer.writerow('')
#         for metric in mean_qscore_matrix:
#             writer.writerow(metric)
#
#         writer.writerow('')
#         for metric in sequence_length_matrix:
#             writer.writerow(metric)
#
#         writer.writerow('')
#         for metric in nucleotide_count_matrix:
#             writer.writerow(metric)
#
#         writer.writerow('')
#         for metric in nucleotide_proportion_matrix:
#             writer.writerow(metric)
#
#         writer.writerow('')
#         writer.writerows(general_information_list)


# def log_file_tsv(fast5_data , basecall_stat, result_directory):
#     '''
#     Creation of a tsv file summarising the statistics got
#     :param fast5_data: tuple containing informations from a raw FAST5 file
#     :param basecall_stat: basecalling_stat_plotter instance
#     :param result_directory: result directory
#     '''
#     version, flowcell_id, hostname, minion_run_id, protocol_run_id = fast5_data
#     num_called_template, mean_qscore_template = basecall_stat.stat_generation()
#     counter_template, total_nucleotide_template = basecall_stat.counter()
#     occupancy_pore = basecall_stat.occupancy_pore()
#
#     with open(result_directory + 'dataframe.csv', 'a') as tsv_file:
#
#         writer = csv.writer(tsv_file, delimiter='\t')
#
#         for index, element in num_called_template.iteritems():
#             writer.writerow(['num.called.template_'+index, element])
#
#         for index, element in num_called_template.iteritems():
#             writer.writerow(["mean.qscore.template_"+index, np.round(element, decimals=2)])
#
#         for nucleotide, count in counter_template.items():
#
#             writer.writerow(["nucleotide.{}".format(nucleotide),float(count)])
#             if nucleotide == 'total':
#                 continue
#             calcul = float(count) / float(total_nucleotide_template)
#             writer.writerow(["nucleotide.{}.proportion".format(nucleotide), np.round(calcul, decimals=2)])
#
#         for index, value in occupancy_pore.items():
#             writer.writerow(["channel.occupancy.{}".format(index),value])
#
#
#         writer.writerow(["number.of.reads",float(len(basecall_stat.sequence_length_template))])
#         writer.writerow(["flowcell.serial.number",flowcell_id])
#         writer.writerow(["minknown.version",version])
#         writer.writerow(["hostname",hostname])
#         writer.writerow(["minion.serial.number",numMinion])
#         writer.writerow(["run.id",run_id])
