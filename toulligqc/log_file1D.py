import os
import numpy as np
import csv

def log_file1D(fast5_data , basecall_stat, result_directory):
    '''
    Create a log file like where different information about the minion run are printed
    :param fast5_data: tuple containing informations from a raw FAST5 file
    :param basecall_stat: basecalling_stat_plotter instance
    :param result_directory: result directory
    '''
    version, flowcell_id, hostname, numMinion, run_id = fast5_data

    num_called_template, mean_qscore_template = basecall_stat.stat_generation()

    counter_template, total_nucleotide_template = basecall_stat.counter()

    occupancy_pore = basecall_stat.occupancy_pore()
    completeName = os.path.join(result_directory+'statistics/', "run_statistics_file.txt")


    with open(completeName, 'w') as file_data:

        for index, element in num_called_template.iteritems():
            file_data.write("num.called.template.{}={}\n".format(index, element))

        for index, element in mean_qscore_template.iteritems():
            file_data.write("mean.qscore.template.{}={}\n".format(index, np.round(element, decimals=2)))

        for nucleotide, count in counter_template.items():
            file_data.write("nucleotide.{}.template={}\n".format(nucleotide,count))
            if nucleotide == 'total':
                continue
            calcul = float(count) / float(total_nucleotide_template)
            file_data.write("nucleotide.{}.proportion={}\n".format(nucleotide, np.round(calcul, decimals=2)))


        file_data.write("total.number.of.sequence={}\n".format(basecall_stat.fast5_tot))

        for index, value in occupancy_pore.items():
            file_data.write("channel.occupancy.{}={}\n".format(index, value))

        file_data.write("Number.of.reads={}\n".format(len(basecall_stat.sequence_length_template)))
        file_data.write("flowcell.serial.number={}\n".format(flowcell_id))
        file_data.write("minknown.version={}\n".format(version))
        file_data.write("hostname={}\n".format(hostname))
        file_data.write("minion.serial.number={}\n".format(numMinion))
        file_data.write(("run.id={}\n".format(run_id)))


        for index, element in basecall_stat.statistics_read_size().iteritems():
            file_data.write("Read.fastq.length.{}={}\n".format(index, np.round(element, decimals=2)))

def log_file_tsv(fast5_data , basecall_stat, result_directory):
    '''
    Creation of a tsv file summarising the statistics got
    :param fast5_data: tuple containing informations from a raw FAST5 file
    :param basecall_stat: basecalling_stat_plotter instance
    :param result_directory: result directory
    '''
    version, flowcell_id, hostname, numMinion, run_id = fast5_data
    num_called_template, mean_qscore_template = basecall_stat.stat_generation()
    counter_template, total_nucleotide_template = basecall_stat.counter()
    occupancy_pore = basecall_stat.occupancy_pore()

    with open(result_directory + 'dataframe.csv', 'a') as tsv_file:

        writer = csv.writer(tsv_file, delimiter='\t')

        for index, element in num_called_template.iteritems():
            writer.writerow(['num.called.template_'+index, element])

        for index, element in num_called_template.iteritems():
            writer.writerow(["mean.qscore.template_"+index, np.round(element, decimals=2)])

        for nucleotide, count in counter_template.items():

            writer.writerow(["nucleotide.{}".format(nucleotide),float(count)])
            if nucleotide == 'total':
                continue
            calcul = float(count) / float(total_nucleotide_template)
            writer.writerow(["nucleotide.{}.proportion".format(nucleotide), np.round(calcul, decimals=2)])

        for index, value in occupancy_pore.items():
            writer.writerow(["channel.occupancy.{}".format(index),value])


        writer.writerow(["number.of.reads",float(len(basecall_stat.sequence_length_template))])
        writer.writerow(["flowcell.serial.number",flowcell_id])
        writer.writerow(["minknown.version",version])
        writer.writerow(["hostname",hostname])
        writer.writerow(["minion.serial.number",numMinion])
        writer.writerow(["run.id",run_id])
