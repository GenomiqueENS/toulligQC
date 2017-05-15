import os

def log_file1D(fast5_data , basecall_stat):
    """
    Create a log file like aozan where different information are printed
    """

    version, flowcell_id, hostname, numMinion, run_id = fast5_data

    #Retrieve the dataframe with statitstics such as the quartile or std
    #Retrieve the dictionary from albacore summary log

    num_called_template, mean_qscore_template = basecall_stat.stat_generation()

    counter_template, total_nucleotide_template = basecall_stat.counter()

    occupancy_pore = basecall_stat.occupancy_pore()

    completeName = os.path.join('/home/ferrato/Documents/fast5', "fichier_aozan.txt")

    with open(completeName, 'w') as file_data:

        for index, element in num_called_template.iteritems():
            file_data.write("num.called.template.{}={}\n".format(index, element))

        for index, element in num_called_template.iteritems():
            file_data.write("mean.qscore.template.{}={}\n".format(index, element))

        for nucleotide, count in counter_template.items():
            file_data.write("nucleotide.{}.template={}\n".format(nucleotide,count))
            if nucleotide == 'total':
                continue
            calcul = float(count) / float(total_nucleotide_template)
            file_data.write("nucleotide.{}.proportion={}\n".format(nucleotide, calcul))


        file_data.write("total.number.of.sequence={}\n".format(basecall_stat.fast5_tot))

        for index, value in occupancy_pore.items():
            file_data.write("pore.occupancy.{}={}\n".format(index, value))


        file_data.write("flowcell.serial.number={}\n".format(flowcell_id))
        file_data.write("minknown.version={}\n".format(version))
        file_data.write("hostname={}\n".format(hostname))
        file_data.write("minion.serial.number={}\n".format(numMinion))
        file_data.write(("run.id={}\n".format(run_id)))

        for index, element in basecall_stat.statistics_read_size().iteritems():
            file_data.write("Read.fastq.length.{}={}\n".format(index, element))
#log_file1D("/home/ferrato/shares-net/sequencages/nanopore/test_alabacore/save/workspace", "/home/ferrato/shares-net/sequencages/nanopore/test_alabacore/save/sequencing_summary.txt")



