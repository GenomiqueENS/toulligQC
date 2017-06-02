from Bio.SeqIO.QualityIO import FastqGeneralIterator
import os
import bz2
import csv
import re
import glob
import matplotlib

matplotlib.use('Agg')
from matplotlib import pyplot as plt
from collections import Counter
import pandas as pd
import configparser


def get_MinknowVersion(h5py_file):
    """
    Gets the Minknow version from fast5 file
    """
    version = list(h5py_file['/UniqueGlobalKey/tracking_id'].attrs.items())
    version_d = {key: value.decode('utf-8') for key, value in version}
    return version_d['version']


# print(get_MinknowVersion())

def get_FlowcellId(h5py_file):
    """
    Gets the flowcell id from fast5 file
    """
    flowcell_id = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
    flowcell_id_dico = {key: value.decode('utf-8') for key, value in flowcell_id}
    return flowcell_id_dico['flow_cell_id']


def get_Hostname(h5py_file):
    """
    Gets the hostname from fast5 file
    """
    host_name = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
    host_name_dico = {key: value.decode('utf-8') for key, value in host_name}
    return host_name_dico['hostname']


def get_MinIONRunId(h5py_file):
    """
    Gets the id of Minion run
    """
    numMinION = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
    numMinION_dico = {key: value.decode('utf-8') for key, value in numMinION}
    return numMinION_dico['device_id']


def get_ProtocolRunId(h5py_file):
    """
    Gets the run id protocol from fast 5 file
    """
    protocol_run_id = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
    protocol_run_id_dico = {key: value.decode('utf-8') for key, value in protocol_run_id}
    return protocol_run_id_dico['protocol_run_id']


def get_Barcodes(file_list=''):
    """
    Gets the barcode from a file given in input
    """
    if file_list:
        barcode_set = set()
        for file in file_list:
            pattern = re.search(r'barcode(\d{2})', file)

            if pattern:
                barcode = 'barcode{}'.format(pattern.group(1))
                barcode_set.add(barcode)
        return list(barcode_set)

    else:
        configParser = configparser.ConfigParser()

        try:
            configFilePath = r'/configpass/docker_config.txt'
            configParser.read(configFilePath)
            barcode_design_file = "/design.file.directory/design.csv"

        except:
            configFilePath = r'config.txt'
            configParser.read(configFilePath)
            barcode_design_file = configParser.get('config', 'design.file.directory') + "design.csv"

        barcode_set = set()

        with open(barcode_design_file) as csvfile:
            design_file = csv.reader(csvfile, delimiter='\t')

            for row in design_file:
                pattern = re.search(r'BC(\d{2})', row[0])

                if pattern:
                    barcode = 'barcode{}'.format(pattern.group(1))
                    barcode_set.add(barcode)
        return list(barcode_set)


# selection : barcode name list
def get_FastqSeq(barcode_selection, run_name, barcode_present, file_list=''):
    """
    Gets the fastq sequence for a barcode selection
    """
    configParser = configparser.ConfigParser()

    try:
        configFilePath = r'/configpass/docker_config.txt'
        configParser.read(configFilePath)
        path_bz2_directory = "/bz2.fastq.directory/" + run_name


    except:
        configFilePath = r'config.txt'
        configParser.read(configFilePath)
        path_bz2_directory = configParser.get('config', 'bz2.fastq.directory') + run_name

    global_length_array = []
    if configFilePath == r'/configpass/docker_config.txt':
        path_bz2_file = '/design.file.directory/fastq_sequence.txt'
    else:
        path_bz2_file = 'fastq_sequence.txt'
        
    # path_bz2_directory = input('path to bz2 directory:')

    if os.path.isfile(path_bz2_file):
        open(path_bz2_file, 'w').close()

    if not os.path.exists('images'):
        os.makedirs('images')

    if os.path.exists('images'):
        for file in os.listdir('images'):
            os.remove('images/' + file)

    if not os.path.exists('statistics'):
        os.makedirs('statistics')

    if os.path.exists('statistics'):
        for file in os.listdir('statistics'):
            os.remove('statistics/' + file)

    if barcode_present == 'y':

        for bz2_file in glob.glob("{}/*.bz2".format(path_bz2_directory)):

            template_nucleotide_counter = Counter()
            total_nucs_template = 0
            selected_barcoded_sample_fastq_length_array = []

            for selected_barcode in barcode_selection:
                if selected_barcode in bz2_file:
                    print(bz2_file)
                    uncompressedData = bz2.BZ2File(bz2_file, 'rb').read()
                    uncomp = uncompressedData.decode('utf-8')
                    file = open(path_bz2_file, 'w')
                    file.write(uncomp)

                    with open(path_bz2_file) as in_handle:

                        for title, seq, qual in FastqGeneralIterator(in_handle):
                            selected_barcoded_sample_fastq_length_array.append(len(seq))
                            global_length_array.append(len(seq))

                            for template_nucleotide in seq:
                                template_nucleotide_counter[template_nucleotide] += 1
                                total_nucs_template += 1
                    if configFilePath == r'/configpass/docker_config.txt':
                        completeName = os.path.join('/design.file.directory/statistics/', selected_barcode)
                    else:
                        completeName = os.path.join('statistics/', selected_barcode)
                        
                    barcode_file = open(completeName, 'w')

                    selected_barcode_fastq_size = pd.Series(selected_barcoded_sample_fastq_length_array)
                    selected_barcode_fastq_size_statistics = pd.Series.describe(selected_barcode_fastq_size)

                    for index, value in selected_barcode_fastq_size_statistics.iteritems():
                        barcode_file.write("Read.fastq.length.{}={}\n".format(index, value))

                    for nucleotide, nucleotide_proportion in template_nucleotide_counter.items():
                        barcode_file.write("nucleotide.{}.template={}\n".format(nucleotide, nucleotide_proportion))
                        if nucleotide == 'total':
                            continue
                        calcul = float(nucleotide_proportion) / float(total_nucs_template)
                        barcode_file.write("nucleotide.{}.proportion={}\n".format(nucleotide, calcul))
                    barcode_file.close()

                    plt.boxplot(selected_barcoded_sample_fastq_length_array, showfliers=False)
                    plt.title('Read length boxplot for {}'.format(selected_barcode))
                    plt.savefig('images/image_{}.png'.format(selected_barcode))
                    plt.close()
                    file.close()
        return global_length_array

    elif file_list:
        if barcode_present == 'y':
            template_nucleotide_counter = Counter()
            total_nucs_template = 0
            selected_barcoded_sample_fastq_length_array = []
            for file in file_list:
                file = path_bz2_directory + '/' + file
                for selected_barcode in barcode_selection:
                    if selected_barcode in file:
                        print(file)
                        uncompressedData = bz2.BZ2File(file, 'rb').read()
                        uncomp = uncompressedData.decode('utf-8')
                        file = open(path_bz2_file, 'w')
                        file.write(uncomp)

                        with open(path_bz2_file) as in_handle:

                            for title, seq, qual in FastqGeneralIterator(in_handle):
                                selected_barcoded_sample_fastq_length_array.append(len(seq))
                                global_length_array.append(len(seq))

                                for template_nucleotide in seq:
                                    template_nucleotide_counter[template_nucleotide] += 1
                                    total_nucs_template += 1
                        if configFilePath == r'/configpass/docker_config.txt':
                            completeName = os.path.join('/design.file.directory/statistics/', selected_barcode)
                        else:
                            completeName = os.path.join('statistics/', selected_barcode)
                        barcode_file = open(completeName, 'w')

                        selected_barcode_fastq_size = pd.Series(selected_barcoded_sample_fastq_length_array)
                        selected_barcode_fastq_size_statistics = pd.Series.describe(selected_barcode_fastq_size)

                        for index, value in selected_barcode_fastq_size_statistics.iteritems():
                            barcode_file.write("Read.fastq.length.{}={}\n".format(index, value))

                        for nucleotide, nucleotide_proportion in template_nucleotide_counter.items():
                            barcode_file.write("nucleotide.{}.template={}\n".format(nucleotide, nucleotide_proportion))
                            if nucleotide == 'total':
                                continue
                            calcul = float(nucleotide_proportion) / float(total_nucs_template)
                            barcode_file.write("nucleotide.{}.proportion={}\n".format(nucleotide, calcul))
                        barcode_file.close()

                        plt.boxplot(selected_barcoded_sample_fastq_length_array, showfliers=False)
                        plt.title('Read length boxplot for {}'.format(selected_barcode))
                        plt.savefig('images/image_{}.png'.format(selected_barcode))
                        plt.close()
                        file.close()
            return global_length_array

        else:

            for file in file_list:
                file = path_bz2_directory + '/' + file
                template_nucleotide_counter = Counter()
                total_nucs_template = 0

                uncompressedData = bz2.BZ2File(file, 'rb').read()
                uncomp = uncompressedData.decode('utf-8')
                file = open(path_bz2_file, 'w')
                file.write(uncomp)

                with open(path_bz2_file) as in_handle:

                    for title, seq, qual in FastqGeneralIterator(in_handle):
                        global_length_array.append(len(seq))

                        for template_nucleotide in seq:
                            template_nucleotide_counter[template_nucleotide] += 1
                            total_nucs_template += 1
            return template_nucleotide_counter, total_nucs_template, global_length_array


    else:
        for bz2_file in glob.glob("{}/*.bz2".format(path_bz2_directory)):

            template_nucleotide_counter = Counter()
            total_nucs_template = 0

            uncompressedData = bz2.BZ2File(bz2_file, 'rb').read()
            uncomp = uncompressedData.decode('utf-8')
            file = open(path_bz2_file, 'w')
            file.write(uncomp)

            with open(path_bz2_file) as in_handle:

                for title, seq, qual in FastqGeneralIterator(in_handle):
                    global_length_array.append(len(seq))

                    for template_nucleotide in seq:
                        template_nucleotide_counter[template_nucleotide] += 1
                        total_nucs_template += 1
        return template_nucleotide_counter, total_nucs_template, global_length_array







