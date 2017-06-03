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


def get_Barcodes(selected_file_list=''):
    """
    Gets the barcode from a file given in input
    """
    #delete barcode doublon in selected file list 
    if selected_file_list:
        barcode_set = set()
        for file in selected_file_list:
            pattern = re.search(r'barcode(\d{2})', file)

            if pattern:
                barcode = 'barcode{}'.format(pattern.group(1))
                barcode_set.add(barcode)
        return list(barcode_set)

    else:
        configParser = configparser.ConfigParser()
        
        if os.path.isfile('/configpass/docker_config.txt'):
            barcode_design_file = "/design.file.directory/design.csv"

        else:
            configFilePath = r'config.txt'
            configParser.read(configFilePath)

            barcode_design_file = configParser.get('config', 'design.file.directory')
            if barcode_design_file.endswith('/'):
                barcode_design_file = barcode_design_file+"design.csv"
            else:
                barcode_design_file = barcode_design_file+"/"+"design.csv"
        barcode_set = set()

        with open(barcode_design_file) as csvfile:
            design_file = csv.reader(csvfile, delimiter='\t')

            for row in design_file:
                pattern = re.search(r'BC(\d{2})', row[0]):

                if pattern:
                    barcode = 'barcode{}'.format(pattern.group(1))
                    barcode_set.add(barcode)
        return list(barcode_set)

def get_FastqSeq_barcoded(barcode_selection, run_name, barcode_present, selected_file_list=''):

    """
    Gets the fastq sequence for a barcode selection
    """
    configParser = configparser.ConfigParser()

    
    if os.path.isfile('/configpass/docker_config.txt'):
        bz2_fastq_directory = "/bz2.fastq.directory/" + run_name

    else:
        configFilePath = r'config.txt'
        configParser.read(configFilePath)
        bz2_fastq_directory = configParser.get('config', 'bz2.fastq.directory') + run_name
        if bz2_fastq_directory.endswith('/'):
            bz2_fastq_directory = bz2_fastq_directory+run_name
        else:
            bz2_fastq_directory = bz2_fastq_directory+'/'+run_name 

    global_length_array = []

    #fastq_sequence.txt represent the file where the fastq sequence extracted from bz2 files are written

    if os.path.isfile('/configpass/docker_config.txt'): 
        fastq_file = '/working.directory/fastq_sequence.txt'
    else:
        fastq_file = 'fastq_sequence.txt'
        
    # bz2_fastq_directory = input('path to bz2 directory:')

    if os.path.isfile(fastq_file):
        open(fastq_file, 'w').close()

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


    #represent a directory including a set of files barcoded 

    if not selected_file_list:

        for bz2_fastq_file in glob.glob("{}/*.bz2".format(bz2_fastq_directory)):

            template_nucleotide_counter = Counter()
            total_nucs_template = 0
            selected_barcoded_sample_fastq_length_array = []

            for selected_barcode in barcode_selection:
                if selected_barcode in bz2_fastq_file:
                    print(bz2_fastq_file)
                    uncompressedData = bz2.BZ2File(bz2_fastq_file, 'rb').read()
                    uncomp = uncompressedData.decode('utf-8')
                    file = open(fastq_file, 'w')
                    file.write(uncomp)

                    with open(fastq_file) as in_handle:

                        for title, seq, qual in FastqGeneralIterator(in_handle):
                            selected_barcoded_sample_fastq_length_array.append(len(seq))
                            global_length_array.append(len(seq))

                            for template_nucleotide in seq:
                                template_nucleotide_counter[template_nucleotide] += 1
                                total_nucs_template += 1

                    if os.path.isfile('/configpass/docker_config.txt'): 
                        completeName = os.path.join('/working.directory/statistics/', selected_barcode)
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


    #represent selected file list with barcodes

    else selected_file_list:
        if barcode_present == 'y':
            template_nucleotide_counter = Counter()
            total_nucs_template = 0
            selected_barcoded_sample_fastq_length_array = []
            for file in selected_file_list:
                file = bz2_fastq_directory + '/' + file
                for selected_barcode in barcode_selection:
                    if selected_barcode in file:
                        print(file)
                        uncompressedData = bz2.BZ2File(file, 'rb').read()
                        uncomp = uncompressedData.decode('utf-8')
                        file = open(fastq_file, 'w')
                        file.write(uncomp)

                        with open(fastq_file) as in_handle:

                            for title, seq, qual in FastqGeneralIterator(in_handle):
                                selected_barcoded_sample_fastq_length_array.append(len(seq))
                                global_length_array.append(len(seq))

                                for template_nucleotide in seq:
                                    template_nucleotide_counter[template_nucleotide] += 1
                                    total_nucs_template += 1

                        if os.path.isfile('/configpass/docker_config.txt'):
                            completeName = os.path.join('/working.directory/statistics/', selected_barcode)

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

def get_FastqSeq_without_barcode(run_name, selected_file_list=''):
    """
    Gets the fastq sequence
    """
    configParser = configparser.ConfigParser()

    
    if os.path.isfile('/configpass/docker_config.txt'):
        bz2_fastq_directory = "/bz2.fastq.directory/" + run_name

    else:
        configFilePath = r'config.txt'
        configParser.read(configFilePath)
        bz2_fastq_directory = configParser.get('config', 'bz2.fastq.directory') 
        if bz2_fastq_directory.endswith('/'):
            bz2_fastq_directory+run_name
        else:
            bz2_fastq_directory+'/'+run_name

    global_length_array = []

    if os.path.isfile('/configpass/docker_config.txt'): 
        fastq_file = '/working.directory/fastq_sequence.txt'
    else:
        fastq_file = 'fastq_sequence.txt'
        
    # bz2_fastq_directory = input('path to bz2 directory:')

    if os.path.isfile(fastq_file):
        open(fastq_file, 'w').close()

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
    #selected file list without barcode
    if selected_file_list:

        for file in selected_file_list:
            file = bz2_fastq_directory + '/' + file
            template_nucleotide_counter = Counter()
            total_nucs_template = 0

            uncompressedData = bz2.BZ2File(file, 'rb').read()
            uncomp = uncompressedData.decode('utf-8')
            file = open(fastq_file, 'w')
            file.write(uncomp)

            with open(fastq_file) as in_handle:

                for title, seq, qual in FastqGeneralIterator(in_handle):
                    global_length_array.append(len(seq))

                    for template_nucleotide in seq:
                        template_nucleotide_counter[template_nucleotide] += 1
                        total_nucs_template += 1
        return template_nucleotide_counter, total_nucs_template, global_length_array

    #directory with samples without barcode(a set of files not selected)
    else:
        for bz2_fastq_file in glob.glob("{}/*.bz2".format(bz2_fastq_directory)):

            template_nucleotide_counter = Counter()
            total_nucs_template = 0

            uncompressedData = bz2.BZ2File(bz2_fastq_file, 'rb').read()
            uncomp = uncompressedData.decode('utf-8')
            file = open(fastq_file, 'w')
            file.write(uncomp)

            with open(fastq_file) as in_handle:

                for title, seq, qual in FastqGeneralIterator(in_handle):
                    global_length_array.append(len(seq))

                    for template_nucleotide in seq:
                        template_nucleotide_counter[template_nucleotide] += 1
                        total_nucs_template += 1
        return template_nucleotide_counter, total_nucs_template, global_length_array






