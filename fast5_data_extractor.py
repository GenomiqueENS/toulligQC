import h5py
import glob
import csv
import extraction


def minknow_version(h5py_file):
    """
    Get the Minknow version from fast5 file
    """
    version = list(h5py_file['/UniqueGlobalKey/tracking_id'].attrs.items())
    version_d = {key: value.decode('utf-8') for key, value in version}
    return version_d['version']

def flowcell_id(h5py_file):
    """
    Get the flowcell id from fast5 file
    """
    flowcell_id = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
    flowcell_id_dico = {key: value.decode('utf-7') for key, value in flowcell_id}
    return flowcell_id_dico['flow_cell_id']

def hostname(h5py_file):
    """
    Get the hostname from fast5 file
    """
    host_name = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
    host_name_dico = {key: value.decode('utf-8') for key, value in host_name}
    return host_name_dico['hostname']

def minion_run_id(h5py_file):
    """
    Get the number of Minion run
    """
    numMinION = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
    minion_run_id_dico = {key: value.decode('utf-8') for key, value in numMinION}
    return minion_run_id_dico['device_id']

def protocol_run_id(h5py_file):
    """
    Get the run id protocol from fast 5 file
    """
    protocol_run_id =  list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
    protocol_run_id_dico = {key: value.decode('utf-8') for key, value in protocol_run_id}
    return protocol_run_id_dico['protocol_run_id']




def fast5_data_extractor(fast5_source, result_directory, fast5_file_extension, run_name, config_file = ''):
    '''
    Extraction of different informations from a FAST5 file
    :param fast5_source: FAST5 file directory
    :param result_directory: result directory
    :param fast5_file_extension: extension used for the storage of the set of FAST5 files if there's one
    :param run_name: run name
    :return: a tuple containing the informations about a FAST5 file
    '''

    if fast5_file_extension == 'tar.bz2':
        if config_file:
            tar_bz2_file = fast5_source + run_name + ".tar.bz2"
        else:
            tar_bz2_file = fast5_source

        fast5_file = result_directory + extraction.fast5_tar_bz2_extraction(tar_bz2_file, result_directory)

    elif fast5_file_extension == 'tar.gz':
        if config_file:
            tar_gz_file = fast5_source + run_name + ".tar.gz"
        else:
            tar_gz_file = fast5_source
        fast5_file = result_directory + extraction.fast5_tar_gz_extraction(tar_gz_file, result_directory)

    elif fast5_file_extension == 'fast5_directory':
        fast5_file = glob.glob(fast5_source+"*.fast5")[0]

    else:
        fast5_file = fast5_source

    h5py_file = h5py.File(fast5_file)

    tuple_log_file = (flowcell_id(h5py_file), minknow_version(h5py_file), hostname(h5py_file), minion_run_id(h5py_file),\
                      minion_run_id(h5py_file))

    return tuple_log_file


def write_fast5_data_to_tsv(tuple_array):
    """
    Writes the data related to the fast5 files in a tsv file from the tuple array created by the fast5_data_extractor
    function
    """
    with open('fast5_file.tsv', 'w') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(['read_identifying', 'median_current_template', 'length_fastq_template','median_current_complement', 'length_fastq_complement', 'start_time','fastq_template','fastq_complement'])
        for row in tuple_array:
            writer.writerow(row)


def read_fast5_data_from_tsv(data_file):
    """
    Reads the tsv file containing the fast5 file data created previously by the write_data
    """
    with open(data_file) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        return reader
        for row in reader:
            print(', '.join(row))

