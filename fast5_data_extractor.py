import h5py
import glob
import csv
import getter1D
import parser
import extraction


def fast5_data_extractor(fast5_file_directory, result_directory, dico_extension):
    """
    Creates a dataframe from a collection of fast5 files. 
    Needs the fast5 file directory as input (where fast5 files are stored) and returns a tuple with a set of information
    about the fast5 files.
    """

    run_name, selected_file, is_docker, is_barcode = parser.get_args()
    if dico_extension['fast5_file_extension'] == 'tar.bz2':
        tar_bz2_file = fast5_file_directory + run_name + ".tar.bz2"
        fast5_file = result_directory + extraction.fast5_tar_bz2_extraction(tar_bz2_file, result_directory)

    elif dico_extension['fast5_file_extension'] == 'tar.gz':
        tar_gz_file = fast5_file_directory + run_name + ".tar.bz2"
        fast5_file = result_directory + extraction.fast5_tar_gz_extraction(tar_gz_file, result_directory)

    else:
        fast5_file = glob.glob(fast5_file_directory+"*.fast5")[0]

    print(fast5_file)
    h5py_file = h5py.File(fast5_file)

    # version
    version = getter1D.get_MinknowVersion(h5py_file)

    # flowcell_id
    flowcell_id = getter1D.get_FlowcellId(h5py_file)

    # hostname
    hostname = getter1D.get_Hostname(h5py_file)

    # numMinion
    numMinion = getter1D.get_MinIONRunId(h5py_file)

    # run_id
    run_id = getter1D.get_ProtocolRunId(h5py_file)

    tuple_log_file = (flowcell_id, version , hostname,numMinion,run_id)

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

