import h5py
import glob
import numpy as np

from collections import Counter
import csv
import getter1D
import os
import subprocess
import tarfile
import parser
def fast5_data_extractor(fast5_file_directory):
    """
    Creates a dataframe from a collection of fast5 files. 
    Needs the fast5 file directory as input (where fast5 files are stored) and returns a tuple with a set of information
    about the fast5 files.
    """
    run_name,selected_file,is_docker,is_barcode = parser.get_args()
    dico_path = parser.file_path_initialization()
    fast5_file = glob.glob(dico_path['design_file_directory']+'*.fast5')[0]
    if not fast5_file:
        print('Pas de fichier fast5')
        return 0
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
