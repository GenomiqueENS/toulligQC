import h5py
import glob
import numpy as np

from collections import Counter
import csv
import getter1D


def fast5_data_extractor(fast5_file_directory):
    """Create a dataframe from collections of fast5 files
        :param fast5_file_directory: directory where fast5 files are stored
        :return : tuple with different informations about fast5 files
        """

    for fast5_file in glob.glob("{}/*.fast5".format(fast5_file_directory)):

        try:
            h5py_file = h5py.File(fast5_file)
            # version
            version = getter1D.get_MinknowVersion(h5py_file)

            # flowcell_id
            flowcell_id = getter1D.getFlowcellId(h5py_file)

            # hostname
            hostname = getter1D.get_Hostname(h5py_file)

            # numMinion
            numMinion = getter1D.getNumMinION(h5py_file)

            # run_id
            run_id = getter1D.getProtocolRunId(h5py_file)

        except:
            continue

        break

    tuple_log_file = (flowcell_id, version , hostname,numMinion,run_id)

    return tuple_log_file


def write_data(tuple_array):
    """
    Write data in a tsv file
    """
    with open('fast5_file.tsv', 'w') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(['read_identifying', 'median_current_template', 'length_fastq_template','median_current_complement', 'length_fastq_complement', 'start_time','fastq_template','fastq_complement'])
        for row in tuple_array:
            writer.writerow(row)


def read_data(data_file):
    """
    Read the file created previously
    """
    with open(data_file) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        return reader
        for row in reader:
            print(', '.join(row))

