from Bio.SeqIO.QualityIO import FastqGeneralIterator

import csv
import re
import matplotlib
matplotlib.use('Agg')
import parser


def get_MinknowVersion(h5py_file):
    """
    Get the Minknow version from fast5 file
    """
    version = list(h5py_file['/UniqueGlobalKey/tracking_id'].attrs.items())
    version_d = {key: value.decode('utf-8') for key, value in version}
    return version_d['version']
#print(get_MinknowVersion())

def get_FlowcellId(h5py_file):
    """
    Get the flowcell id from fast5 file
    """
    flowcell_id = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
    flowcell_id_dico = {key: value.decode('utf-7') for key, value in flowcell_id}
    return flowcell_id_dico['flow_cell_id']

def get_Hostname(h5py_file):
    """
    Get the hostname from fast5 file
    """
    host_name = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
    host_name_dico = {key: value.decode('utf-8') for key, value in host_name}
    return host_name_dico['hostname']

def get_MinIONRunId(h5py_file):
    """
    Get the number of Minion run
    """
    numMinION = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
    numMinION_dico = {key: value.decode('utf-8') for key, value in numMinION}
    return numMinION_dico['device_id']

def get_ProtocolRunId(h5py_file):
    """
    Get the run id protocol from fast 5 file
    """
    protocol_run_id =  list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
    protocol_run_id_dico = {key: value.decode('utf-8') for key, value in protocol_run_id}
    return protocol_run_id_dico['protocol_run_id']

def get_barcode():
    """
    Get the barcode from a file given in input
    """
    dico_path = parser.file_path_initialization()
    barcode_file = dico_path['design_file_directory']+"design.csv"
    #barcode_file = input('Enter the name of the file where the barcode is stored:')
    set_doublon = set()

    with open(barcode_file) as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t')

        for row in spamreader:
            pattern = re.search(r'BC(\d{2})', row[0])

            if pattern:
                barcode = 'barcode{}'.format(pattern.group(1))
                set_doublon.add(barcode)
    return list(set_doublon)
