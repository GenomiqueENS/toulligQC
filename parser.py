import argparse
import configparser
import os

def get_args():
    parser = argparse.ArgumentParser()
    # Add more options if you like
    parser.add_argument("-n", "--run_name", action='store', dest="run_name", help="Name of the sample file",default=True)
    parser.add_argument("-b","--barcode", action='store_true',dest='is_barcode',help="Barcode usage",default=False)
    parser.add_argument("-c", "--config_file", action='store', dest='config_file', help="Configuration file", default=False)
    parser.add_argument('-f', '--arg', nargs='+', type=str, help="Path to directory without config file in the same order that the config file")
    argument_list = vars(parser.parse_args())
    file_list = argument_list["arg"]
    argument_value=parser.parse_args()
    run_name = argument_value.run_name
    config_file = argument_value.config_file
    is_barcode = argument_value.is_barcode
    return run_name, config_file, is_barcode, file_list

def file_path_initialization(config_file, is_barcode, run_name):

    dico_path = {}
    configFilePath = config_file
    config = configparser.ConfigParser()
    config.read(configFilePath)
    print(config.get('config', 'result.directory'))
    dico_path['result_directory'] = config.get('config', 'result.directory')
    dico_path['basecall_log'] = config.get('config', 'log.file')
    dico_path['fastq_directory'] = config.get('config', 'fastq.directory')
    dico_path['fast5_directory'] = config.get('config', 'fast5.directory')

    if is_barcode:
        dico_path['design_file_directory'] = config.get('config', 'design.file.directory')

    for key, value in dico_path.items():
        if value.endswith('/'):
            continue
        else:
            dico_path[key] = value + '/'

    dico_path['result_directory'] = dico_path['result_directory'] + run_name + '/'

    return dico_path

def config_file_initialization(config_file, is_barcode, run_name, list_file=''):
    if list_file:
        dico_path = {}
        dico_path['result_directory'] = list_file[0]
        dico_path['basecall_log'] = list_file[1]
        dico_path['fastq_directory'] = list_file[2]
        dico_path['fast5_directory'] = list_file[3]

        if is_barcode:
            dico_path['design_file_directory'] = list_file[4]

    elif config_file:
        dico_path = file_path_initialization(config_file, is_barcode, run_name)

    elif os.path.isfile('~/.toulligqc'):
        dico_path = file_path_initialization('~/.toulligqc',is_barcode, run_name)

    else:
        print('Error, not a config file')
        return 0

    return dico_path


def extension(config_file, list_file):

    dico_extension = {}
    if list_file:
        dico_extension['fast5_file_extension'] = list_file[5]
        dico_extension['fastq_file_extension'] = list_file[6]
    else:
        config = configparser.ConfigParser()
        config.read(config_file)
        dico_extension['fast5_file_extension'] = config.get('extension', 'fast5.file.extension')
        dico_extension['fastq_file_extension'] = config.get('extension', 'fastq.file.extension')

    return dico_extension