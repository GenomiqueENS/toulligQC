import argparse
import configparser

def get_args():
    parser = argparse.ArgumentParser()
    # Add more options if you like
    parser.add_argument("-n", action='store', dest="run_name", help="name of the sample file",default=True)
    parser.add_argument("-d",action='store_true',dest='docker',help="docker use",default=False)
    parser.add_argument("-b",action='store_true',dest='is_barcode',help="docker use",default=False)
    parser.add_argument('-s','--arg',nargs='+',type=str)
    argument_list=vars(parser.parse_args())
    argument_value=parser.parse_args()
    run_name = argument_value.run_name
    selected_file = argument_list["arg"]
    is_docker = argument_value.docker
    is_barcode = argument_value.is_barcode
    return run_name,selected_file,is_docker,is_barcode

def file_path_initialization():
    dico_path = {}
    run_name,selected_file,_ ,is_barcode = get_args()
    configFilePath = r'config.txt'

    config = configparser.ConfigParser()
    config.read(configFilePath)

    dico_path['result_directory'] = config.get('config', 'result.directory')
    dico_path['design_file_directory'] = config.get('config', 'design.file.directory')
    dico_path['basecall_log'] = config.get('config', 'log.file')
    dico_path['fastq_directory'] = config.get('config', 'fastq.directory')
    dico_path['fast5_directory'] = config.get('config', 'fast5.directory')
    for key, value in dico_path.items():
        if value.endswith('/'):
           continue
        else:
            dico_path[key] = value+'/'

    return dico_path

