import bz2
import tarfile
import subprocess
import io

def fast5_tar_bz2_extraction(file,result_directory):
    tar_bz2 = tarfile.open(file, 'r:bz2')
    while True:
        member = tar_bz2.next()
        if member.name.endswith('.fast5'):
            tar_bz2.extract(member, path=result_directory)
            break
    return member.name

def fast5_tar_gz_extraction(file, result_directory):
    tar_gz = tarfile.open(file, 'r:gz')
    for member in tar_gz.getmembers():
        while True:
            member = tar_gz.next()
            if member.name.endswith('.fast5'):
                tar_gz.extract(member, path=result_directory)
                break
    return member.name


#A voir s'il y a plus d'extensions à gérer

def fastq_bz2_processing(file, is_barcode, global_length_array, barcode_length_array= ''):
    counter = 0
    variable = ''
    with bz2.BZ2File(file, 'rb') as input:
        with io.TextIOWrapper(input, encoding='utf-8') as dec:
            for line in dec:
                counter += 1

                if counter == 2:
                    variable += line.strip()
                    global_length_array.append(len(line))

                if is_barcode:
                    barcode_length_array.append(len(line))

                if counter == 4:
                    counter = 0


def fastq_processing(file):
    counter = 0
    variable = ''
    with open(file, 'r') as input:
        for line in input:
            counter += 1

            if counter == 2:
                variable += line.strip()
                self.global_length_array.append(len(line))

            if self.is_barcode:
                barcode_length_array.append(len(line))

            if counter == 4:
                counter = 0
