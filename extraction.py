import tarfile

def fast5_tar_bz2_extraction(tar_bz2_file,result_directory):
    '''
    Extraction of a FAST5 file from a set of FAST5 files
    :param tar_bz2_file: tar bz2 file containing the set of the raw FAST5 files
    :param result_directory: result directory
    :return: a FAST5 file
    '''
    tar_bz2 = tarfile.open(tar_bz2_file, 'r:bz2')
    while True:
        member = tar_bz2.next()
        if member.name.endswith('.fast5'):
            tar_bz2.extract(member, path=result_directory)
            break
    return member.name

def fast5_tar_gz_extraction(tar_gz_file, result_directory):
    '''
    Extraction of a FAST5 file from a set of FAST5 files
    :param tar_gz_file: tar gz file containing the set of the raw FAST5 files
    :param result_directory: result directory
    :return: a FAST5 file
    '''
    tar_gz = tarfile.open(tar_gz_file, 'r:gz')
    for member in tar_gz.getmembers():
        while True:
            member = tar_gz.next()
            if member.name.endswith('.fast5'):
                tar_gz.extract(member, path=result_directory)
                break
    return member.name




