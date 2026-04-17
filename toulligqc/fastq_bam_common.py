#
#                  ToulligQC development code
#
# This code may be freely distributed and modified under the
# terms of the GNU General Public License version 3 or later
# and CeCILL. This should be distributed with the code. If you
# do not have a copy, see:
#
#      http://www.gnu.org/licenses/gpl-3.0-standalone.html
#      http://www.cecill.info/licences/Licence_CeCILL_V2-en.html
#
# Copyright for this code is held jointly by the Genomic platform
# of the Institut de Biologie de l'École Normale Supérieure and
# the individual authors.
#
# For more information on the ToulligQC project and its aims,
# visit the home page at:
#
#      https://github.com/GenomiqueENS/toulligQC
#

import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed

from tqdm import tqdm


def extract_headerTag(header, tagGroup, tag, defaultValue=None):

    if tagGroup not in header:
        if defaultValue is not None:
            return defaultValue
        else:
            raise KeyError(tagGroup)

    first_entry = header[tagGroup][0]

    if tag not in first_entry:
        if defaultValue is not None:
            return defaultValue
        else:
            raise KeyError(tag)

    return first_entry[tag]


def batch_iterator(iterator, batch_size):
    batch = []
    i = 0
    for entry in iterator:
        i += 1
        batch.append(entry.to_string())
        if i == batch_size:
            yield batch
            batch = []
            i = 0
    if len(batch):
        yield batch


def multiprocessing_submit(
    func,
    iterator,
    n_process=mp.cpu_count() - 1,
    pbar=True,
    pbar_update=500,
    *arg,
    **kwargs,
):

    executor = ProcessPoolExecutor(n_process)

    max_queue = n_process * 2
    if pbar:
        pbar = tqdm(unit="read", desc="Processed")

    futures = {}
    n_job_in_queue = 0
    while True:
        while n_job_in_queue < max_queue:
            i = next(iterator, None)
            if not i:
                break
            futures[executor.submit(func, i, *arg, **kwargs)] = None
            n_job_in_queue += 1

        job = next(as_completed(futures), None)

        # no more job
        if job is None:
            break
        else:
            n_job_in_queue -= 1
            pbar.update(pbar_update)
            yield job
            del futures[job]
