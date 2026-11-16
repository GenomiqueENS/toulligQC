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
from collections.abc import Callable, Iterator
from concurrent.futures import ProcessPoolExecutor, as_completed

from tqdm import tqdm


def extract_headerTag(header: dict, tagGroup: str, tag: str, first_only=False):

    if tagGroup not in header:
        return None

    result_list = []
    for v in header[tagGroup]:
        if tag in v and v[tag] not in result_list:
            result_list.append(v[tag])

    if len(result_list) == 0:
        return None

    if first_only:
        return result_list[0]
    return ", ".join(result_list)


def extract_headerTag_to_dict(header: dict, tagGroup: str, tag: str) -> dict:

    result = {}

    if tagGroup not in header:
        return result

    tag_group_dict = header[tagGroup][0]
    if tag not in tag_group_dict:
        return result

    value = tag_group_dict[tag]
    for v in value.split(" "):
        k, v = v.split("=", 1)
        result[k] = v

    return result


def batch_iterator(iterator, batch_size: int) -> Iterator[list[str]]:
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
    func: Callable,
    iterator,
    n_process: int = mp.cpu_count() - 1,
    pbar: bool = True,
    pbar_update: int = 500,
    *arg,
    **kwargs,
) -> Iterator:

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
