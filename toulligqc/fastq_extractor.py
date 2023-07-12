import multiprocessing as mp
from math import log
import os
import numpy as np
from tqdm import tqdm
import pandas as pd
import gzip
import time
from datetime import datetime
from toulligqc.sequencing_summary_common import log_task
from concurrent.futures import ProcessPoolExecutor, as_completed
from toulligqc.sequencing_summary_common import describe_dict
from toulligqc.sequencing_summary_common import set_result_value
from toulligqc.sequencing_summary_common import add_image_to_result
from toulligqc.sequencing_summary_common import check_result_values
from toulligqc.sequencing_summary_common import count_boolean_elements
from toulligqc.sequencing_summary_common import get_result_value
from toulligqc.sequencing_summary_common import series_cols_boolean_elements
from toulligqc import plotly_graph_generator as pgg
from toulligqc.common import is_numpy_1_24


class fastqExtractor:

    def __init__(self, config_dictionary):
        self.config_file_dictionary = config_dictionary
        self.fastq = config_dictionary['fastq'].split('\t')
        self.images_directory = config_dictionary['images_directory']
        self.threshold_Qscore = int(config_dictionary['threshold'])
        self.batch_size = int(config_dictionary['batch_size'])
        self.thread = int(config_dictionary['thread'])
        self.rich = False
        self.runid, self.sampleid, self.model_version_id = ['Unknow']*3
        if 'quiet' not in config_dictionary or config_dictionary['quiet'].lower() != 'true':
            self.quiet = False
        else:
            self.quiet = True

    def check_conf(self):
        """
        Configuration checking
        :return: nothing
        """
        for fastq in self.fastq:
            if not os.path.isfile(fastq):
                return False, "FASTQ file does not exists: " + fastq
            return True, ""


    def init(self):
        """
        Creation of the dataframe containing all info from fastq
        :return: Panda's Dataframe object
        """
        start_time = time.time()
        self.dataframe_1d = self._load_fastq_data()
        if self.dataframe_1d.empty:
            raise pd.errors.EmptyDataError("Dataframe is empty")
        self.dataframe_dict = {}
        log_task(self.quiet,
                 'Load FASTQ file ({:,.2f} MB used)'.format(self.dataframe_1d.memory_usage(deep=True).sum()/1024/1024),
                 start_time,
                 time.time())


    @staticmethod
    def get_name() -> str:
        """
        Get the name of the extractor.
        :return: the name of the extractor
        """
        return 'fastq'


    @staticmethod
    def get_report_data_file_id() -> str:
        """
        Get the report.data id of the extractor.
        :return: the report.data id
        """
        return 'basecaller.sequencing.summary.1d.extractor'


    def graph_generation(self, result_dict):
        """
        Generation of the different graphs containing in the plotly_graph_generator module
        :return: images array containing the title and the path toward the images
        """
        images = list()

        add_image_to_result(self.quiet, images, time.time(), pgg.read_count_histogram(result_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.read_length_scatterplot(self.dataframe_dict, self.images_directory))
        if self.rich:
            add_image_to_result(self.quiet, images, time.time(), pgg.yield_plot(self.dataframe_1d, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.read_quality_multiboxplot(self.dataframe_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.allphred_score_frequency(self.dataframe_dict, self.images_directory))
        if self.rich:
            add_image_to_result(self.quiet, images, time.time(), pgg.plot_performance(self.dataframe_1d, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.twod_density(self.dataframe_dict, self.images_directory))
        if self.rich:
            add_image_to_result(self.quiet, images, time.time(), pgg.sequence_length_over_time(self.dataframe_dict, self.images_directory))
            add_image_to_result(self.quiet, images, time.time(), pgg.phred_score_over_time(self.dataframe_dict, result_dict, self.images_directory))
        return images


    def clean(self, result_dict):
        """
        Removing dictionary entries that will not be kept in the report.data file
        :return:
        """
        check_result_values(self, result_dict)
        self.dataframe_dict.clear()
        self.dataframe_1d.iloc[0:0]


    def _fastq_read_batch(self, read_batch):
        """
        parse each line of FASTQ quality line:
        return: [read length, mean Qscore, type of read (pass or fail)]
        """
        #rich = False if isinstance(read_batch[0], str) else True
        fastq_lines = []
        if self.rich:
            for read in read_batch:
                name = read[0]
                qscore = self._avg_qual(read[1])
                passes_filtering = True if qscore > self.threshold_Qscore else False
                start_time, ch = self._extract_info_from_name(name)
                fastq_lines.append((len(read[1]), qscore, passes_filtering, start_time, ch))
        else:
            for read in read_batch:
                qscore = self._avg_qual(read)
                passes_filtering = True if qscore > self.threshold_Qscore else False
                fastq_lines.append((len(read), qscore, passes_filtering))
        return fastq_lines


    def _fastq_batch_generator(self):
        """
        read FASTQ file in small batch
        yeald : list of lines (quality line): batch of n size
        """
        #if len(self.fastq) == 1:
            #open_fn = gzip.open if self.fastq.endswith('.gz') else open
        #else:
        #open_fn = gzip.open if self.fastq[0].endswith('.gz') else open
        for fastq in self.fastq:
            open_fn = gzip.open if fastq.endswith('.gz') else open
            with open_fn(fastq, 'rt') as fastq:
                batch = []
                line_num = 0
                if self.rich:
                    for line in fastq:
                        line_num += 1
                        if line_num % 4 == 1:
                            name = line.rstrip()
                        elif line_num % 4 == 0:
                            batch.append((name, line.rstrip()))
                        if len(batch) == self.batch_size:
                            yield batch
                            batch = []
                else:       
                    for line in fastq:
                        line_num += 1
                        if line_num % 4 == 0:
                            batch.append(line.rstrip())
                            if len(batch) == self.batch_size:
                                yield batch
                                batch = []
                if len(batch) > 0:
                    yield batch


    def extract(self, result_dict):
        """
        Get Phred score (Qscore) and Length details (frequencies, ratios, yield and statistics) per type read (pass or fail)
        :param result_dict:
        """
        start_time = time.time()
        self._fill_series_dict(self.dataframe_dict, self.dataframe_1d)

        self._set_result_dict_telemetry_value(result_dict, "run.id", self.runid)
        self._set_result_dict_telemetry_value(result_dict, "sample.id", self.sampleid)
        self._set_result_dict_telemetry_value(result_dict, "model.file", self.model_version_id)

        set_result_value(self, result_dict, "read.count", len(self.dataframe_1d))
        set_result_value(self, result_dict, "read.pass.count",
                         count_boolean_elements(self.dataframe_1d, 'passes_filtering', True))
        set_result_value(self, result_dict, "read.fail.count",
                         count_boolean_elements(self.dataframe_1d, 'passes_filtering', False))
        total_reads = get_result_value(self, result_dict, "read.count")

        # Ratios
        set_result_value(self, result_dict, "read.pass.ratio",
                         (get_result_value(self, result_dict, "read.pass.count") / total_reads))
        set_result_value(self, result_dict, "read.fail.ratio",
                         (get_result_value(self, result_dict, "read.fail.count") / total_reads))

        # Frequencies
        set_result_value(self, result_dict, "read.count.frequency", 100)

        read_pass_frequency = (get_result_value(self,
                                                result_dict, "read.pass.count") / total_reads) * 100
        set_result_value(self,
                         result_dict, "read.pass.frequency", read_pass_frequency)

        read_fail_frequency = (get_result_value(self,
                                                result_dict, "read.fail.count") / total_reads) * 100
        set_result_value(self,
                         result_dict, "read.fail.frequency", read_fail_frequency)

        # Yield, n50, run time
        set_result_value(self, result_dict, "yield", sum(self.dataframe_dict["all.reads.sequence.length"]))

        set_result_value(self, result_dict, "n50", self._compute_NXX(50))
        set_result_value(self, result_dict, "l50", self._compute_LXX(50))

        if self.rich:
            set_result_value(self, result_dict, "run.time", max(self.dataframe_1d['start_time']))
            # Get channel occupancy statistics and store each value into result_dict
            for index, value in self._occupancy_channel().items():
                set_result_value(self,
                                result_dict, "channel.occupancy.statistics." + index, value)
        
        # Get statistics about all reads length and store each value into result_dict
        sequence_length_statistics = self.dataframe_dict["all.reads.sequence.length"].describe()

        for index, value in sequence_length_statistics.items():
            set_result_value(self,
                             result_dict, "all.read.length." + index, value)

        # Add statistics (without count) about read pass/fail length in the result_dict
        describe_dict(self, result_dict, self.dataframe_dict["pass.reads.sequence.length"],
                      "pass.reads.sequence.length")
        describe_dict(self, result_dict, self.dataframe_dict["fail.reads.sequence.length"],
                      "fail.reads.sequence.length")

        # Get Qscore statistics without count value and store them into result_dict
        qscore_statistics = self.dataframe_1d['mean_qscore'].describe().drop(
            "count")

        for index, value in qscore_statistics.items():
            set_result_value(self,
                             result_dict, "all.read.qscore." + index, value)

        # Add statistics (without count) about read pass/fail qscore in the result_dict
        describe_dict(self, result_dict, self.dataframe_dict["pass.reads.mean.qscore"], "pass.reads.mean.qscore")
        describe_dict(self, result_dict, self.dataframe_dict["fail.reads.mean.qscore"], "fail.reads.mean.qscore")

        log_task(self.quiet, 'Extract info from FASTQ file', start_time, time.time())       


    def _load_fastq_data(self):
        """
        Load FASTQ dataframe
        :return: a Pandas Dataframe object
        """
        run_info = self.check_fastq()

        if run_info:
            self.rich = True
            self.runid, self.sampleid, self.model_version_id = run_info   
        else:
            self.rich = False

        read_batchs = self._fastq_batch_generator()
        rst_futures = multiprocessing_submit(self._fastq_read_batch,
                                                        read_batchs, 
                                                        n_process=self.thread, 
                                                        pbar_update=self.batch_size)
        fq_data = []
        
        for _, f in enumerate(rst_futures):
            fq_data.extend(f.result())

        columns = ['sequence_length', 'mean_qscore', 'passes_filtering', 'start_time', 'channel'] if self.rich \
                else ['sequence_length', 'mean_qscore', 'passes_filtering']

        fq_data = pd.DataFrame(fq_data, columns=columns)

        fq_data['sequence_length'] = fq_data['sequence_length'].astype(np.uint32)
        fq_data['mean_qscore'] = fq_data['mean_qscore'].astype(np.float32)
        fq_data['passes_filtering'] = fq_data['passes_filtering'].astype(np.bool)

        if self.rich:
            fq_data["start_time"] = fq_data["start_time"] - fq_data["start_time"].min()
            fq_data['channel'] = fq_data['channel'].astype(np.int16)
            fq_data['start_time'] = fq_data['start_time'].astype(np.float64)

        return fq_data 


    def _fill_series_dict(self, df_dict, df):
        """
        """
        for read_type in ['pass', 'fail']:
            read_type_bool = True if read_type == 'pass' else False

            df_dict[read_type + '.reads.sequence.length'] = series_cols_boolean_elements(df,
                                                                                         'sequence_length',
                                                                                         'passes_filtering',
                                                                                         read_type_bool)
            df_dict[read_type + '.reads.mean.qscore'] = series_cols_boolean_elements(df,
                                                                                     'mean_qscore',
                                                                                     'passes_filtering',
                                                                                     read_type_bool)
        df_dict["all.reads.sequence.length"] = df['sequence_length']
        df_dict["all.reads.mean.qscore"] = df['mean_qscore']
        if self.rich:
            df_dict["all.reads.start.time"] = df['start_time']


    def _compute_LXX(self, x):
        """Compute LXX value of total sequence length"""
        data = self.dataframe_dict["all.reads.sequence.length"].dropna().values
        data.sort()
        half_sum = data.sum() * x / 100
        cum_sum = 0
        count = 0
        for v in data:
            cum_sum += v
            count += 1
            if cum_sum >= half_sum:
                return count
            

    def _compute_NXX(self, x):
        """Compute NXX value of total sequence length"""
        data = self.dataframe_dict["all.reads.sequence.length"].dropna().values
        data.sort()
        half_sum = data.sum() * x / 100
        cum_sum = 0
        for v in data:
            cum_sum += v
            if cum_sum >= half_sum:
                return int(v)


    def _occupancy_channel(self):
        """
        Statistics about the channels of the flowcell
        :return: pd.Series object containing statistics about the channel occupancy without count value
        """
        total_reads_per_channel = pd.value_counts(self.dataframe_1d["channel"])
        return pd.DataFrame.describe(total_reads_per_channel)


    def _extract_info_from_name(self, name):
        """
        
        """
        metadata = dict(x.split("=") for x in name.split(" ")[1:])
        start_time = self._timeISO_to_float(metadata['start_time'])

        return  start_time, metadata['ch']


    def _timeISO_to_float(self, iso_datetime):
        """
        
        """
        dt = datetime.strptime(iso_datetime, '%Y-%m-%dT%H:%M:%SZ')
        unix_timestamp = dt.timestamp()
        return unix_timestamp


    def _avg_qual(self, quals):
        """
        Estimates mean quality Phred score
        return: float
        """
        if quals:
            qscore =  -10 * log(sum([10**((ord(q)-33) / -10) for q in quals]) / len(quals), 10)
            return round(qscore, 2)
        else:
            return None 


    def check_fastq(self):
        """
        """
        #fastq = self.fastq if len(self.fastq) == 1 else self.fastq[0]
        open_fn = gzip.open if self.fastq[0].endswith('.gz') else open

        with open_fn(self.fastq[0], 'rt') as fq:
            first_line = fq.readline().strip('\n')
        try:
            metadata = dict(x.split("=") for x in first_line.split(" ")[1:])
            if 'model_version_id' not in metadata:
                metadata['model_version_id'] = 'Unknow'
            return metadata['runid'] , metadata['sampleid'] , metadata['model_version_id'] 
        except:
            return None
        
    def _set_result_dict_telemetry_value(self, result_dict, key, new_value):
        """
        """
        final_key = "sequencing.telemetry.extractor." + key
        current_value = None

        if final_key in result_dict:
            current_value = result_dict[final_key]
            if len(current_value) == 0:
                current_value = None

        if new_value is None:
            new_value = current_value

        result_dict[final_key] = new_value


def multiprocessing_submit(func, iterator, n_process=mp.cpu_count()-1 ,pbar = True, pbar_update = 500,  *arg, **kwargs):
    executor = ProcessPoolExecutor(n_process)

    max_queue = n_process * 2
    if pbar:
        pbar = tqdm(unit = 'read', desc='Processed')

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


def batch_iterator(iterator, batch_size):
    batch = []
    i=0
    for entry in iterator:
        i += 1
        batch.append(entry)
        if i == batch_size:
            yield batch
            batch = []
            i = 0
    if len(batch):
        yield batch