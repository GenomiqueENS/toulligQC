from math import log
import os
import numpy as np
import pandas as pd
import time
import pysam
from datetime import datetime
from toulligqc.extractor_common import log_task
from toulligqc.extractor_common import describe_dict
from toulligqc.extractor_common import set_result_value
from toulligqc.extractor_common import add_image_to_result
from toulligqc.extractor_common import check_result_values
from toulligqc.extractor_common import count_boolean_elements
from toulligqc.extractor_common import get_result_value
from toulligqc.extractor_common import set_result_dict_telemetry_value
from toulligqc.extractor_common import fill_series_dict
from toulligqc.extractor_common import timeISO_to_float
from toulligqc.common_statistics import compute_NXX, compute_LXX, occupancy_channel, avg_qual
from toulligqc.fastq_bam_common import multiprocessing_submit, extract_headerTag
from toulligqc.fastq_bam_common import batch_iterator
from toulligqc import plotly_graph_generator as pgg


class uBAM_Extractor:
    def __init__(self, config_dictionary):
        self.config_file_dictionary = config_dictionary
        self.ubam = config_dictionary['bam'].split('\t')
        self.images_directory = config_dictionary['images_directory']
        self.threshold_Qscore = int(config_dictionary['threshold'])
        self.batch_size = int(config_dictionary['batch_size'])
        self.thread = int(config_dictionary['thread'])
        self.header = dict()
        if 'quiet' not in config_dictionary or config_dictionary['quiet'].lower() != 'true':
            self.quiet = False
        else:
            self.quiet = True

    def check_conf(self):
        """
        Configuration checking
        :return: nothing
        """
        for uBAM in self.ubam:
            if not os.path.isfile(uBAM):
                return False, "BAM file does not exists: " + uBAM
            return True, ""


    def init(self):
        """
        Creation of the dataframe containing all info from uBAM
        :return: Panda's Dataframe object
        """
        start_time = time.time()
        self.dataframe_1d = self._load_uBAM_data()
        if self.dataframe_1d.empty:
            raise pd.errors.EmptyDataError("Dataframe is empty")
        self.dataframe_dict = {}
        log_task(self.quiet,
                 'Load BAM file ({:,.2f} MB used)'.format(self.dataframe_1d.memory_usage(deep=True).sum()/1024/1024),
                 start_time,
                 time.time())


    def clean(self, result_dict):
        """
        Removing dictionary entries that will not be kept in the report.data file
        :return:
        """
        check_result_values(self, result_dict)
        self.dataframe_dict.clear()
        self.dataframe_1d.iloc[0:0]


    @staticmethod
    def get_name() -> str:
        """
        Get the name of the extractor.
        :return: the name of the extractor
        """
        return 'ubam'


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
        add_image_to_result(self.quiet, images, time.time(), pgg.yield_plot(self.dataframe_1d, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.read_quality_multiboxplot(self.dataframe_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.allphred_score_frequency(self.dataframe_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.plot_performance(self.dataframe_1d, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.twod_density(self.dataframe_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.sequence_length_over_time(self.dataframe_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.phred_score_over_time(self.dataframe_dict, result_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.speed_over_time(self.dataframe_dict, self.images_directory))
        return images


    def extract(self, result_dict):
        """
        Get Phred score (Qscore) and Length details (frequencies, ratios, yield and statistics) per type read (pass or fail)
        :param result_dict:
        """
        start_time = time.time()
        fill_series_dict(self.dataframe_dict, self.dataframe_1d)

        set_result_dict_telemetry_value(result_dict, "run.id", self.header["run_id"])
        set_result_dict_telemetry_value(result_dict, "sample.id", self.header["sample_id"])
        set_result_dict_telemetry_value(result_dict, "model.file", self.header["model_version_id"])
        set_result_dict_telemetry_value(result_dict, "software.name", self.header["basecaller"])
        set_result_dict_telemetry_value(result_dict, "software.version", self.header["basecaller_version"])
        set_result_dict_telemetry_value(result_dict, "flowcell.id", self.header["flow_cell_id"])
        set_result_dict_telemetry_value(result_dict, "basecalling.date", self.header["run_date"])
        set_result_dict_telemetry_value(result_dict, "pass.threshold.qscore", str(self.threshold_Qscore))

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

        set_result_value(self, result_dict, "n50", compute_NXX(self.dataframe_dict, 50))
        set_result_value(self, result_dict, "l50", compute_LXX(self.dataframe_dict, 50))

        set_result_value(self, result_dict, "run.time", max(self.dataframe_1d['start_time']))
        # Get channel occupancy statistics and store each value into result_dict
        for index, value in occupancy_channel(self.dataframe_1d).items():
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

        log_task(self.quiet, 'Extract info from uBAM file', start_time, time.time())       


    def _load_uBAM_data(self):
        """
        Load uBAM dataframe
        :return: a Pandas Dataframe object
        """  
        self._get_header()
        uBAM_chunks = self._uBAM_batch_generator()
        rst_futures = multiprocessing_submit(self._uBAM_batch_reader,
                                                        uBAM_chunks, 
                                                        n_process=self.thread, 
                                                        pbar_update=self.batch_size)
        uBAM_data = []
        
        for _, f in enumerate(rst_futures):
            uBAM_data.extend(f.result())

        columns = ['sequence_length', 'mean_qscore', 'passes_filtering', 'start_time', 'channel', 'duration']
        
        uBAM_data = pd.DataFrame(uBAM_data, columns=columns)

        uBAM_data['sequence_length'] = uBAM_data['sequence_length'].astype(np.uint32)
        uBAM_data['mean_qscore'] = uBAM_data['mean_qscore'].astype(np.float32)
        uBAM_data['passes_filtering'] = uBAM_data['passes_filtering'].astype(np.bool)
        uBAM_data["start_time"] = uBAM_data["start_time"] - uBAM_data["start_time"].min()
        uBAM_data['channel'] = uBAM_data['channel'].astype(np.int16)
        uBAM_data['start_time'] = uBAM_data['start_time'].astype(np.float64)
        uBAM_data['duration'] = uBAM_data['duration'].astype(np.float32)
        return uBAM_data 


    def _uBAM_batch_reader(self, uBAM_chunk):
        """
        parse each line of uBAM quality line:
        return: [read length, mean Qscore, type of read (pass or fail)]
        """
        #def process_bam_chunk(bam_chunk):
        rec_data = []
        for rec in uBAM_chunk:
            rec_dict = self._process_record(rec)
            rec_data.append(rec_dict)
        return rec_data


    def _uBAM_batch_generator(self):
        """
        read uBAM file in small chunk
        yield : list of lines (quality line): batch of n size
        """
        for ubam in self.ubam:
            samfile = pysam.AlignmentFile(ubam, "rb", check_sq=False)
            bam_batch = batch_iterator(samfile, batch_size=self.batch_size)
            for batch in bam_batch:
                yield batch


    def _timeISO_to_float(self, iso_datetime, format):
        """
        """
        dt = datetime.strptime(iso_datetime, format)
        unix_timestamp = dt.timestamp()
        return unix_timestamp


    def _get_header(self):
        samfile = pysam.AlignmentFile(self.ubam[0], "rb", check_sq=False)
        header = samfile.header.to_dict()
        run_id, model_version_id =  extract_headerTag(header,'RG','ID').split('_', 1)
        self.header = {
        "run_id" : run_id,
        "run_date" : extract_headerTag(header, 'RG', 'DT'),
        "sample_id" : extract_headerTag(header,'RG','SM'),
        "basecaller" : extract_headerTag(header,'PG','PN'),
        "basecaller_version" : extract_headerTag(header,'PG','VN'),
        "model_version_id" : model_version_id,
        "flow_cell_id" : extract_headerTag(header,'RG','PU')
        }


    def _process_record(self, rec):
        """
        extract QC info from BAM record
        return : dict of QC info
        """
        tags = rec.split("\t")
        iso_start_time = tags[17].split(':',2)[2]
        qual = avg_qual(tags[10])
        passes_filtering = True if qual > self.threshold_Qscore else False
        data = [
            len(tags[9]), # read length
            qual, # AVG Qscore
            passes_filtering, # Passing filter
            timeISO_to_float(iso_start_time, '%Y-%m-%dT%H:%M:%S.%f%z'), # start time
            tags[16].split(':',2)[2], # Channel
            tags[12].split(':',2)[2] # Duration
        ]
        return data