
import os
import numpy as np
import pandas as pd
import gzip
import time
from toulligqc.extractor_common import log_task
from toulligqc.extractor_common import describe_dict
from toulligqc.extractor_common import set_result_value
from toulligqc.extractor_common import check_result_values
from toulligqc.extractor_common import add_image_to_result
from toulligqc.extractor_common import count_boolean_elements
from toulligqc.extractor_common import get_result_value
from toulligqc.extractor_common import fill_series_dict
from toulligqc.extractor_common import set_result_dict_telemetry_value
from toulligqc.extractor_common import timeISO_to_float
from toulligqc.extractor_common import extract_barcode_info
from toulligqc.common_statistics import compute_NXX, compute_LXX, occupancy_channel, avg_qual
from toulligqc.fastq_bam_common import multiprocessing_submit
from toulligqc.common import is_numpy_1_24
from toulligqc import plotly_graph_generator as pgg


class fastqExtractor:

    def __init__(self, config_dictionary):
        self.config_dictionary = config_dictionary
        self.fastq = config_dictionary['fastq'].split('\t')
        self.images_directory = config_dictionary['images_directory']
        self.threshold_Qscore = int(config_dictionary['threshold'])
        self.batch_size = int(config_dictionary['batch_size'])
        self.thread = int(config_dictionary['thread'])
        self.rich = False
        self.runid, self.sampleid, self.model_version_id = ['Unknow']*3
        self.is_barcode = False
        if config_dictionary['barcoding'] == 'True':
            self.is_barcode = True
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

        # Add missing categories
        if 'barcode_arrangement' in self.dataframe_1d.columns:
            self.dataframe_1d['barcode_arrangement'] = self.dataframe_1d['barcode_arrangement'].cat.add_categories([0, 
                                                                                                    'other barcodes',
                                                                                                    'passes_filtering'])
        self.dataframe_1d = self.dataframe_1d.fillna(0)
        self.barcode_selection = self.config_dictionary['barcode_selection']
        

        log_task(self.quiet,
                 'Load FASTQ file ({:,.2f} MB used)'.format(self.dataframe_1d.memory_usage(deep=True).sum()/1024/1024),
                 start_time,
                 time.time())
    
    def clean(extractor, result_dict):
        """
        Removing dictionary entries that will not be kept in the report.data file
        :return:
        """

        # Check values in result_dict (avoid Series and Dataframes)
        check_result_values(extractor, result_dict)

        # Clear dictionary for Series and Dataframe
        extractor.dataframe_dict.clear()

        # Clear DataFrame
        extractor.dataframe_1d.iloc[0:0]


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
        if self.is_barcode:
            add_image_to_result(self.quiet, images, time.time(), pgg.barcode_percentage_pie_chart_pass(self.dataframe_dict,
                                                                                                       self.barcode_selection,
                                                                                                       self.images_directory))

            read_fail = self.dataframe_dict["read.fail.barcoded"]
            if not (len(read_fail) == 1 and read_fail["other barcodes"] == 0):
                add_image_to_result(self.quiet, images, time.time(), pgg.barcode_percentage_pie_chart_fail(self.dataframe_dict,
                                                                                                      self.barcode_selection,
                                                                                                      self.images_directory))

            add_image_to_result(self.quiet, images, time.time(), pgg.barcode_length_boxplot(self.dataframe_dict,
                                                                                            self.images_directory))

            add_image_to_result(self.quiet, images, time.time(), pgg.barcoded_phred_score_frequency(self.dataframe_dict,
                                                                                                    self.images_directory))
        return images


    def extract(self, result_dict):
        """
        Get Phred score (Qscore) and Length details (frequencies, ratios, yield and statistics) per type read (pass or fail)
        :param result_dict:
        """
        start_time = time.time()
        fill_series_dict(self.dataframe_dict, self.dataframe_1d)

        set_result_dict_telemetry_value(result_dict, "run.id", self.runid)
        set_result_dict_telemetry_value(result_dict, "sample.id", self.sampleid)
        set_result_dict_telemetry_value(result_dict, "model.file", self.model_version_id)

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

        if self.rich:
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
        if self.is_barcode:
            extract_barcode_info(self, result_dict,
                                 self.barcode_selection,
                                 self.dataframe_dict,
                                 self.dataframe_1d)

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
        rst_futures = multiprocessing_submit(self._fastq_batch_reader,
                                                        read_batchs, 
                                                        n_process=self.thread, 
                                                        pbar_update=self.batch_size)
        fq_df = []
        
        for _, f in enumerate(rst_futures):
            fq_df.extend(f.result())

        columns = ['sequence_length', 'mean_qscore', 'passes_filtering']
        if self.rich:
            columns.extend(['start_time', 'channel'])
        if self.is_barcode:
            columns.append('barcode_arrangement')

        fq_df = pd.DataFrame(fq_df, columns=columns)

        fq_df['sequence_length'] = fq_df['sequence_length'].astype(np.uint32)
        fq_df['mean_qscore'] = fq_df['mean_qscore'].astype(np.float32)
        fq_df['passes_filtering'] = fq_df['passes_filtering'].astype(np.bool_ if is_numpy_1_24 else np.bool)

        if self.rich:
            fq_df["start_time"] = fq_df["start_time"] - fq_df["start_time"].min()
            fq_df['start_time'] = fq_df['start_time'].astype(np.float64)
            fq_df['channel'] = fq_df['channel'].astype(np.int16)
        if self.is_barcode:
            fq_df['barcode_arrangement'] = fq_df['barcode_arrangement'].astype("category")
        return fq_df 


    def _fastq_batch_generator(self):
        """
        read FASTQ file in small batch
        yield : list of lines (quality line): batch of n size
        """
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


    def _fastq_batch_reader(self, read_batch):
        """
        parse each line of FASTQ quality line:
        return: [read length, mean Qscore, type of read (pass or fail)]
        """
        fastq_lines = []
        if self.rich:
            for read in read_batch:
                name = read[0]
                qscore = avg_qual(read[1])
                passes_filtering = True if qscore > self.threshold_Qscore else False
                if self.is_barcode:
                    start_time, ch, barcode = self._extract_info_from_name(name)
                    fastq_lines.append((len(read[1]), qscore, passes_filtering, start_time, ch, barcode))
                else:
                    start_time, ch = self._extract_info_from_name(name)
                    fastq_lines.append((len(read[1]), qscore, passes_filtering, start_time, ch))
        else:
            for read in read_batch:
                if len(read)>0:
                    qscore = avg_qual(read)
                    passes_filtering = True if qscore > self.threshold_Qscore else False
                    fastq_lines.append((len(read), qscore, passes_filtering))
        return fastq_lines


    def check_fastq(self):
        """
        """
        open_fn = gzip.open if self.fastq[0].endswith('.gz') else open

        with open_fn(self.fastq[0], 'rt') as fq:
            first_line = fq.readline().strip('\n')
        metadata = dict(x.split("=") for x in first_line.split(" ")[1:])
        if 'barcode' not in metadata:
            self.is_barcode = False
        if 'model_version_id' not in metadata:
                metadata['model_version_id'] = 'Unknow'
        try:
            return metadata['runid'] , metadata['sampleid'] , metadata['model_version_id'] 
        except:
            return None


    def _extract_info_from_name(self, name):
        """
        """
        metadata = dict(x.split("=") for x in name.split(" ")[1:])
        start_time = timeISO_to_float(metadata['start_time'], '%Y-%m-%dT%H:%M:%SZ')
        if self.is_barcode:
            return  start_time, metadata['ch'], metadata['barcode']
        return  start_time, metadata['ch']