from math import log
import os
import numpy as np
import pandas as pd
import time
import pysam
from collections import defaultdict
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
from toulligqc.extractor_common import extract_barcode_info
from toulligqc.common_statistics import compute_NXX, compute_LXX, occupancy_channel, avg_qual
from toulligqc.fastq_bam_common import multiprocessing_submit, extract_headerTag
from toulligqc.fastq_bam_common import batch_iterator
from toulligqc.common import is_numpy_1_24
from toulligqc import plotly_graph_generator as pgg


class uBAM_Extractor:
    def __init__(self, config_dictionary):
        self.config_dictionary = config_dictionary
        self.ubam = config_dictionary['bam'].split('\t')
        self.images_directory = config_dictionary['images_directory']
        self.threshold_Qscore = int(config_dictionary['threshold'])
        self.batch_size = int(config_dictionary['batch_size'])
        self.thread = int(config_dictionary['thread'])
        self.header = dict()
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
        self.dataframe = self._load_uBAM_file()
        if self.dataframe.empty:
            raise pd.errors.EmptyDataError("Dataframe is empty")
        self.dataframe_dict = {}
        
        # Add missing categories
        if 'barcode_arrangement' in self.dataframe.columns:
            self.dataframe['barcode_arrangement'] = self.dataframe['barcode_arrangement'].cat.add_categories([0,
                                                                                        'other barcodes',
                                                                                        'passes_filtering'])
        
        # Replace all NaN values by 0 to avoid data manipulation errors when columns are not the same length
        self.dataframe = self.dataframe.fillna(0)
            
        self.barcode_selection = self.config_dictionary['barcode_selection']

        log_task(self.quiet,
                 'Load BAM file ({:,.2f} MB used)'.format(self.dataframe.memory_usage(deep=True).sum()/1024/1024),
                 start_time,
                 time.time())


    def clean(self, result_dict):
        """
        Removing dictionary entries that will not be kept in the report.data file
        :return:
        """
        check_result_values(self, result_dict)
        self.dataframe_dict.clear()
        self.dataframe.iloc[0:0]


    @staticmethod
    def get_name() -> str:
        """
        Get the name of the extractor.
        :return: the name of the extractor
        """
        return 'uBAM'


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
        add_image_to_result(self.quiet, images, time.time(), pgg.yield_plot(self.dataframe, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.read_quality_multiboxplot(self.dataframe_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.allphred_score_frequency(self.dataframe_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.plot_performance(self.dataframe, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.twod_density(self.dataframe_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.sequence_length_over_time(self.dataframe_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.phred_score_over_time(self.dataframe_dict, result_dict, self.images_directory))
        add_image_to_result(self.quiet, images, time.time(), pgg.speed_over_time(self.dataframe_dict, self.images_directory))
        if self.is_barcode:
            if "barcode_alias" in self.config_dictionary:
                barcode_alias = self.config_dictionary['barcode_alias']
            else:
                barcode_alias = None 
            add_image_to_result(self.quiet, images, time.time(), pgg.barcode_percentage_pie_chart_pass(self.dataframe_dict,
                                                                                                       self.barcode_selection,
                                                                                                       self.images_directory,
                                                                                                       barcode_alias))

            read_fail = self.dataframe_dict["read.fail.barcoded"]
            if not (len(read_fail) == 1 and read_fail["other barcodes"] == 0):
                add_image_to_result(self.quiet, images, time.time(), pgg.barcode_percentage_pie_chart_fail(self.dataframe_dict,
                                                                                                      self.barcode_selection,
                                                                                                      self.images_directory,
                                                                                                      barcode_alias))

            add_image_to_result(self.quiet, images, time.time(), pgg.barcode_length_boxplot(self.dataframe_dict,
                                                                                            self.images_directory,
                                                                                            barcode_alias))

            add_image_to_result(self.quiet, images, time.time(), pgg.barcoded_phred_score_frequency(self.dataframe_dict,
                                                                                                    self.images_directory,
                                                                                                    barcode_alias))
        return images


    def extract(self, result_dict):
        """
        Get Phred score (Qscore) and Length details (frequencies, ratios, yield and statistics) per type read (pass or fail)
        :param result_dict:
        """
        start_time = time.time()
        fill_series_dict(self.dataframe_dict, self.dataframe)

        set_result_dict_telemetry_value(result_dict, "run.id", self.header["run_id"])
        set_result_dict_telemetry_value(result_dict, "sample.id", self.header["sample_id"])
        set_result_dict_telemetry_value(result_dict, "model.file", self.header["model_version_id"])
        set_result_dict_telemetry_value(result_dict, "software.name", self.header["basecaller"])
        set_result_dict_telemetry_value(result_dict, "software.version", self.header["basecaller_version"])
        set_result_dict_telemetry_value(result_dict, "flowcell.id", self.header["flow_cell_id"])
        set_result_dict_telemetry_value(result_dict, "basecalling.date", self.header["run_date"])
        set_result_dict_telemetry_value(result_dict, "pass.threshold.qscore", str(self.threshold_Qscore))

        set_result_value(self, result_dict, "read.count", len(self.dataframe))
        set_result_value(self, result_dict, "read.pass.count",
                         count_boolean_elements(self.dataframe, 'passes_filtering', True))
        set_result_value(self, result_dict, "read.fail.count",
                         count_boolean_elements(self.dataframe, 'passes_filtering', False))
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

        set_result_value(self, result_dict, "run.time", max(self.dataframe['start_time']))
        # Get channel occupancy statistics and store each value into result_dict
        for index, value in occupancy_channel(self.dataframe).items():
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
        qscore_statistics = self.dataframe['mean_qscore'].describe().drop(
            "count")

        for index, value in qscore_statistics.items():
            set_result_value(self,
                             result_dict, "all.read.qscore." + index, value)

        # Add statistics (without count) about read pass/fail qscore in the result_dict
        describe_dict(self, result_dict, self.dataframe_dict["pass.reads.mean.qscore"], "pass.reads.mean.qscore")
        describe_dict(self, result_dict, self.dataframe_dict["fail.reads.mean.qscore"], "fail.reads.mean.qscore")
        if self.is_barcode:
            extract_barcode_info(self, result_dict,
                                 self.barcode_selection,
                                 self.dataframe_dict,
                                 self.dataframe)
                                 
        log_task(self.quiet, 'Extract info from uBAM file', start_time, time.time())       


    def _load_uBAM_file(self):
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
        uBAM_df = []
        
        for _, f in enumerate(rst_futures):
            uBAM_df.extend(f.result())

        columns = ['sequence_length', 'mean_qscore', 'passes_filtering', 'start_time', 'channel', 'duration']
        if self.is_barcode:
            columns.append('barcode_arrangement')

        uBAM_df = pd.DataFrame(uBAM_df, columns=columns)

        uBAM_df['sequence_length'] = uBAM_df['sequence_length'].astype(np.uint32)
        uBAM_df['mean_qscore'] = uBAM_df['mean_qscore'].astype(np.float32)
        uBAM_df['passes_filtering'] = uBAM_df['passes_filtering'].astype(np.bool_ if is_numpy_1_24 else np.bool)
        uBAM_df["start_time"] = uBAM_df["start_time"] - uBAM_df["start_time"].min()
        uBAM_df['channel'] = uBAM_df['channel'].astype(np.int16)
        uBAM_df['start_time'] = uBAM_df['start_time'].astype(np.float64)
        uBAM_df['duration'] = uBAM_df['duration'].astype(np.float32)
        if self.is_barcode:
            uBAM_df['barcode_arrangement'] = uBAM_df['barcode_arrangement'].astype("category")
        return uBAM_df 


    def _uBAM_batch_reader(self, uBAM_chunk):
        """
        parse each line of uBAM quality line:
        return: [read length, mean Qscore, type of read (pass or fail)]
        """
        #def process_bam_chunk(bam_chunk):
        rec_data = []
        record_count = 0
        for rec in uBAM_chunk:
            record_count += 1
            rec_dict = self._process_record(rec, record_count)
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


    def _get_header(self):
        sam_file = pysam.AlignmentFile(self.ubam[0], "rb", check_sq=False)
        header = sam_file.header.to_dict()
        run_id, model_version_id = extract_headerTag(header, 'RG','ID',
                                                     'Unknown_Unknown').split('_', 1)
        self.header = {
            "run_id": run_id,
            "run_date": extract_headerTag(header, 'RG', 'DT', 'Unknown'),
            "sample_id": extract_headerTag(header, 'RG', 'SM', 'Unknown'),
            "basecaller": extract_headerTag(header, 'PG', 'PN', 'Unknown'),
            "basecaller_version": extract_headerTag(header, 'PG', 'VN', 'Unknown'),
            "model_version_id": model_version_id,
            "flow_cell_id": extract_headerTag(header, 'RG', 'PU', 'Unknown')
        }


    def _process_record(self, rec, record_count):
        """
        extract QC info from BAM record
        return : dict of QC info
        """
        fields = rec.split("\t")

        # Parse optional fields
        attributes = {}
        for t in fields[11:]:
            k, t, v = t.split(':', 2)
            attributes[k] = v

        iso_start_time = attributes.get('st', None)
        qual = avg_qual(fields[10])
        passes_filtering = True if qual > self.threshold_Qscore else False
        data = [
            len(fields[9]), # read length
            qual, # AVG Qscore
            passes_filtering, # Passing filter
            float(record_count) if iso_start_time is None else timeISO_to_float(iso_start_time, '%Y-%m-%dT%H:%M:%S.%f%z'), # start time
            attributes.get('ch', '1'),  # Channel
            attributes.get('du', '1')  # Duration
        ]
        if self.is_barcode:
            data.append(attributes.get('BC', 'unclassified'))
        return data
