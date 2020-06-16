import sys, os, re
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../toulligqc")
from toulligqc import sequencing_summary_extractor as sse
import unittest
from unittest.mock import patch, Mock, MagicMock
import config as cfg
import pandas as pd
import pandas.util.testing as testing
import numpy as np
from distutils import util

####################################################################################
# Tests of the SequencingSummaryExtractor class with several configuration cases : #
#   - Whole standard config : sequencing summary file + barcoding files            #
#   - Sequencing_summary file with barcode info embedded                           #
#   - Sequencing_summary file only (no barcodes)                                   #
#   - Barcoding files only (no sequencing_summary)                                 #
#   - Random files                                                                 #
#   - Directory case                                                               #
#   - No files                                                                     # 
####################################################################################

class TestSequencingSummaryExtractorWholeConfig (unittest.TestCase):

    """ Test SequencingSummaryExtractor class with complete configuration files """

    @classmethod 
    def setUpClass(cls):
        """
        Setup for comparing two dataframes (10 first rows):
        - number of columns
        - correct names of columns
        - correct merging of dataframes
        - values returned
        """

        cls.config = cfg.whole_config
        
        cols = [
            'channel',
            'start_time',
            'passes_filtering',
            'sequence_length_template',
            'mean_qscore_template',
            'barcode_arrangement']
        
        data_test = np.array([[20, 8.92175, True, 257, 14.379014999999999, "barcode12"],
                              [61, 9.2585, False, 388, 4.8182160000000005, "unclassified"],
                              [382, 9.2825, True, 414, 12.683063, "barcode12"],
                              [149, 8.76625, True, 850, 11.174918, "barcode12"],
                              [277, 9.00525, True, 632, 7.292675, "barcode12"],
                              [131, 9.17225, True, 852, 12.534144, "barcode08"],
                              [430, 10.3745, True, 381, 10.167657, "unclassified"],
                              [90, 8.7955, True, 1114, 12.7405, "barcode07"],
                              [123, 11.04125, True, 402, 11.276638, "unclassified"],
                              [313, 10.12025, True, 761, 13.403101000000001, "barcode12"],
                              [204, 10.84925, True, 583, 11.399572000000001, "barcode10"]])

        cls.expected_df = pd.DataFrame(data=data_test, columns=cols, index=pd.Int64Index(data=pd.RangeIndex(0, 11)))
        
        # Convert explicitly string values of passes_filtering into booleans
        cls.expected_df.passes_filtering.replace({'True': True, 'False': False}, inplace=True)
        
        cls.expected_df = cls.expected_df.astype({
            'channel': np.int16,
            'start_time': np.float,
            'passes_filtering': np.bool,
            'sequence_length_template': np.int16,
            'mean_qscore_template': np.float,
            'barcode_arrangement': object
        })


    #Tests for check_conf method
    def test_check_conf_whole_config(self):
        """
        Test for configuration with a sequencing_summary_file and barcoding files (pass + fail)
        """
        actual = sse.SequencingSummaryExtractor(self.config).check_conf()
        self.assertEqual((True, ""), actual)


    # Tests for _load_sequencing_summary_data()
    def test_load_sequencing_summary_data_columns(self):
        """
        Test if dataframe return by _load_sequencing_summary_data has 
        the expected number and names of columns 
        """

        actual_df_col = sse.SequencingSummaryExtractor(self.config)._load_sequencing_summary_data().columns

        assert self.expected_df.columns.equals(actual_df_col)
        assert len(self.expected_df.columns) == len(actual_df_col)


    def test_load_sequencing_summary_data_assert_equal(self):
        """
        Test correct merging of dataframes in _load_sequencing_summary_data
        """

        actual_df = sse.SequencingSummaryExtractor(self.config)._load_sequencing_summary_data().head(11)
        
        testing.assert_frame_equal(self.expected_df, actual_df, check_exact=True)


    def test_load_sequencing_summary_data_col_types (self):
        """
        Test column types of dataframe returned by load_sequencing_summary_data
        """

        actual_df = sse.SequencingSummaryExtractor(self.config)._load_sequencing_summary_data()
        testing.is_bool(actual_df['passes_filtering'].dtype)
        
        # test if all columns are filled with number values except passes_filtering & barcode_arrangement
        actual_df_types = actual_df.drop(['passes_filtering', 'barcode_arrangement'], axis=1, inplace=False).dtypes
        testing.is_number(actual_df_types)
        
        
    def test_load_sequencing_summary_data_barcode_data (self):
        """
        Test if values of barcode_arrangement column are correct
        """
        actual_df = sse.SequencingSummaryExtractor(self.config)._load_sequencing_summary_data()
        random_values = actual_df['barcode_arrangement'].sample(n=10)

        for index, value in random_values.iteritems():
            # when barcode_arrangement values are NaN type, do nothing
            if pd.isna(value):
                continue
            assert re.match("(^barcode)|(^unclassified)", value)

        

    #Test for init method
    def test_init_whole_config(self):
        """
        Test if instance variables created in init method are correct
        """
        # Rename expected_df columns 
        self.expected_df.rename(columns={'sequence_length_template': 'sequence_length',
                                           'mean_qscore_template': 'mean_qscore'}, inplace=True)
        
        # Create instance of SequencingSummaryExtractor class
        instance = sse.SequencingSummaryExtractor(self.config)
        # Launch init method
        instance.init()
        # Get instances variables of the dataframe_1d
        actual_df_channel = instance.channel_df.head(11)
        actual_df_passes_filtering = instance.passes_filtering_df.head(11)
        actual_df_sequence_length = instance.sequence_length_df.head(11)
        actual_df = instance.dataframe_1d.head(11)
        
        # Test all instances variables VS expected
        testing.assert_series_equal(self.expected_df['channel'], actual_df_channel)
        testing.assert_series_equal(self.expected_df['passes_filtering'], actual_df_passes_filtering)
        testing.assert_series_equal(self.expected_df['sequence_length'], actual_df_sequence_length)
        testing.assert_frame_equal(self.expected_df, actual_df)


    #Test for extract method
    def test_extract_whole_config(self):
        """Test if correct values are put in result_dict in extract method"""
        
        # Rename expected_df columns 
        self.expected_df.rename(columns={'sequence_length_template': 'sequence_length',
                                           'mean_qscore_template': 'mean_qscore'}, inplace=True)
        instance = sse.SequencingSummaryExtractor(self.config)
        instance.init()

        actual_dict = {}
        expected_dict = {}
        # create patched dict for assertions
        with patch.dict(expected_dict, {'basecaller.sequencing.summary.1d.extractor.read.fail.sorted':
            sorted(instance.dataframe_1d['start_time'].loc[instance.passes_filtering_df == bool(False)] / 3600)}):
            
            instance.extract(actual_dict)
            # compare fastq_entries value
            read_fail_sorted = 'basecaller.sequencing.summary.1d.extractor.read.fail.sorted'
            self.assertEqual(actual_dict[read_fail_sorted], expected_dict[read_fail_sorted])
            
            # compare read_pass_length_mean value
            read_pass_length_mean = 'basecaller.sequencing.summary.1d.extractor.read.pass.length.mean'
            expected_dict.update({read_pass_length_mean: pd.DataFrame.mean(
                instance.sequence_length_df[instance.passes_filtering_df == True])})
            
            assert actual_dict[read_pass_length_mean] == expected_dict[read_pass_length_mean]
            

class TestSequencingSummaryExtractorOnlySequencingSummary (unittest.TestCase):

    """ 
    Test SequencingSummaryExtractor class with only sequencing summary file 
    """

    @classmethod 
    def setUpClass(cls):

        cls.config = cfg.only_seq_summary_config
        cls.missing_data = cfg.missing_data_config
    
    
    def test_load_sequencing_summary_data_only_ss(self):
        """Test if read_id is not in the dataframe column"""
        
        actual_df = sse.SequencingSummaryExtractor(self.config)._load_sequencing_summary_data()
        assert 'read_id' not in actual_df.columns


    def test_extract_only_sequencing_summary(self):
        """
        Test if correct values are in the result_dict
        """
        
        actual_dict = {}
        
        instance = sse.SequencingSummaryExtractor(self.config)
        instance.init()
        instance.extract(actual_dict)
        
        # values to compare
        read_count = len(instance.dataframe_1d)
        read_pass_count = len(instance.dataframe_1d.loc[instance.dataframe_1d['passes_filtering'] == True])
        read_fail_count = len(instance.dataframe_1d.loc[instance.dataframe_1d['passes_filtering'] == False])
        read_pass_ratio = read_pass_count/read_count
        read_fail_ratio = read_fail_count/read_count
        read_pass_frequency = read_pass_ratio * 100
        yield_count = sum(instance.sequence_length_df)
        channel_max = pd.DataFrame.max(pd.value_counts(instance.dataframe_1d['channel']))
        mean_length = pd.DataFrame.mean(instance.sequence_length_df)
        read_pass_length_min = pd.DataFrame.min(instance.sequence_length_df[instance.dataframe_1d['passes_filtering'] == True])
        read_fail_qscore = instance.qscore_df.loc[instance.dataframe_1d['passes_filtering'] == False]
        read_fail_qscore_50 = pd.Series.quantile(read_fail_qscore, q=0.5)
        
        self.assertEqual(read_count, actual_dict['basecaller.sequencing.summary.1d.extractor.read.count'])
        self.assertEqual(read_pass_count, actual_dict['basecaller.sequencing.summary.1d.extractor.read.pass.count'])
        self.assertEqual(read_fail_count, actual_dict['basecaller.sequencing.summary.1d.extractor.read.fail.count'])
        self.assertEqual(read_pass_ratio, actual_dict['basecaller.sequencing.summary.1d.extractor.read.pass.ratio'])
        self.assertEqual(read_fail_ratio, actual_dict['basecaller.sequencing.summary.1d.extractor.read.fail.ratio'])
        self.assertEqual(read_pass_frequency, actual_dict['basecaller.sequencing.summary.1d.extractor.read.pass.frequency'])
        self.assertEqual(yield_count, actual_dict['basecaller.sequencing.summary.1d.extractor.yield'])
        self.assertEqual(channel_max, actual_dict['basecaller.sequencing.summary.1d.extractor.channel.occupancy.statistics.max'])
        self.assertEqual(mean_length, actual_dict['basecaller.sequencing.summary.1d.extractor.all.read.length.mean'])
        self.assertEqual(read_pass_length_min, actual_dict['basecaller.sequencing.summary.1d.extractor.read.pass.length.min'])
        testing.assert_series_equal(read_fail_qscore, actual_dict['basecaller.sequencing.summary.1d.extractor.read.fail.qscore'])
        self.assertEqual(read_fail_qscore_50, actual_dict['basecaller.sequencing.summary.1d.extractor.read.fail.qscore.50%'])

       
    def test_extract_barcode_info_not_called(self) :
        """Test that barcode_info method is not called for this config"""
        
        mock = MagicMock()
        mock.extract()
        mock._extract_barcode_info.assert_not_called()
        

    def test_empty_dataframe_raise_exception(self):
        """Test for raising EmptyDataError when loading missing data"""
        
        with self.assertRaises(pd.errors.EmptyDataError) as context:
            sse.SequencingSummaryExtractor(self.missing_data).init()
            self.assertTrue("Dataframe is empty", str(context))


class TestSequencingSummaryExtractorSequencingSummaryBarcodes (unittest.TestCase):

    """ 
    Test SequencingSummaryExtractor class with only sequencing summary file 
    """

    @classmethod 
    def setUpClass(cls):

        cls.config = cfg.seq_summary_with_barcodes_config
        

    def test_barcoding_methods_called(self):
        """Test that barcode methods are called"""

        mock = MagicMock()
        mock.extract()
        mock._extract_barcode_info()
        mock._get_barcode_selection
        mock._extract_barcode_info.assert_called_once()

            

class TestSequencingSummaryExtractorBarcodingFiles (unittest.TestCase):

    """ Test SequencingSummaryExtractor class with only barcoding files """

    @classmethod 
    def setUpClass(cls):

        cls.config = cfg.only_barcoding_config


    def test_check_conf_only_barcodes(self):
        """Test case of only barcoding files"""

        actual = sse.SequencingSummaryExtractor(self.config).check_conf()
        self.assertEqual((False, "No sequencing summary file has been found"), actual)


class TestSequencingSummaryExtractorDirectory (unittest.TestCase):

    """ Test SequencingSummaryExtractor class with a directory """
   
    @classmethod 
    def setUpClass(cls):
        
        cls.config = cfg.directory_config

        
    def test_check_conf_with_directory(self):
        """ Test if _is_sequencing_summary_file method returns a FileNotFoundError catched by check_conf when passing a directory"""
        
        with self.assertRaises(FileNotFoundError) as context:
            
            #Call check_conf which calls internally is_sequencing_summary_file method
            sse.SequencingSummaryExtractor(self.config).check_conf()
            sse.SequencingSummaryExtractor(self.config)._is_sequencing_summary_file(self.config['sequencing_summary_source'])
            self.assertTrue((False, "No such file or directory " + self.config.get('sequencing_summary_source')) in context.exception)
    


class TestSequencingSummaryExtractorNoFiles(unittest.TestCase):

    """ Test SequencingSummaryExtractor class with a directory """

    @classmethod 
    def setUpClass(cls):

        cls.config = cfg.no_seq_summary_source_config
        
    def test_check_conf_no_file(self):
        """
        Test check conf method with no sequencing_summary_source
        """
        
        actual = sse.SequencingSummaryExtractor(self.config).check_conf()
        self.assertEqual((False, 'No file has been defined'), actual)


    def test_load_sequencing_summary_data_should_raise_exception(self):
        """
        Test load_sequencing_summary_data method with no sequencing_summary_file and no barcoding one
        """
        
        with self.assertRaises(FileNotFoundError) as context:
            
            sse.SequencingSummaryExtractor(self.config)._load_sequencing_summary_data()
            self.assertTrue("Sequencing summary file not found", str(context))


