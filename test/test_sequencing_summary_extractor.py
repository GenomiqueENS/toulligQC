import sys, os, re
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../toulligqc")
from toulligqc import sequencing_summary_extractor as sse
import unittest
import config as cfg
import pandas as pd
import pandas.util.testing as testing
import numpy as np
from distutils import util

#TODO: put some tests files and include them in test/ instead of using local config tests files
#TODO: my_data_path = os.path.join(THIS_DIR, os.pardir, 'data_folder/data.csv')
#TODO: my_data_path = os.path.join(THIS_DIR, 'testdata.csv')

###########################################################################
# Tests with whole configuration (sequencing_summary and barcoding files) # 
###########################################################################

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
        
        data_test = np.array([[20, 8.92175, 1.22075, 2441, True, 1564, 257, 14.379014999999999, "barcode12"],
                              [61, 9.2585, 1.31325, 2626, False, 2447,
                                  388, 4.8182160000000005, "unclassified"],
                              [382, 9.2825, 1.15525, 2310, True,
                                  2152, 414, 12.683063, "barcode12"],
                              [149, 8.76625, 2.0585, 4117, True,
                               4041, 850, 11.174918, "barcode12"],
                              [277, 9.00525, 1.851, 3702, True,
                                  2952, 632, 7.292675, "barcode12"],
                              [131, 9.17225, 2.43525, 4870, True,
                               4318, 852, 12.534144, "barcode08"],
                              [430, 10.3745, 1.03325, 2066, True,
                               1933, 381, 10.167657, "unclassified"],
                              [90, 8.7955, 2.91875, 5837, True,
                                  5707, 1114, 12.7405, "barcode07"],
                              [123, 11.04125, 1.28325, 2566, True,
                               2419, 402, 11.276638, "unclassified"],
                              [313, 10.12025, 1.989, 3978, True, 3861,
                               761, 13.403101000000001, "barcode12"],
                              [204, 10.84925, 1.63625, 3272, True, 2816, 583, 11.399572000000001, "barcode10"]])

        cls.expected_df = pd.DataFrame(data=data_test, columns=[
            'channel',
            'start_time',
            'duration',
            'num_events',
            'passes_filtering',
            'num_events_template',
            'sequence_length_template',
            'mean_qscore_template',
            'barcode_arrangement'], index=pd.Int64Index(data=pd.RangeIndex(0, 11)))
        
        # Convert explicitly string values of passes_filtering into booleans
        cls.expected_df.passes_filtering.replace({'True': True, 'False': False}, inplace=True)
        
        cls.expected_df = cls.expected_df.astype({
            'channel': np.int16,
            'start_time': np.float,
            'duration': np.float,
            'num_events': np.int16,
            'num_events_template': np.int16,
            'passes_filtering': np.bool,
            'sequence_length_template': np.int16,
            'mean_qscore_template': np.float,
            'barcode_arrangement': object
        })


    #Tests for check_conf()
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
            assert re.match("(^barcode)|(^unclassified)", value)

        

###########################################
# Tests with only sequencing_summary_file # 
###########################################

class TestSequencingSummaryExtractorOnlySequencingSummary (unittest.TestCase):

    """ 
    Test SequencingSummaryExtractor class with only sequencing summary file 
    """

    @classmethod 
    def setUpClass(cls):

        cls.config = cfg.only_seq_summary_config
    

    

# #TODO:
###################################
# Tests with only barcoding files # 
###################################

class TestSequencingSummaryExtractorBarcodingFiles (unittest.TestCase):

    """ Test SequencingSummaryExtractor class with only barcoding files """

    @classmethod 
    def setUpClass(cls):

        pass




#########################################
# Tests with directory instead of files # 
#########################################

class TestSequencingSummaryExtractorDirectory (unittest.TestCase):

    """ Test SequencingSummaryExtractor class with a directory """
   
    @classmethod 
    def setUpClass(cls):
        
        cls.config = cfg.directory_config

        
    def test_init_should_raise_Value_Error(self):
        """Test if init method returns a ValueError when passing a directory"""

        with self.assertRaises(IsADirectoryError) as context:
            sse.SequencingSummaryExtractor(self.config).__init__()
            self.assertTrue("ValueError: The sequencing summary file must be a file path not a directory path" in context.exception)

    #TODO: change this method, bc it tests only __init__ which will be refactored
    def test_check_conf_with_directory(self):
        """ Test if check_conf method returns a ValueError when passing a directory"""
        
        with self.assertRaises(IsADirectoryError) as context:
            
            sse.SequencingSummaryExtractor(self.config).check_conf()
            self.assertTrue((False, "The path is a directory : " + self.config.get('sequencing_summary_source')) in context.exception)
    


##################################
# Tests with no sequencing files # 
##################################

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


