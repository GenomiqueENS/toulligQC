import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../toulligqc")
from toulligqc import sequencing_summary_extractor as sse
import unittest
from unittest.mock import patch
import pandas as pd
from pandas.util.testing import assert_frame_equal
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
    def setUpClass(self):
        
        config= {
        "barcoding": "True",
        "result_directory": "/home/karine/Bureau/ToulligQC2_B2020/data/data_small/output_debugging/",
        "dpi": "100",
        "sequencing_summary_source": "/home/karine/Bureau/ToulligQC2_B2020/data/data_small/sequencing_summary_small.txt\t/home/karine/Bureau/ToulligQC2_B2020/data/data_small/barcoding_summ_pass_small.txt\t/home/karine/Bureau/ToulligQC2_B2020/data/data_small/barcoding_summ_fail_small.txt"
        }

        self.sse_instance = sse.SequencingSummaryExtractor(config)

    #Tests for check_conf()
    def test_check_conf_whole_config(self):
        """
        Test for configuration with a sequencing_summary_file and barcoding files (pass + fail)
        """
        actual = sse.SequencingSummaryExtractor.check_conf(self.sse_instance)
        self.assertEqual((True, ""), actual)


    # Tests for _load_sequencing_summary_data()

    def test_load_sequencing_summary_data_dataframes(self):
        """
        Compare two Dataframes :
        - number of columns
        - correct names of columns
        - correct merging of dataframes
        - values returned
        """

        actual_df = sse.SequencingSummaryExtractor._load_sequencing_summary_data(
            self.sse_instance).head(11)

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

        expected_df = pd.DataFrame(data=data_test, columns=[
            'channel',
            'start_time',
            'duration',
            'num_events',
            'passes_filtering',
            'num_events_template',
            'sequence_length_template',
            'mean_qscore_template',
            'barcode_arrangement'], index=pd.Int64Index(data=pd.RangeIndex(0, 11)))

        
        expected_df = expected_df.astype({
            'channel': np.int16,
            'start_time': np.float,
            'duration': np.float,
            'num_events': np.int16,
            'num_events_template': np.int16,
            'passes_filtering': str, # set to str, because np.bool converts all string to True
            'sequence_length_template': np.int16,
            'mean_qscore_template': np.float,
            'barcode_arrangement': object
        })
        # transform strings to bool type in column passes_filtering 
        for values in expected_df['passes_filtering']:
            values = bool(util.strtobool(values))
            
        print(expected_df.dtypes)
        print(actual_df.dtypes)
        
        assert_frame_equal(expected_df, actual_df, check_names=True, check_like=True, check_dtype=False)


        


# #TODO:
# ###########################################
# # Tests with only sequencing_summary_file # 
# ###########################################

# class TestSequencingSummaryExtractorOnlySequencingSummary (unittest.TestCase):

#     """ 
#     Test SequencingSummaryExtractor class with only sequencing summary file 
#     """

#     @classmethod 
#     def setUpClass(self):

#         pass

    

# #TODO:
# ###################################
# # Tests with only barcoding files # 
# ###################################

# class TestSequencingSummaryExtractorBarcodingFiles (unittest.TestCase):

#     """ Test SequencingSummaryExtractor class with only barcoding files """

#     @classmethod 
#     def setUpClass(self):

#         pass




#########################################
# Tests with directory instead of files # 
#########################################

class TestSequencingSummaryExtractorDirectory (unittest.TestCase):

    """ Test SequencingSummaryExtractor class with a directory """

    @classmethod 
    def setUpClass(self):
        
        self.config= {
        "barcoding": "True",
        "result_directory": "/home/karine/Bureau/ToulligQC2_B2020/data/data_small/output_debugging/",
        "dpi": "100",
        "sequencing_summary_source": "/home/karine/Bureau/ToulligQC2_B2020/data/data_small/"
        }

        
    def test_init_should_raise_Value_Error(self):
        """Test if init method returns a ValueError when passing a directory"""
        
        with self.assertRaises(ValueError) as context:
            
            sse_instance = sse.SequencingSummaryExtractor(self.config)
            self.assertTrue("ValueError: The sequencing summary file must be a file path not a directory path" in context.exception)


    def test_check_conf_with_directory(self):
        """ Test if check_conf method returns a ValueError when passing a directory"""
        
        with self.assertRaises(ValueError) as context:
            
            sse_instance = sse.SequencingSummaryExtractor(self.config)
            self.assertTrue((False, "The sequencing summary file is not a file: " + "/home/karine/Bureau/ToulligQC2_B2020/data/data_small/") in context.exception)

    
#TODO:
##################################
# Tests with no sequencing files # 
##################################

class TestSequencingSummaryExtractorNoFiles(unittest.TestCase):

    """ Test SequencingSummaryExtractor class with a directory """

    @classmethod 
    def setUpClass(self):

        config= {
        "barcoding": "True",
        "result_directory": "/home/karine/Bureau/ToulligQC2_B2020/data/data_small/output_debugging/",
        "dpi": "100",
        "sequencing_summary_source": "\t" #TODO: test this config
        }
        self.sse_instance = sse.SequencingSummaryExtractor(config)
    
    def test_check_conf_no_file(self):
        """
        Test check conf method with no sequencing_summary_source
        """
        
        actual = sse.SequencingSummaryExtractor.check_conf(self.sse_instance)
        #TODO: test fails because when sequencing_summary_source is empty (ie ""), len(seq_summ_files) != 0
        #TODO: and condition "if not os.path.isfile(f):" is met
        self.assertEqual((False, 'The sequencing summary file is not a file: '), actual)


    def test_load_sequencing_summary_data_should_raise_exception(self):
        """
        Test case of no sequencing_summary_file and no barcoding one
        """

        actual = sse.SequencingSummaryExtractor._load_sequencing_summary_data(self.sse_instance)
        
        # ValueError("Sequencing summary file not found nor barcoding summary file(s)")
        #self.assertRaises(ValueError, actual)
        
        with self.assertRaises(ValueError) as context:
            actual_dataframe = sse.SequencingSummaryExtractor._load_sequencing_summary_data(sse_instance)
            
            self.assertTrue("Sequencing summary file not found nor barcoding summary file(s)" in context.exception)

