import sys, os
#sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "../toulligqc")
import test_sequencing_summary_extractor as test_sse
import unittest

# initialization of the test loader
loader = unittest.TestLoader()

# add tests to the test loader
standard_suite = loader.loadTestsFromTestCase(test_sse.TestSequencingSummaryExtractorBarcodingFiles)
no_file_suite = unittest.TestLoader().loadTestsFromTestCase(test_sse.TestSequencingSummaryExtractorOnlySequencingSummary)
alltests = unittest.TestSuite([standard_suite, no_file_suite])

# initialize a runner, pass it your suite and run it
runner = unittest.TextTestRunner(verbosity=3).run(alltests)

# when running through command line : 
# if __name__ == '__main__':
#     unittest.main()