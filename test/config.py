"""
Listing of all configurations required for unit tests
"""

whole_config= {
"barcoding": "True",
"barcode_selection" : set(["BC01, BC02, BC03, BC04, BC05, BC06, BC07"]),
"result_directory": "/home/karine/src/python/toulligqc2/data/output/",
"dpi": "100",
"sequencing_summary_source": "/home/karine/src/python/toulligqc2/data/sequencing_summary_small.txt\t/home/karine/src/python/toulligqc2/data/barcoding_summ_pass_small.txt\t/home/karine/src/python/toulligqc2/data/barcoding_summ_fail_small.txt"
}


no_seq_summary_source_config= {
"barcoding": "False",
"result_directory": "/home/karine/src/python/toulligqc2/data/output/",
"dpi": "100",
"sequencing_summary_source": ""
}

directory_config= {
"barcoding": "False",
"result_directory": "/home/karine/src/python/toulligqc2/data/output/",
"dpi": "100",
"sequencing_summary_source": "/home/karine/src/python/toulligqc2/data/"
}

only_seq_summary_config= {
"barcoding": "True",
"barcode_selection" : set(["BC01, BC02, BC03, BC04, BC05, BC06, BC07"]),
"result_directory": "/home/karine/src/python/toulligqc2/data/output/",
"dpi": "100",
"sequencing_summary_source": "/home/karine/src/python/toulligqc2/data/sequencing_summary_small.txt"
}

seq_summary_with_barcodes_config= {
"barcoding": "True",
"barcode_selection" : set(["BC01, BC02, BC03, BC04, BC05, BC06, BC07"]),
"result_directory": "/home/karine/src/python/toulligqc2/data/output/",
"dpi": "100",
"sequencing_summary_source": "/home/karine/src/python/toulligqc2/data/sequencing_summary_barcode_small.txt"
}

only_barcoding_config= {
"barcoding": "True",
"barcode_selection" : set(["BC01, BC02, BC03, BC04, BC05, BC06, BC07"]),
"result_directory": "/home/karine/src/python/toulligqc2/data/output/",
"dpi": "100",
"sequencing_summary_source": "/home/karine/src/python/toulligqc2/data/barcoding_summ_pass_small.txt\t/home/karine/src/python/toulligqc2/data/barcoding_summ_fail_small.txt"
}

random_file_config= {
"barcoding": "False",
"result_directory": "/home/karine/src/python/toulligqc2/data/output/",
"dpi": "100",
"sequencing_summary_source": "/home/karine/src/python/toulligqc2/data/test.txt"
}

