from pathlib import Path
# assign path to all sequencing_summary and barcode data files
path = Path(__file__).parents[1] / "test_data/sequencing_summary/"
path = str(path)

"""
Listing of all configurations required for unit tests
"""


whole_config = {
    "barcoding": "True",
    "barcode_selection": ["BC01, BC02, BC03, BC04, BC05, BC06, BC07"],
    "result_directory": path,
    "dpi": "100",
    "sequencing_summary_source": path + "/sequencing_summary_small.txt\t"\
    + path + "/barcoding_summ_pass_small.txt\t"\
    + path + "/barcoding_summ_fail_small.txt"
}


no_seq_summary_source_config= {
"barcoding": "False",
"result_directory": path,
"dpi": "100",
"sequencing_summary_source": ""
}

directory_config= {
"barcoding": "False",
"result_directory": path,
"dpi": "100",
"sequencing_summary_source": path
}

only_seq_summary_config= {
"barcoding": "False",
"barcode_selection" : set(["BC01, BC02, BC03, BC04, BC05, BC06, BC07"]),
"result_directory": path,
"dpi": "100",
"sequencing_summary_source": path + "/Albacore-2.3.1_basecall-1D-RNA_sequencing_summary.txt"
}

seq_summary_with_barcodes_config= {
"barcoding": "True",
"barcode_selection" : ["BC01, BC02, BC03, BC04, BC05, BC06, BC07"],
"result_directory": path,
"dpi": "100",
"sequencing_summary_source": path + "/Guppy-2.2.4-basecall-1D-DNA_sequencing_summary+barcode.txt"
}

only_barcoding_config= {
"barcoding": "True",
"barcode_selection" : set(["BC01, BC02, BC03, BC04, BC05, BC06, BC07"]),
"result_directory": path,
"dpi": "100",
"sequencing_summary_source": path + "/barcoding_summ_pass_small.txt\t"\
    + path + "/barcoding_summ_fail_small.txt"
}

random_file_config= {
"barcoding": "False",
"result_directory": path,
"dpi": "100",
"sequencing_summary_source": path + "/random_file.txt"
}

missing_data_config = {
"barcoding": "False",
"result_directory": path,
"dpi": "100",
"sequencing_summary_source": path + "/sequencing_summary_with_missing_data.txt"
}