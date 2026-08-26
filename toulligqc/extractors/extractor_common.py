#
#                  ToulligQC development code
#
# This code may be freely distributed and modified under the
# terms of the GNU General Public License version 3 or later
# and CeCILL. This should be distributed with the code. If you
# do not have a copy, see:
#
#      http://www.gnu.org/licenses/gpl-3.0-standalone.html
#      http://www.cecill.info/licences/Licence_CeCILL_V2-en.html
#
# Copyright for this code is held jointly by the Genomic platform
# of the Institut de Biologie de l'École Normale Supérieure and
# the individual authors.
#
# For more information on the ToulligQC project and its aims,
# visit the home page at:
#
#      https://github.com/GenomiqueENS/toulligQC
#

# This module contains common methods for sequencing summary modules.

import bz2
import gzip
import sys
import time
from datetime import datetime

import pandas as pd

from toulligqc.core import common


def set_result_value(extractor, result_dict: dict, key: str, value) -> None:
    """Set a key/value pair in the result_dict.

    Args:
        extractor: Extractor whose report data file id namespaces the key.
        result_dict: Dictionary gathering the extracted statistics.
        key: String entry to add to result_dict.
        value: int, float, list, pd.Series or pd.DataFrame value of the
            corresponding key.
    """

    _check_result_key_value(key, value)
    result_dict[extractor.get_report_data_file_id() + "." + key] = value


def get_result_value(extractor, result_dict: dict, key: str):
    """Return the value associated with a result_dict key.

    Args:
        extractor: Extractor whose report data file id namespaces the key.
        result_dict: Dictionary gathering the extracted statistics.
        key: String entry to look up in result_dict.

    Returns:
        The value associated with the namespaced key.

    Raises:
        KeyError: If the key is not present in result_dict.
    """
    if (extractor.get_report_data_file_id() + "." + key) not in result_dict.keys():
        raise KeyError("Key {key} not found").__format__(key)
    return result_dict.get(extractor.get_report_data_file_id() + "." + key)


def check_result_values(extractor, result_dict: dict) -> None:
    prefix = extractor.get_report_data_file_id() + "."
    for key, value in result_dict.items():
        if key.startswith(prefix):
            _check_result_key_value(key, value)


def _check_result_key_value(key, value) -> None:
    if not isinstance(key, str):
        raise TypeError(f"Invalid type for key: {type(key)}")

    if (
        not isinstance(value, int)
        and not isinstance(value, float)
        and not isinstance(value, str)
    ):
        raise TypeError(f"Invalid type for the value of the key {key}: {type(value)} ")


def describe_dict(extractor, result_dict: dict, function, entry: str) -> None:
    """Set describe statistics for a key.

    Fills the result_dict with the mean, min, max, median and percentiles
    (without the count value) computed from the given series.

    Args:
        extractor: Extractor whose report data file id namespaces the key.
        result_dict: Dictionary gathering the extracted statistics.
        function: Series returning the values to describe.
        entry: Entry to put in result_dict completed with the statistics.
    """
    stats = pd.Series.describe(function).drop("count")
    for key, value in stats.items():
        set_result_value(extractor, result_dict, entry + "." + key, value)


def count_boolean_elements(
    dataframe: pd.DataFrame, column_name: str, boolean: bool
) -> int:
    """Return the number of values of a column filtered by a boolean.

    Args:
        dataframe: DataFrame to filter.
        column_name: Name of the dataframe column.
        boolean: Boolean value to filter on.

    Returns:
        The number of rows whose column value equals the boolean.
    """
    return len(dataframe.loc[dataframe[column_name] == bool(boolean)])


def series_cols_boolean_elements(
    dataframe: pd.DataFrame, column_name1: str, column_name2: str, boolean: bool
) -> pd.Series:
    """Return the values of a column filtered by a boolean on another column.

    Args:
        dataframe: The dataframe_1d to filter.
        column_name1: Column whose values are returned.
        column_name2: Column used for the boolean filter.
        boolean: Boolean value used to filter the rows.

    Returns:
        A pd.Series object with the selected values.
    """
    return dataframe[column_name1].loc[dataframe[column_name2] == bool(boolean)]


def df_cols_boolean_elements(
    dataframe: pd.DataFrame, column_name1: list, column_name2: str, boolean: bool
) -> pd.Series:
    """Return the values of columns filtered by a boolean on another column.

    Args:
        dataframe: The dataframe_1d to filter.
        column_name1: Columns whose values are returned.
        column_name2: Column used for the boolean filter.
        boolean: Boolean value used to filter the rows.

    Returns:
        A pd.Series object with the selected values.
    """
    return dataframe[column_name1].loc[dataframe[column_name2] == bool(boolean)]


def sorted_series_boolean_elements_divided(
    dataframe: pd.DataFrame,
    column_name1: str,
    column_name2: str,
    boolean: bool,
    denominator: int,
) -> pd.Series:
    """Return a sorted series of column values filtered by a boolean and divided.

    Args:
        dataframe: The dataframe_1d to filter.
        column_name1: Column whose values are returned.
        column_name2: Column used for the boolean filter.
        boolean: Boolean value used to filter the rows.
        denominator: Number to divide the values by.

    Returns:
        A sorted pd.Series of the filtered values divided by the denominator.
    """
    return (
        dataframe[column_name1].loc[dataframe[column_name2] == bool(boolean)]
        / denominator
    ).sort_values()


def fill_series_dict(df_dict: dict, df: pd.DataFrame) -> None:
    """Fill a dict with the read length, qscore, time and duration series.

    Args:
        df_dict: Dictionary populated with the per-read-type series.
        df: DataFrame holding the sequencing summary columns.
    """
    for read_type in ["pass", "fail"]:
        read_type_bool = True if read_type == "pass" else False

        # Read length series
        df_dict[read_type + ".reads.sequence.length"] = series_cols_boolean_elements(
            df, "sequence_length", "passes_filtering", read_type_bool
        )

        # Read qscore series
        df_dict[read_type + ".reads.mean.qscore"] = series_cols_boolean_elements(
            df, "mean_qscore", "passes_filtering", read_type_bool
        )

    # Read length series
    df_dict["all.reads.sequence.length"] = df["sequence_length"]

    # Mean QScore
    df_dict["all.reads.mean.qscore"] = df["mean_qscore"]

    # Time series
    if "start_time" in df:
        df_dict["all.reads.start.time"] = df["start_time"]

    # Duration series
    if "duration" in df:
        df_dict["all.reads.duration"] = df["duration"]


def extract_barcode_info(
    extractor,
    result_dict: dict,
    barcode_selection: list,
    dataframe_dict: dict,
    df: pd.DataFrame,
) -> None:
    """Gather all barcode info for graphs: pass/fail reads and per-barcode frequency.

    Args:
        extractor: Extractor whose report data file id namespaces the keys.
        result_dict: Dictionary gathering the extracted statistics.
        barcode_selection: List of selected barcodes.
        dataframe_dict: Dictionary populated with the barcode dataframes.
        df: DataFrame holding the sequencing summary columns.
    """
    # Add values unclassified and other to barcode list
    if "unclassified" not in barcode_selection:
        barcode_selection.append("unclassified")

    # If the barcode_arrangement column contains a barcode kit id
    mask = df["barcode_arrangement"].str.startswith(("SQK", "VQK"))

    if mask.any():
        df["barcode_arrangement"] = df["barcode_arrangement"].astype(str)
        df.loc[mask, "barcode_arrangement"] = df.loc[
            mask, "barcode_arrangement"
        ].str.extract(r"[SV]QK-.+_(.+)$")[0]

    # Create keys barcode.arrangement, and read.pass/fail.barcode in dataframe_dict with all values of
    # column barcode_arrangement when reads are passed/failed
    dataframe_dict["barcode.arrangement"] = df["barcode_arrangement"]

    # Print warning message if a barcode is unknown
    barcodes_found = set(df["barcode_arrangement"].unique())
    for element in barcode_selection:
        if element not in barcodes_found and element != "other barcodes":
            sys.stderr.write(
                f"\033[93mWarning:\033[0m The barcode {element} doesn't exist in input data\n"
            )

    # Get barcodes frequency by Bases
    df_base_pass_barcode = series_cols_boolean_elements(
        df, ["barcode_arrangement", "sequence_length"], "passes_filtering", True
    )

    dataframe_dict["base.pass.barcoded"] = _barcode_bases(
        extractor,
        barcode_selection,
        result_dict,
        "base.pass.barcoded",
        df_base_pass_barcode,
    )

    df_base_fail_barcode = series_cols_boolean_elements(
        df, ["barcode_arrangement", "sequence_length"], "passes_filtering", False
    )

    dataframe_dict["base.fail.barcoded"] = _barcode_bases(
        extractor,
        barcode_selection,
        result_dict,
        "base.fail.barcoded",
        df_base_fail_barcode,
    )

    # Get barcodes frequency by read type
    series_read_pass_barcode = series_cols_boolean_elements(
        df, "barcode_arrangement", "passes_filtering", True
    )

    dataframe_dict["read.pass.barcoded"] = _barcode_frequency(
        extractor,
        barcode_selection,
        result_dict,
        "read.pass.barcoded",
        series_read_pass_barcode,
    )

    series_read_fail_barcode = series_cols_boolean_elements(
        df, "barcode_arrangement", "passes_filtering", False
    )

    dataframe_dict["read.fail.barcoded"] = _barcode_frequency(
        extractor,
        barcode_selection,
        result_dict,
        "read.fail.barcoded",
        series_read_fail_barcode,
    )

    read_pass_barcoded_count = get_result_value(
        extractor, result_dict, "read.pass.barcoded.count"
    )
    read_fail_barcoded_count = get_result_value(
        extractor, result_dict, "read.fail.barcoded.count"
    )

    # Add key "read.pass.barcoded.frequency"
    total_reads = get_result_value(extractor, result_dict, "read.count")
    set_result_value(
        extractor,
        result_dict,
        "read.pass.barcoded.frequency",
        (read_pass_barcoded_count / total_reads) * 100,
    )

    # Add key "read.fail.barcoded.frequency"
    set_result_value(
        extractor,
        result_dict,
        "read.fail.barcoded.frequency",
        (read_fail_barcoded_count / total_reads) * 100,
    )

    # Replaces all rows with unused barcodes (ie not in barcode_selection) in column barcode_arrangement with the 'other' value

    df.loc[
        ~df["barcode_arrangement"].isin(barcode_selection), "barcode_arrangement"
    ] = "other barcodes"

    if "other barcodes" not in barcode_selection:
        barcode_selection.append("other barcodes")

    # Create dataframes filtered by barcodes and read quality
    for index_barcode, barcode in enumerate(barcode_selection):
        barcode_all_reads_df = df[df["barcode_arrangement"] == barcode]
        barcode_pass_reads_df = barcode_all_reads_df.loc[df["passes_filtering"]]
        barcode_fail_reads_df = barcode_all_reads_df.loc[~df["passes_filtering"]]

        # Add all barcode statistics to result_dict based on values of selected dataframes
        _barcode_stats(
            extractor,
            result_dict,
            barcode_all_reads_df,
            barcode_pass_reads_df,
            barcode_fail_reads_df,
            barcode,
        )

    # Add filtered dataframes (all info by barcode and by length or qscore) to dataframe_dict
    _barcode_selection_dataframe(
        dataframe_dict,
        df,
        "sequence_length",
        "barcode_selection_sequence_length_dataframe",
        "length",
    )
    _barcode_selection_dataframe(
        dataframe_dict,
        df,
        "mean_qscore",
        "barcode_selection_sequence_phred_dataframe",
        "qscore",
    )


def _barcode_selection_dataframe(
    dataframe_dict: dict,
    df: pd.DataFrame,
    df_column_name: str,
    df_key_name: str,
    melted_column_name: str,
) -> None:
    """Build a per-barcode, per-read-type dataframe for length or qscore boxplots.

    Groups all reads by barcode and read type (pass/fail) for read length or
    Phred score info, then reshapes the dataframe from wide to long format so it
    displays barcode, read type and read length or Phred score per read.

    Args:
        dataframe_dict: Dictionary populated with the resulting dataframe.
        df: The dataframe_1d holding the reads.
        df_column_name: Name of the dataframe_1d column used for the new
            barcode selection dataframe.
        df_key_name: String name of the entry to put in dataframe_dict.
        melted_column_name: Value (qscore or length) used for renaming the
            melted dataframe column.
    """
    # Count total number of rows
    nrows = df.shape[0]
    # Create a new dataframe with 3 columns : 'passes_filtering', 'barcode_arrangement' and the column name parameter
    filtered_df = df.filter(
        items=["passes_filtering", df_column_name, "barcode_arrangement"]
    )

    # Reshape dataframe with new MultiIndex : numbered index of df length + passes filtering index and then shape data by barcode
    barcode_selection_dataframe = filtered_df.set_index(
        [pd.RangeIndex(start=0, stop=nrows), "passes_filtering"], drop=True
    ).pivot(columns="barcode_arrangement")

    # Remove the column parameter index
    barcode_selection_dataframe.columns.droplevel(level=0)

    # Remove sequence_length Multindex to only have barcode_arrangement column labels
    barcode_selection_dataframe.columns = barcode_selection_dataframe.columns.droplevel(
        level=0
    )

    # Reset index to have all labels in the same level
    barcode_selection_dataframe.reset_index(level="passes_filtering", inplace=True)

    # Add final dataframe to dataframe_dict
    dataframe_dict[df_key_name] = barcode_selection_dataframe


def _barcode_stats(
    extractor,
    result_dict: dict,
    barcode_selected_dataframe: pd.DataFrame,
    barcode_selected_read_pass_dataframe: pd.DataFrame,
    barcode_selected_read_fail_dataframe: pd.DataFrame,
    barcode_name: str,
) -> None:
    """Put per-barcode length and qscore statistics in result_dict.

    Uses the describe method to compute statistics for each selected dataframe
    (all.read, read.pass and read.fail). The count statistic is not included for
    qscore.

    Args:
        extractor: Extractor whose report data file id namespaces the keys.
        result_dict: Dictionary gathering the extracted statistics.
        barcode_selected_dataframe: All reads for the barcode.
        barcode_selected_read_pass_dataframe: Passing reads for the barcode.
        barcode_selected_read_fail_dataframe: Failing reads for the barcode.
        barcode_name: Name of the barcode.
    """
    df_dict = {
        "all.read.": barcode_selected_dataframe,
        "read.pass.": barcode_selected_read_pass_dataframe,
        "read.fail.": barcode_selected_read_fail_dataframe,
    }

    for df_name, df in df_dict.items():  # df_dict.items = all.read/read.pass/read.fail
        for stats_index, stats_value in df["sequence_length"].describe().items():
            key_to_result_dict = (
                df_name + barcode_name.replace(" ", ".") + ".length." + stats_index
            )
            set_result_value(extractor, result_dict, key_to_result_dict, stats_value)

        for stats_index, stats_value in (
            df["mean_qscore"].describe().drop("count").items()
        ):
            key_to_result_dict = df_name + barcode_name + ".qscore." + stats_index
            set_result_value(extractor, result_dict, key_to_result_dict, stats_value)


def _barcode_frequency(
    extractor,
    barcode_selection: list,
    result_dict: dict,
    entry: str,
    df_filtered: pd.Series,
) -> pd.Series:
    """Count reads per barcode and compute their frequencies.

    Counts reads by values of barcode_selection, computes the sum of counts by
    barcode_selection and the sum of unclassified counts, regroups all non-used
    barcodes under the index "other" and computes all frequency values for each
    number of barcoded reads.

    Args:
        extractor: Extractor whose report data file id namespaces the keys.
        barcode_selection: List of selected barcodes.
        result_dict: Dictionary gathering the extracted statistics.
        entry: Entry about barcoded counts.
        df_filtered: Series of barcode arrangements to count.

    Returns:
        A pd.Series with all barcode (used, non-used and unclassified) counts.
    """
    # Regroup all barcoded read in Series
    all_barcode_count = df_filtered.value_counts()

    # Retain only existing barcodes from barcode_selection list
    barcodes_found = set(df_filtered.unique())
    barcode_selection_existing = [x for x in barcode_selection if x in barcodes_found]

    # Sort by list of barcode_selection
    count_sorted = all_barcode_count.sort_index()[barcode_selection_existing]
    # Replace all NaN values to zero
    count_sorted.fillna(0, inplace=True)

    # Compute sum of all used barcodes without barcode 'unclassified'
    # set_result_value(extractor, result_dict, entry + '.count', sum(count_sorted.drop("unclassified")))
    if "unclassified" in count_sorted.index:
        set_result_value(
            extractor,
            result_dict,
            entry + ".count",
            sum(count_sorted.drop("unclassified")),
        )
    else:
        set_result_value(extractor, result_dict, entry + ".count", sum(count_sorted))

    # Replace entry name ie read.pass/fail.barcode with read.pass/fail.non.used.barcodes.count
    non_used_barcodes_count_key = entry.replace(".barcoded", ".non.used.barcodes.count")

    # Compute all reads of barcodes that are not in the barcode_selection list
    other_barcode_count = sum(all_barcode_count) - sum(count_sorted)
    set_result_value(
        extractor, result_dict, non_used_barcodes_count_key, other_barcode_count
    )

    # Create Series for all non-used barcode counts and rename index array with "other"
    other_all_barcode_count = pd.Series(other_barcode_count, index=["other barcodes"])

    # Append Series of non-used barcode counts to the Series of barcode_selection counts
    count_sorted = pd.concat([count_sorted, other_all_barcode_count]).sort_index()

    # Compute frequency for all barcode counts and save into dataframe_dict
    for barcode in count_sorted.to_dict():
        frequency_value = count_sorted[barcode] * 100 / sum(count_sorted)
        set_result_value(
            extractor,
            result_dict,
            entry.replace(".barcoded", ".") + barcode + ".frequency",
            frequency_value,
        )

    return count_sorted


def _barcode_bases(
    extractor,
    barcode_selection: list,
    result_dict: dict,
    entry: str,
    df_filtered: pd.DataFrame,
) -> pd.Series:
    """Count bases per barcode and compute their frequencies.

    Counts bases by values of barcode_selection, computes the sum of counts by
    barcode_selection and the sum of unclassified counts, regroups all non-used
    barcodes under the index "other" and computes all frequency values for each
    number of barcoded bases.

    Args:
        extractor: Extractor whose report data file id namespaces the keys.
        barcode_selection: List of selected barcodes.
        result_dict: Dictionary gathering the extracted statistics.
        entry: Entry about barcoded counts.
        df_filtered: DataFrame with ``barcode_arrangement`` and
            ``sequence_length`` columns.

    Returns:
        A pd.Series with all barcode (used, non-used and unclassified) base counts.
    """
    # Regroup all barcoded and sum all read lengths in df
    all_barcode_count = df_filtered.groupby("barcode_arrangement", observed=False)[
        "sequence_length"
    ].sum()

    # Retain only existing barcodes from barcode_selection list
    barcodes_found = set(df_filtered["barcode_arrangement"].unique())
    barcode_selection_existing = [x for x in barcode_selection if x in barcodes_found]

    # Sort by list of barcode_selection
    count_sorted = all_barcode_count.sort_index()[barcode_selection_existing]
    # Replace all NaN values to zero
    count_sorted.fillna(0, inplace=True)

    # Compute sum of all used barcodes without barcode 'unclassified'
    # set_result_value(extractor, result_dict, entry + '.count', sum(count_sorted.drop("unclassified")))
    if "unclassified" in count_sorted.index:
        set_result_value(
            extractor,
            result_dict,
            entry + ".count",
            sum(count_sorted.drop("unclassified")),
        )
    else:
        set_result_value(extractor, result_dict, entry + ".count", sum(count_sorted))
    # Replace entry name ie read.pass/fail.barcode with read.pass/fail.non.used.barcodes.count
    non_used_barcodes_count_key = entry.replace(".barcoded", ".non.used.barcodes.count")

    # Compute all reads of barcodes that are not in the barcode_selection list
    other_barcode_count = sum(all_barcode_count) - sum(count_sorted)
    set_result_value(
        extractor, result_dict, non_used_barcodes_count_key, other_barcode_count
    )

    # Create Series for all non-used barcode counts and rename index array with "other"
    other_all_barcode_count = pd.Series(other_barcode_count, index=["other barcodes"])

    # Append Series of non-used barcode counts to the Series of barcode_selection counts
    count_sorted = pd.concat([count_sorted, other_all_barcode_count]).sort_index()

    # Compute frequency for all barcode counts and save into dataframe_dict
    for barcode in count_sorted.to_dict():
        frequency_value = count_sorted[barcode] * 100 / sum(count_sorted)
        set_result_value(
            extractor,
            result_dict,
            entry.replace(".barcoded", ".") + barcode + ".frequency",
            frequency_value,
        )

    return count_sorted


def log_task(quiet: bool, msg: str, start_time: float, end_time: float) -> None:
    if not quiet:
        delta = end_time - start_time
        print(f"  - {msg} in {common.format_duration(delta)}")


def add_image_to_result(
    quiet: bool, image_list: list, start_time: float, image: tuple
) -> None:
    end_time = time.time()
    log_task(quiet, f'Creation of image "{image[0]}"', start_time, end_time)
    image_list.append(image)


def timeISO_to_float(iso_datetime: str, format: str) -> float:
    """Convert an ISO datetime string to a Unix timestamp.

    Args:
        iso_datetime: Datetime string to parse.
        format: Primary datetime format to try before falling back to
            ``%Y-%m-%dT%H:%M:%SZ``.

    Returns:
        The corresponding Unix timestamp as a float.
    """
    try:
        dt = datetime.strptime(iso_datetime, format)
    except ValueError:
        format = "%Y-%m-%dT%H:%M:%SZ"
        dt = datetime.strptime(iso_datetime, format)
    unix_timestamp = dt.timestamp()
    return unix_timestamp


def read_first_line_file(filename: str) -> str:
    """Load the first line of a file.

    Supports plain, gzip- and bzip2-compressed files.

    Args:
        filename: Name of the file to load.

    Returns:
        The first line of the file.

    Raises:
        FileNotFoundError: If the file cannot be opened.
    """
    try:
        if filename.endswith(".gz"):
            with gzip.open(filename, "rt") as f:
                return f.readline()
        elif filename.endswith(".bz2"):
            with bz2.open(filename, "rt") as f:
                return f.readline()
        else:
            with open(filename) as f:
                return f.readline()
    except OSError:
        raise FileNotFoundError


def set_result_dict_telemetry_value(result_dict: dict, key: str, new_value) -> None:
    """Set a telemetry value in result_dict, keeping any existing value.

    Args:
        result_dict: Dictionary gathering the extracted statistics.
        key: Telemetry key (without the extractor prefix).
        new_value: New value to store; if None, the current value is kept.
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


def pd_read_sequencing_summary(file: str, cols: list, data_type: dict) -> pd.DataFrame:
    try:
        return pd.read_csv(file, sep="\t", usecols=cols, dtype=data_type)
    except (ValueError, KeyError):
        # Handle case where 'passes_filtering' column doesn't exist
        del data_type["passes_filtering"]
        cols.remove("passes_filtering")
        return pd.read_csv(file, sep="\t", usecols=cols, dtype=data_type)
