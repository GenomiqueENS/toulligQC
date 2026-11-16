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

# Creation of a text file containing statistics retrieved from the result_dict dictionary.
#  The information by modules is retained in a key-value form.
# The prefix of a key being the report data file id of the module.


def add_values_to_unwritten_key(result_dict: dict, values: list) -> None:
    """Register keys that must not be written to the report.data file.

    Args:
        result_dict: Dictionary gathering the extracted statistics.
        values: List of keys to add to the unwritten keys list.
    """
    return result_dict["unwritten.keys"].extend(values)


def statistics_generator(config_dictionary, result_dict: dict) -> None:
    """Write the report.data statistics file.

    Creates a text file where the information and statistics about the MinION
    run are printed in key-value form.

    Args:
        config_dictionary: Configuration dictionary holding the output path.
        result_dict: Dictionary gathering the extracted statistics.
    """

    if config_dictionary["data_report_path"] is None:
        return

    with open(config_dictionary["data_report_path"], "w") as file_data:
        for key, value in result_dict.items():
            if key not in result_dict["unwritten.keys"]:
                file_data.write(f"{key}={value}\n")
