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

# Extraction of run information from the sequencing_telemetry.js file

import bz2
import gzip
import json
import os.path

from toulligqc.configuration import ToulligqcConf


class SequencingTelemetryExtractor:
    def __init__(self, config_dictionary: ToulligqcConf) -> None:
        self.config_file_dictionary = config_dictionary
        self.telemetry_source = config_dictionary["sequencing_telemetry_source"]

        if os.path.isdir(self.telemetry_source):
            self.telemetry_file = self.telemetry_source + "/sequencing_telemetry.js"
        else:
            self.telemetry_file = self.telemetry_source

    def check_conf(self) -> tuple[bool, str]:
        """Check the configuration.

        Returns:
            A tuple ``(is_valid, message)`` where the message describes the
            error when the configuration is invalid.
        """

        if not os.path.isfile(self.telemetry_file):
            return False, "Telemetry file does not exists: " + self.telemetry_file

        return True, ""

    def init(self) -> None:
        """Initialize the extractor (no-op for telemetry)."""
        return

    @staticmethod
    def get_name() -> str:
        """Get the name of the extractor.

        Returns:
            The name of the extractor.
        """
        return "Sequencing telemetry"

    @staticmethod
    def get_report_data_file_id() -> str:
        """Get the report.data id of the extractor.

        Returns:
            The report.data id.
        """
        return "sequencing.telemetry.extractor"

    def graph_generation(self, result_dict: dict) -> list:
        """Generate graphs.

        Args:
            result_dict: Dictionary gathering the extracted statistics.

        Returns:
            An empty list (this extractor produces no graph).
        """
        return []

    def clean(self, result_dict: dict) -> None:
        """Clean up the extractor (no-op for telemetry).

        Args:
            result_dict: Dictionary gathering the extracted statistics.
        """
        return

    def extract(self, result_dict: dict) -> None:
        """Extract data from the sequencing_telemetry.js file.

        Args:
            result_dict: Dictionary which gathers all the extracted information
                that will be reported in the report.data file.
        """

        array = _load_json(self.telemetry_file)

        result_dict[self.get_report_data_file_id() + ".source"] = self.telemetry_file

        self._set_result_dict_value(
            result_dict, ".flowcell.id", array, "tracking_id", "flow_cell_id"
        )
        self._set_result_dict_value(
            result_dict, ".minknow.version", array, "tracking_id", "version"
        )
        self._set_result_dict_value(
            result_dict, ".hostname", array, "tracking_id", "hostname"
        )
        self._set_result_dict_value(
            result_dict, ".operating.system", array, "tracking_id", "operating_system"
        )
        self._set_result_dict_value(
            result_dict, ".run.id", array, "tracking_id", "run_id"
        )
        self._set_result_dict_value(
            result_dict, ".protocol.run.id", array, "tracking_id", "protocol_run_id"
        )
        self._set_result_dict_value(
            result_dict, ".protocol.group.id", array, "tracking_id", "protocol_group_id"
        )
        self._set_result_dict_value(
            result_dict, ".sample.id", array, "tracking_id", "sample_id"
        )
        self._set_result_dict_value(
            result_dict, ".exp.start.time", array, "tracking_id", "exp_start_time"
        )
        self._set_result_dict_value(
            result_dict, ".device.id", array, "tracking_id", "device_id"
        )
        self._set_result_dict_value(
            result_dict, ".device.type", array, "tracking_id", "device_type"
        )
        self._set_result_dict_value(
            result_dict,
            ".distribution.version",
            array,
            "tracking_id",
            "distribution_version",
        )
        self._set_result_dict_value(
            result_dict,
            ".flow.cell.product.code",
            array,
            "tracking_id",
            "flow_cell_product_code",
        )
        self._set_result_dict_value(
            result_dict, ".basecalling.date", array, "tracking_id", "time_stamp"
        )

        self._set_result_dict_value(
            result_dict, ".software.name", array, "software", "name"
        )
        self._set_result_dict_value(
            result_dict, ".software.version", array, "software", "version"
        )
        self._set_result_dict_value(
            result_dict, ".software.analysis", array, "software", "analysis"
        )

        if "albacore_opts" in array[0]:
            self._set_result_dict_value(
                result_dict, ".kit.version", array, "albacore_opts", "kit"
            )
            self._set_result_dict_value(
                result_dict, ".flowcell.version", array, "albacore_opts", "flowcell"
            )
            self._set_result_dict_value(
                result_dict,
                ".model.file",
                array,
                "albacore_opts",
                "local_bc_temp_model",
            )

        if "opts" in array[0]:
            self._set_result_dict_value(
                result_dict, ".kit.version", array, "opts", "kit"
            )
            self._set_result_dict_value(
                result_dict,
                ".sequencing.kit.version",
                array,
                "context_tags",
                "sequencing_kit",
            )
            self._set_result_dict_value(
                result_dict, ".barcode.kits.version", array, "opts", "barcode_kits"
            )
            self._set_result_dict_value(
                result_dict, ".flowcell.version", array, "opts", "flowcell"
            )
            self._set_result_dict_value(
                result_dict, ".model.file", array, "opts", "model_file"
            )
            self._set_result_dict_value(
                result_dict, ".pass.threshold.qscore", array, "opts", "min_qscore"
            )
            self._set_result_dict_value(
                result_dict,
                ".selected.speed.bases.per.second",
                array,
                "context_tags",
                "selected_speed_bases_per_second",
            )
            self._set_result_dict_value(
                result_dict,
                ".sample.frequency",
                array,
                "context_tags",
                "sample_frequency",
            )

    def _set_result_dict_value(
        self, result_dict: dict, key: str, array, dict_name: str, dict_key: str
    ) -> None:

        final_key = self.get_report_data_file_id() + key
        current_value = None
        new_value = None

        if final_key in result_dict:
            current_value = result_dict[final_key]
            if len(current_value) == 0:
                current_value = None

        if dict_name in array[0] and dict_key in array[0][dict_name]:
            new_value = array[0][dict_name][dict_key]

        if new_value is None:
            new_value = current_value

        if new_value is None:
            new_value = ""

        result_dict[final_key] = new_value


def _load_json(filename: str):
    """Load a JSON file, handling gzip- and bzip2-compressed files.

    Args:
        filename: Name of the file to load.

    Returns:
        The parsed JSON object.
    """
    if filename.endswith(".gz"):
        print("Load json gzip")
        with gzip.open(filename, "rb") as f:
            return json.loads(f.read(), encoding="utf-8")
    elif filename.endswith(".bz2"):
        print("Load json bzip2")
        with bz2.open(filename, "rb") as f:
            return json.loads(f.read(), encoding="utf-8")
    else:
        with open(filename) as f:
            return json.load(f)
