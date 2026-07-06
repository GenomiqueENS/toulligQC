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

import datetime
import os
import platform as pf
import sys
import tempfile as tp

from toulligqc.core.configuration import ToulligqcConf


class ToulligqcInfoExtractor:
    """Report ToulligQC run info and the list of selected extractors.

    The list of extractors called from the command line depends on the
    ``parse_args`` function (toulligqc.py) that fills the config_dictionary.

    Attributes:
        _config_dictionary: The configuration dictionary (configuration.py).
        _extractors_list: List of extractor modules that will be executed.
    """

    def __init__(self, config_dictionary: ToulligqcConf, extractors_list: list) -> None:
        self._config_dictionary = config_dictionary
        self._extractors_list = extractors_list
        self._debug = False

        if (
            "debug" in config_dictionary
            and config_dictionary["debug"].lower() == "true"
        ):
            self._debug = True

    @staticmethod
    def get_name() -> str:
        """Get the name of the extractor.

        Returns:
            The name of the extractor.
        """
        return "Toulligqc info"

    @staticmethod
    def get_report_data_file_id() -> str:
        """Get the report.data id of the extractor.

        Returns:
            The report.data id.
        """
        return "toulligqc.info.extractor"

    def check_conf(self) -> tuple[bool, str]:
        """Check the configuration.

        Returns:
            A tuple ``(is_valid, message)``; always valid for this extractor.
        """
        return True, ""

    def init(self) -> None:
        """Initialize the extractor (no-op)."""
        return

    def extract(self, result_dict: dict) -> None:
        """Extract config details and the list of extractors into result_dict.

        Args:
            result_dict: Dictionary gathering the extracted statistics.
        """

        result_dict["unwritten.keys"] = ["unwritten.keys"]

        # Add ToulligQC info
        self._toulligqc_info(result_dict)

        # Add system and python information
        if self._debug:
            self._system_and_python_info(result_dict)

        # Add QC info
        self._qc_info(result_dict)

        # Add the list of used extractors
        result_dict["toulligqc.info.extractors"] = []
        for e in self._extractors_list:
            result_dict["toulligqc.info.extractors"].append(e.get_report_data_file_id())

    def graph_generation(self, result_dict: dict) -> list:
        """Generate graphs.

        Args:
            result_dict: Dictionary gathering the extracted statistics.

        Returns:
            An empty list (this extractor produces no graph).
        """
        return []

    def clean(self, result_dict: dict) -> None:
        return

    @staticmethod
    def _system_and_python_info(result_dict: dict) -> dict:
        """Fill result_dict with the OS parameters and environment variables.

        Args:
            result_dict: Dictionary which gathers all the extracted information
                that will be reported in the report.data file.

        Returns:
            The result_dict dictionary populated with system and Python info.
        """
        result_dict["toulligqc.info.system.hostname"] = os.uname()[1]
        result_dict["toulligqc.info.system.username"] = os.environ.get("USERNAME")
        result_dict["toulligqc.info.system.user.home"] = os.environ["HOME"]
        result_dict["toulligqc.info.system.temporary.directory"] = tp.gettempdir()
        result_dict["toulligqc.info.system.operating.system"] = pf.processor()

        # Environment variables
        for name, value in os.environ.items():
            result_dict["toulligqc.info.system.env." + name] = value

        # Python info
        result_dict["toulligqc.info.python.version"] = pf.python_version()
        result_dict["toulligqc.info.python.implementation"] = pf.python_implementation()

        # Python dependencies versions
        for name, module in sorted(sys.modules.items()):
            if hasattr(module, "__version__"):
                result_dict["toulligqc.info.python.dependancy." + name + ".version"] = (
                    module.__version__
                )
            elif hasattr(module, "VERSION"):
                result_dict["toulligqc.info.python.dependancy." + name + ".version"] = (
                    module.VERSION
                )

        return result_dict

    def _toulligqc_info(self, result_dict: dict) -> None:

        result_dict["toulligqc.info.version"] = self._config_dictionary["app.version"]
        result_dict["toulligqc.info.start.time"] = (
            datetime.datetime.now().astimezone().replace(microsecond=0).isoformat()
        )
        result_dict["toulligqc.info.report.name"] = self._config_dictionary[
            "report_name"
        ]
        result_dict["toulligqc.info.executable.path"] = sys.argv[0]
        result_dict["toulligqc.info.command.line"] = sys.argv

    def _qc_info(self, result_dict: dict) -> None:

        result_dict["toulligqc.info.html.report.path"] = self._config_dictionary.get(
            "html_report_path", "Undefined"
        )
        result_dict["toulligqc.info.data.report.path"] = self._config_dictionary.get(
            "data_report_path", "Undefined"
        )
        result_dict["toulligqc.info.image.directory"] = self._config_dictionary.get(
            "images_directory", "Undefined"
        )
        result_dict["toulligqc.info.barcode.option"] = "False"

        if self._config_dictionary["barcoding"].lower() == "true":
            result_dict["toulligqc.info.barcode.option"] = "True"
            result_dict["toulligqc.info.barcode.selection"] = self._config_dictionary[
                "barcode_selection"
            ]
