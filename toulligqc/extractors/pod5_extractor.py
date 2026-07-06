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

# Extraction of the information about the Pod5 files

import os
import shutil
import sys
import tarfile
import tempfile

import pod5 as p5

from toulligqc.core.common import find_file_in_directory, set_result_dict_value
from toulligqc.core.configuration import ToulligqcConf


class Pod5Extractor:
    """Extract run information from a Pod5 file.

    Attributes:
        pod5_source: Pod5 file or directory.
        report_name: Report name.
        pod5_file_extension: Extension used for the storage of the set of Pod5
            files, if any.
    """

    def __init__(self, config_dictionary: ToulligqcConf) -> None:
        self.config_file_dictionary = config_dictionary
        self.pod5_source = config_dictionary["pod5_source"]
        self.file_to_process = None
        self.report_name = config_dictionary["report_name"]
        self.pod5_file_extension = ""
        self.pod5_file = ""
        self.get_report_data_file_id()

    def check_conf(self) -> tuple[bool, str]:
        """Check the configuration.

        Returns:
            A tuple ``(is_valid, message)`` where the message describes the
            error when the configuration is invalid.
        """

        if not os.path.exists(self.pod5_source):
            return (
                False,
                "The input file or directory for Pod5 file does not exists: "
                + self.pod5_source,
            )

        if os.path.isdir(self.pod5_source):
            file_found = find_file_in_directory(self.pod5_source, "pod5")
            if file_found is None:
                return False, "No Pod5 file found in directory: " + self.pod5_source
            self.file_to_process = file_found
        else:
            self.file_to_process = self.pod5_source

        if self.file_to_process.endswith(".tar"):
            self.pod5_file_extension = "tar"

        elif self.file_to_process.endswith(".tar.gz"):
            self.pod5_file_extension = "tar.gz"

        elif self.file_to_process.endswith(".tar.bz2"):
            self.pod5_file_extension = "tar.bz2"

        elif self.file_to_process.endswith(".pod5"):
            self.pod5_file_extension = "pod5"

        else:
            return (
                False,
                "The file extension for Pod5 input is not supported "
                "(only .pod5, .tar, .tar.gz or .tar.bz2 are supported): "
                + self.pod5_source,
            )

        return True, ""

    def init(self) -> None:
        """Initialize the extractor (no-op for Pod5)."""
        return

    @staticmethod
    def get_name() -> str:
        """Get the name of the extractor.

        Returns:
            The name of the extractor.
        """
        return "Pod5"

    @staticmethod
    def get_report_data_file_id() -> str:
        """Get the report.data id of the extractor.

        Returns:
            The report.data id.
        """
        return "pod5.extractor"

    def extract(self, result_dict: dict) -> None:
        """Extract the different information from the Pod5 files.

        Args:
            result_dict: Dictionary which gathers all the extracted information
                that will be reported in the report.data file.
        """
        p5_file = self._read_pod5()
        run_info_dict = self._get_pod5_items(p5_file)
        tracking_id_dict = run_info_dict.tracking_id
        if len(tracking_id_dict) == 0:
            return

        prefix = "sequencing.telemetry.extractor"
        result_dict[prefix + ".source"] = self.pod5_source
        set_result_dict_value(
            result_dict, prefix + ".flowcell.id", tracking_id_dict, "flow_cell_id"
        )
        set_result_dict_value(
            result_dict, prefix + ".minknow.version", tracking_id_dict, "version"
        )
        set_result_dict_value(
            result_dict, prefix + ".hostname", tracking_id_dict, "hostname"
        )
        set_result_dict_value(
            result_dict,
            prefix + ".operating.system",
            tracking_id_dict,
            "operating_system",
        )
        set_result_dict_value(
            result_dict, prefix + ".run.id", tracking_id_dict, "run_id"
        )
        set_result_dict_value(
            result_dict,
            prefix + ".protocol.run.id",
            tracking_id_dict,
            "protocol_run_id",
        )
        set_result_dict_value(
            result_dict,
            prefix + ".protocol.group.id",
            tracking_id_dict,
            "protocol_group_id",
        )
        set_result_dict_value(
            result_dict, prefix + ".sample.id", tracking_id_dict, "sample_id"
        )
        set_result_dict_value(
            result_dict, prefix + ".exp.start.time", tracking_id_dict, "exp_start_time"
        )
        set_result_dict_value(
            result_dict, prefix + ".device.id", tracking_id_dict, "device_id"
        )
        set_result_dict_value(
            result_dict, prefix + ".device.type", tracking_id_dict, "device_type"
        )
        set_result_dict_value(
            result_dict,
            prefix + ".distribution.version",
            tracking_id_dict,
            "distribution_version",
        )
        set_result_dict_value(
            result_dict,
            prefix + ".flow.cell.product.code",
            tracking_id_dict,
            "flow_cell_product_code",
        )

        context_tags_dict = run_info_dict.context_tags
        if len(context_tags_dict) != 0:
            set_result_dict_value(
                result_dict,
                prefix + ".selected.speed.bases.per.second",
                context_tags_dict,
                "selected_speed_bases_per_second",
            )
            set_result_dict_value(
                result_dict,
                prefix + ".sample.frequency",
                context_tags_dict,
                "sample_frequency",
            )
            set_result_dict_value(
                result_dict,
                prefix + ".sequencing.kit.version",
                context_tags_dict,
                "sequencing_kit",
            )

    def graph_generation(self, result_dict: dict) -> list:
        """Generate graphs.

        Args:
            result_dict: Dictionary gathering the extracted statistics.

        Returns:
            An empty list (this extractor produces no graph).
        """
        return []

    def clean(self, result_dict: dict) -> None:
        """Delete the temporary Pod5 file extracted from the tar archive.

        Also removes dictionary entries that will not be kept in the
        report.data file.

        Args:
            result_dict: Dictionary which gathers all the extracted information
                that will be reported in the report.data file.
        """
        if self.temporary_directory:
            shutil.rmtree(self.temporary_directory, ignore_errors=True)

    def _read_pod5(self) -> p5.Reader:
        """Extract one Pod5 file and open it as a p5.Reader object.

        Returns:
            The opened Pod5 file as a p5.Reader object.
        """
        self.temporary_directory = tempfile.mkdtemp()
        if (
            self.pod5_file_extension == "tar"
            or self.pod5_file_extension == "tar.gz"
            or self.pod5_file_extension == "tar.bz2"
        ):
            self.pod5_file = self._pod5_tar_extraction(
                self.file_to_process, self.pod5_file_extension, self.temporary_directory
            )
        elif self.pod5_file_extension == "pod5" or self.pod5_file_extension == ".pod5":
            self.pod5_file = self.file_to_process
        else:
            err_msg = "There is a problem with the pod5 file or the tar file"
            sys.exit(err_msg)
        p5_file = p5.Reader(self.pod5_file)

        return p5_file

    def _pod5_tar_extraction(
        self, tar_file: str, extension: str, output_directory: str
    ) -> str:
        """Extract a Pod5 file stored in a tar archive.

        Args:
            tar_file: Tar file containing the set of raw Pod5 files.
            extension: Extension of the tar file (``tar``, ``tar.gz`` or
                ``tar.bz2``).
            output_directory: Directory where the Pod5 file is extracted.

        Returns:
            The path to the extracted Pod5 file.
        """

        if extension == "tar.gz":
            tar_mode = "r:gz"
        elif extension == "tar.bz2":
            tar_mode = "r:bz2"
        else:
            tar_mode = "r"

        tf = tarfile.open(name=tar_file, mode=tar_mode)
        while True:
            member = tf.next()
            if member.name.endswith(".pod5"):
                tf.extract(member, path=output_directory)
                break
        return output_directory + "/" + member.name

    def _get_pod5_items(self, h5py_file):
        """Extract run information from the first read of a Pod5 file.

        Args:
            h5py_file: Opened Pod5 file (p5.Reader object).

        Returns:
            The run info of the first read, or an empty dict if the file has
            no reads.
        """

        for read_record in h5py_file.reads():
            return read_record.run_info
        return {}
