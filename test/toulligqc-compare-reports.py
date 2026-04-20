#!/usr/bin/python3

import sys

import pytest

default_places = 5


def parse_value(value_str):
    """Parse a values and return int, float or string."""
    try:
        return int(value_str)
    except ValueError:
        try:
            return float(value_str)
        except ValueError:
            return value_str


def compare_floats(a, b, places=default_places):

    return a == pytest.approx(b, nan_ok=True, rel=1e-4)


def is_skip_key(key):

    if key in (
        "toulligqc.info.version",
        "toulligqc.info.start.time",
        "toulligqc.info.report.name",
        "toulligqc.info.executable.path",
        "toulligqc.info.command.line",
        "toulligqc.info.image.directory",
    ):
        return True

    if key.endswith(".duration"):
        return True

    if key.endswith(".report.path"):
        return True

    return False


def compare_files(file1_path, file2_path):
    """Compare keys/values of two files and returns differences."""
    with open(file1_path) as f1, open(file2_path) as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

    dict1 = {}
    dict2 = {}

    # Parse first file
    for line in lines1:
        line = line.strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        dict1[key.strip()] = parse_value(value.strip())

    # Parse second file
    for line in lines2:
        line = line.strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        dict2[key.strip()] = parse_value(value.strip())

    # Find differences
    differences = []

    # Keys in the first file but not in the second file
    for key in dict1:
        if key not in dict2:
            differences.append((key, dict1[key], None, "Missing key in second file"))

    # Keys in the second file but not in the first file
    for key in dict2:
        if key not in dict1:
            differences.append((key, None, dict2[key], "Missing key in first file"))

    # Keys in both files but not equals
    for key in dict1:
        if is_skip_key(key):
            continue

        if key in dict2:
            val1 = dict1[key]
            val2 = dict2[key]

            if isinstance(val1, float) and isinstance(val2, float):
                if not compare_floats(val1, val2):
                    differences.append((key, val1, val2, "Float values are not equal"))
            elif val1 != val2:
                differences.append((key, val1, val2, "Values are not equal"))

    return differences


def main():

    if len(sys.argv) != 3:
        sys.exit("Syntax: compare-toulligqc-report.py file1 file2")

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]

    differences = compare_files(file1_path, file2_path)

    if len(differences) > 0:
        print("Differences found:")
        for diff in differences:
            key, val1, val2, message = diff
            print(f"Clé: {key}")
            print(f"  File 1: {val1}")
            print(f"  File 2: {val2}")
            print(f"  Message: {message}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
