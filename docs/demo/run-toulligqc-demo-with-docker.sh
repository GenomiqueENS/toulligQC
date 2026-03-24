#!/bin/bash

DATA_URL="http://outils.genomique.biologie.ens.fr/leburon/downloads/toulligqc-example/toulligqc_demo_data.tar.bz2"

# Check if required commands are installed
for c in $(echo curl tar docker); do
    command -v $c > /dev/null
    if [ "$?" -ne 0 ]; then
        echo "Error: $c command is required in the PATH to use this script." >&2
        exit 1
    fi
done

# Get script directory
WORKING_DIR=$(readlink -f ".")

# Download and untar data if necessary
if [ ! -d "$WORKING_DIR/toulligqc_demo_data" ]; then

    if [ ! -f "$WORKING_DIR/toulligqc_demo_data.tar.bz2" ]; then
        echo "* Download demo data"
        curl --silent --remote-name --location "$DATA_URL"
	if [ "$?" -ne 0 ]; then
            echo "Error: Unable to download data" >&2
            exit 1
	fi
    fi

    echo "* Untar demo data"
    tar -xjf "$WORKING_DIR/toulligqc_demo_data.tar.bz2"
    if [ "$?" -ne 0 ]; then
        echo "Error: Unable to untar data" >&2
        exit 1
    fi

    # Remove old demo scripts from demo data archive
    chmod +w "$WORKING_DIR/toulligqc_demo_data/run-toulligqc.sh" "$WORKING_DIR/toulligqc_demo_data/run-toulligqc-with-docker.sh" 2> /dev/null
    rm "$WORKING_DIR/toulligqc_demo_data/run-toulligqc.sh" "$WORKING_DIR/toulligqc_demo_data/run-toulligqc-with-docker.sh" 2> /dev/null
fi

# Get the full data directory path
DATA_DIR="$WORKING_DIR/toulligqc_demo_data"

# Create output directory
if [ ! -d "$WORKING_DIR/output" ]; then
    mkdir "$WORKING_DIR/output"
    if [ $? -ne 0 ]; then
        echo "Error: Unable to create output directory" >&2
        exit 1
    fi
fi

# Create temporary directory
if [ ! -d "$WORKING_DIR/tmp" ]; then
    mkdir "$WORKING_DIR/tmp"
    if [ $? -ne 0 ]; then
        echo "Error: Unable to create temporary directory" >&2
        exit 1
    fi
fi

# Check ToulligQC version
TOULLIGQC_VERSION=$(docker run genomicpariscentre/toulligqc:latest toulligqc --version)
if [ "$?" -ne 0 ]; then
    echo "Error: Unable to get ToulligQC version" >&2
    exit 1
fi
if [ "${TOULLIGQC_VERSION:0:2}" != "2." ]; then
    echo "Error: Invalid ToulligQC version (ToulligQC >= 2.0 is required): $TOULLIGQC_VERSION" >&2
    exit 1
fi

echo "* Launch ToulligQC"
docker run -ti \
           --rm \
           -v "$DATA_DIR:$DATA_DIR" \
           -v "$WORKING_DIR/output:$WORKING_DIR/output" \
           -v "$WORKING_DIR/tmp:/tmp" \
           -e TMPDIR=/tmp \
           --read-only \
	   -u $(id -u):$(id -g) \
           genomicpariscentre/toulligqc:latest \
	   toulligqc \
             --report-name               ToulligQC_Demo_Data \
             --barcoding \
	     --barcodes                  BC01,BC02,BC03,BC04,BC05,BC07 \
             --telemetry-source          "$DATA_DIR/sequencing_telemetry.js" \
             --sequencing-summary-source "$DATA_DIR/sequencing_summary.txt" \
             --sequencing-summary-source "$DATA_DIR/barcoding_summary_pass.txt" \
	     --sequencing-summary-source "$DATA_DIR/barcoding_summary_fail.txt" \
	     --output-directory          "$WORKING_DIR/output"
