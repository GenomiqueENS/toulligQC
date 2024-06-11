FROM ubuntu:24.04

MAINTAINER Laurent Jourdren <jourdren@bio.ens.psl.eu>
ARG VERSION=2.7.0
ARG INSTALL_PACKAGES="git"
RUN apt update && \
    DEBIAN_FRONTEND=noninteractive apt install --yes \
                    $INSTALL_PACKAGES \
                    python3 \
                    python3-pip\
                    python3-matplotlib\
                    python3-plotly\
                    python3-h5py\
                    python3-pandas\
                    python3-numpy\
                    python3-scipy\
                    python3-sklearn \
                    python3-tqdm\
                    python3-pysam && \
    pip3 install --break-system-packages "pod5==0.3.10" "ezcharts==0.7.6" && \
    cd /tmp && \
    git clone https://github.com/GenomicParisCentre/toulligQC && \
    cd toulligQC && \
    git checkout v$VERSION && \
    python3 setup.py build install && \
    rm -rf /tmp/toulligQC && \
    apt remove --yes --purge $INSTALL_PACKAGES && \
    apt autoremove --yes --purge && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*
