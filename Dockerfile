FROM ubuntu:23.10

MAINTAINER Laurent Jourdren <jourdren@bio.ens.psl.eu>
ARG VERSION=2.6
RUN apt update && \
    DEBIAN_FRONTEND=noninteractive apt install --yes \
                    python3 \
                    python3-pip\
                    git\
                    python3-tk\
                    python3-h5py\
                    python3-matplotlib\
                    python3-scipy\
                    python3-pandas\
                    python3-numpy\
                    python3-tqdm\
                    python3-pysam\
                    python3-sklearn\
                    python3-plotly && \
    pip3 install --break-system-packages "pod5==0.3.6" && \
    cd /tmp && \
    git clone https://github.com/GenomicParisCentre/toulligQC && \
    cd toulligQC && \
    git checkout v$VERSION && \
    python3 setup.py build install && \
    apt remove --yes git && \
    apt autoremove --yes && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*
