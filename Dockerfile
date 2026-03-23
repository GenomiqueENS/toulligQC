FROM ubuntu:24.04

ARG VERSION=2.8.4
ARG UV_VERSION=0.9.17
ARG INSTALL_PACKAGES="curl git"
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
                    python3-pysam && \
    pip3 install --break-system-packages "pod5==0.3.10" "ezcharts==0.15.2" && \
    curl --proto '=https' --tlsv1.2 -LsSf https://github.com/astral-sh/uv/releases/download/${UV_VERSION}/uv-installer.sh | sh && \
    cd /tmp && \
    git clone https://github.com/GenomiqueENS/toulligQC.git && \
    cd toulligQC && \
    git checkout v$VERSION && \
    cd /tmp/toulligQC && \
    /root/.local/bin/uv build && \
    cd dist && \
    pip3 install --break-system-packages toulligqc-*.tar.gz && \
    cd / && \
    rm -rf /tmp/toulligQC && \
    rm -rf /root/.local && \
    mkdir -p /root/.local/bin && \
    touch /root/.local/bin/env && \
    apt remove --yes --purge $INSTALL_PACKAGES && \
    apt autoremove --yes --purge && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*
