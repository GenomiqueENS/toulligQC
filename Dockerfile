FROM ubuntu:16.04

MAINTAINER Lionel Ferrato Berberian <ferrato@biologie.ens.fr>

RUN apt-get update && apt-get install -y \
    software-properties-common
RUN add-apt-repository universe


RUN apt-get update && apt-get install -y \
    python3.5 \
    python3-pip\
    git\
    python3-tk -y\
    python3-h5py

RUN pip3 install matplotlib\
   seaborn

RUN apt-get install -y\
    python3-numpy\
    python3-pandas

RUN pip3 install  --upgrade pip setuptools && rm -r /root/.cache 
 
RUN mkdir /scripts
WORKDIR /scripts
RUN git clone https://github.com/lionelus/git_project.git && apt-get remove -y git
#ENTRYPOINT ["/usr/bin/python3", "/scripts/git_project/toulligqc.py"]
#CMD ["-h"]
