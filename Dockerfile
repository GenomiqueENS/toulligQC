FROM ubuntu:16.04

RUN apt-get update && apt-get install -y \
    software-properties-common
RUN add-apt-repository universe


RUN apt-get update
RUN apt-get install -y \
    python3.5 \
    python3-pip\
    git

RUN apt-get install python3-tk -y
RUN apt-get install python3-pip -y 

RUN pip3 install --upgrade pip setuptools && \ 
rm -r /root/.cache

RUN pip3 install matplotlib
RUN pip3 install h5py
RUN pip3 install pandas
RUN pip3 install seaborn
RUN pip3 install numpy
RUN pip3 install PyPDF2
RUN pip3 install python-docx

RUN mkdir /scripts
WORKDIR /scripts
RUN git clone https://github.com/GenomicParisCentre/toulligQC.git


