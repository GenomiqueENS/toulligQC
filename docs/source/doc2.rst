===============
getter1D
===============

This module provided informations about the minion runs and the fastq sequences

*get_MinknowVersion(h5py_file)*
   parameter: fast5 file open with h5py
   Get the Minknow version from fast5 file

*getFlowcellId(h5py_file)*
   parameter: fast5 file open with h5py
   Get the flowcell id from fast5 file

*get_Hostname(h5py_file)*
   parameter: fast5 file open with h5py
   Get the hostname from fast5 file

*getNumMinION(h5py_file)*
   parameter: fast5 file open with h5py
   Get the number of Minion run

*getProtocolRunId(h5py_file)*
   parameter: fast5 file open with h5py
   Get the run id protocol from fast5 file

*get_barcode()*
   parameter: fast5 file open with h5py
   Get the barcode from a file given in input

*get_fastq(selection)*
   parameter: fast5 file open with h5py
   Get the fastq sequence
