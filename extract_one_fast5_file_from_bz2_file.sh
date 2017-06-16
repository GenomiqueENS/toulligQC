#!/bin/bash

f=/home/ferrato/shares-net/sequencages/nanopore/fast5/raw/$1.tar.bz2 ; tar -xf $f --occurrence $(tar -tf $f | grep -e '\.fast5$' | head -n 1) -O >/home/ferrato/ownCloud/toulligQC/fast5_file.fast5 
