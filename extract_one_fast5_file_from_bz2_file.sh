#!/bin/bash

f=/fast5.directory/$1.tar.bz2 ; tar -xf $f --occurrence $(tar -tf $f | grep -e '\.fast5$' | head -n 1) -O > /design.file.directory/fast5_file.fast5 
