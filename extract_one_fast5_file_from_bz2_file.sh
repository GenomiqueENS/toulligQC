#!/bin/bash

f=$1.tar.bz2 ; tar -xf $f --occurrence $(tar -tf $f | grep -e '\.fast5$' | head -n 1) -O >$2 
