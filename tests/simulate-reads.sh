#!/bin/bash

set -xeu

[[ "$#" -ne 1 ]] && DIR=`pwd` || DIR=$1
SDIR=$DIR/simulated_reads
mkdir -p $SDIR

randomreads.sh ref=$DIR/data/all-viral-neighbors.fna out=$SDIR/viral-neighbors-10m.fq reads=10m len=150
