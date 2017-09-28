#!/bin/bash
set -xeu

[[ "$#" -ne 1 ]] && DIR=`pwd` || DIR=$1
SDIR=$DIR/simulated_reads
CDIR=$DIR/classification-results
mkdir -p $CDIR

NAM=viral-neighbors-10m
time $DIR/install/krakenu --threads 4 --db $DIR/dbs/refseq-viral --fastq ~/kraken-hll-test/simulated_reads/$NAM.fq --report-file $CDIR/$NAM.krakenu.report > $CDIR/$NAM.krakenu
