#!/bin/bash
set -xeu

[[ "$#" -ne 1 ]] && DIR=`pwd` || DIR=$1
SDIR=$DIR/simulated_reads
CDIR=$DIR/classification-results
mkdir -p $CDIR

NAM=viral-neighbors-10m
for K in 21 26 31; do
  KFILE=$CDIR/$NAM.k$K.krakenu
  [[ -s $KFILE ]] || time $DIR/install/krakenu --threads 4 --db $DIR/dbs/refseq-viral-k$K --fastq ~/kraken-hll-test/simulated_reads/$NAM.fq --report-file $KFILE.report > $KFILE 2> $KFILE.log
  [[ -s $KFILE.results ]] || $DIR/install/grade_classification  $DIR/dbs/refseq-viral-k$K/taxDB $DIR/data/all-viral-neighbors.map $KFILE > $KFILE.results
  [[ -s $KFILE.results.stats ]] || cut -f 4 $KFILE.results | sort | uniq -c | sort -n > $KFILE.results.stats

done
