#!/bin/bash

set -xeu

[[ "$#" -ne 1 ]] && DIR=`pwd` || DIR=$1
SDIR=$DIR/simulated_reads
CDIR=$DIR/classification-results
mkdir -p $CDIR
mkdir -p $SDIR

run_krakenu_viral() {
  FQ=$1
  NAM=$2
  K=$3
  DAT=$4

  KFILE=$CDIR/$NAM.k$K.krakenu
  [[ -s $KFILE ]] || time $DIR/install/krakenu --threads 4 --db $DIR/dbs/refseq-viral-k$K --fastq $FQ --report-file $KFILE.report > $KFILE 2> $KFILE.log
  [[ "$DAT" == "viral" ]] && SEQMAP=$DIR/dbs/refseq-viral-k$K/seqid2taxid.map || SEQMAP=$DIR/data/all-$DAT.map
  [[ -s $KFILE.results.stats ]] || $DIR/install/grade_classification  $DIR/dbs/refseq-viral-k$K/taxDB $SEQMAP $KFILE $KFILE.results > $KFILE.results.stats
}

run_kraken_viral() {
  FQ=$1
  NAM=$2
  K=$3
  DAT=$4

  KFILE=$CDIR/$NAM.k$K.kraken
  [[ -s $KFILE ]] || time kraken --threads 4 --db $DIR/dbs/refseq-viral-k$K --fastq $FQ > $KFILE 2> $KFILE.log
  [[ "$DAT" == "viral" ]] && SEQMAP=$DIR/dbs/refseq-viral-k$K/seqid2taxid.map || SEQMAP=$DIR/data/all-$DAT.map
  #[[ -s $KFILE.results.stats ]] || 
    $DIR/install/grade_classification  $DIR/dbs/refseq-viral-k$K/taxDB $SEQMAP $KFILE $KFILE.results > $KFILE.results.stats
}



AB=1m
for i in 1 2 3; do
  for dat in viral viral-neighbors bacteria archaea; do
    for len in 75 100 150; do
      NAM=$dat.$AB${len}bp.$i
      FQ=$SDIR/$NAM.fq
      [[ -f $FQ ]] || randomreads.sh -Xmx40g ref=$DIR/data/all-$dat.fna out=$FQ reads=$AB len=$len seed=$i
      for K in 21 26 31; do
        run_krakenu_viral $FQ $NAM $K $dat
        run_kraken_viral $FQ $NAM $K $dat
      done
    done
  done
done
