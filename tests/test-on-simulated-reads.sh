#!/bin/bash

set -eu

[[ "$#" -ne 1 ]] && DIR=`pwd` || DIR=$1
SDIR=$DIR/simulated_reads
CDIR=$DIR/classification-results
CCDIR=$DIR/classification-stats
mkdir -p $CDIR
mkdir -p $CCDIR
mkdir -p $SDIR

[[ `uname` == "Darwin" ]] && THREADS=4 || THREADS=10

run_kraken() {
  local FQ=$1
  local NAM=$2
  local DAT=$3
  local DB_DAT=$4
  local DB_K=$5
  local PROG=$6
  local ALWAYS_SEQMAP=$7;
  local DB=refseq-$DB_DAT-k$K

  mkdir -p $CDIR/against-$DB
  mkdir -p $CCDIR/against-$DB
  local KFILE=$CDIR/against-$DB/$NAM.against-$DB.$PROG
  local KKFILE=$CCDIR/against-$DB/$NAM.against-$DB.$PROG

  if [[ "$PROG" == "kraken" ]]; then 
    CMD="kraken"
  elif [[ "$PROG" == "krakenuniq" ]]; then
    CMD="$DIR/install/krakenuniq --report-file $KFILE.report"
  elif [[ "$PROG" == "krakenhull" ]]; then
    CMD="$DIR/install/krakenuniq --report-file $KFILE.report --uid-mapping"
  else 
    echo "Unknown $PROG"
    return;
  fi

  if [[ ! -s $KFILE ]]; then 
    echo "$CMD --threads $THREADS --db $DIR/dbs/$DB --fastq $FQ --output $KFILE"
    time $CMD --threads $THREADS --db $DIR/dbs/$DB --fastq $FQ --output $KFILE 2>&1 | tee $KFILE.log
  fi
  
  [[ "$DAT" == "$DB_DAT" ]] && SEQMAP=$DIR/dbs/$DB/seqid2taxid.map || SEQMAP=$DIR/data/all-$DAT.map
  [[ "$ALWAYS_SEQMAP" == "ALWAYS_SEQMAP" ]] && SEQMAP=$DIR/dbs/$DB/seqid2taxid.map

  if [[ ! -s "$KKFILE.results.stats" ]]; then
    $DIR/install/grade_classification  $DIR/dbs/$DB/taxDB $SEQMAP $KFILE $KKFILE.results > $KKFILE.results.stats
  else
    echo "$KKFILE.results.stats exist"
  fi
}

AB=1m
for i in 1; do # 2 3
  for dat in viral viral-neighbors bacteria archaea; do
    for len in 100; do ## 75 150
      NAM=$dat.$AB${len}bp.$i
      FQ=$SDIR/$NAM.fq
      [[ -f $FQ ]] || randomreads.sh -Xmx40g ref=$DIR/data/all-$dat.fna out=$FQ reads=$AB len=$len seed=$i
      for K in 31; do
        # run_kraken $FQ $NAM $dat viral $K krakenuniqid
        if [[ `uname` != "Darwin" ]]; then
          run_kraken $FQ $NAM $dat oct2017 $K kraken ALWAYS_SEQMAP
          run_kraken $FQ $NAM $dat oct2017 $K krakenuniq ALWAYS_SEQMAP
          run_kraken $FQ $NAM $dat oct2017 $K krakenuniqid ALWAYS_SEQMAP
        else
          run_kraken $FQ $NAM $dat viral $K kraken
          run_kraken $FQ $NAM $dat viral $K krakenuniq
        fi
      done
    done
  done
done
