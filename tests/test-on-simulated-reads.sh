#!/bin/bash

set -eu

[[ "$#" -ne 1 ]] && DIR=`pwd` || DIR=$1
SDIR=$DIR/simulated_reads
CDIR=$DIR/classification-results
mkdir -p $CDIR
mkdir -p $SDIR

[[ `uname` == "Darwin" ]] && THREADS=4 || THREADS=10

run_kraken() {
  FQ=$1
  NAM=$2
  DAT=$3
  DB_DAT=$4
  DB_K=$5
  PROG=$6
  DB=refseq-$DB_DAT-k$K
  mkdir -p $CDIR/against-$DB
  KFILE=$CDIR/against-$DB/$NAM.against-$DB.$PROG

  if [[ "$PROG" == "kraken" ]]; then 
    CMD="kraken"
  elif [[ "$PROG" == "krakenu" ]]; then
    CMD="$DIR/install/krakenu --report-file $KFILE.report"
  elif [[ "$PROG" == "krakenuid" ]]; then
    CMD="$DIR/install/krakenu --report-file $KFILE.report --uid-mapping"
  else 
    echo "Unknown $PROG"
    return;
  fi

  if [[ ! -s $KFILE ]]; then 
    echo "$CMD --threads $THREADS --db $DIR/dbs/$DB --fastq $FQ --output $KFILE"
    time $CMD --threads $THREADS --db $DIR/dbs/$DB --fastq $FQ --output $KFILE 2>&1 | tee $KFILE.log
  fi
  #[[ "$DAT" == "$DB_DAT" ]] && SEQMAP=$DIR/dbs/$DB/seqid2taxid.map || SEQMAP=$DIR/data/all-$DAT.map
  #[[ -s $KFILE.results.stats ]] || $DIR/install/grade_classification  $DIR/dbs/$DB/taxDB $SEQMAP $KFILE $KFILE.results > $KFILE.results.stats
}

AB=1m
for i in 1 2 3; do
  for dat in viral viral-neighbors bacteria archaea; do
    for len in 75 100 150; do
      NAM=$dat.$AB${len}bp.$i
      FQ=$SDIR/$NAM.fq
      [[ -f $FQ ]] || randomreads.sh -Xmx40g ref=$DIR/data/all-$dat.fna out=$FQ reads=$AB len=$len seed=$i
      for K in 31; do
        run_kraken $FQ $NAM $dat viral $K kraken
        run_kraken $FQ $NAM $dat viral $K krakenu
        #run_kraken $FQ $NAM $dat viral $K krakenuid
      done
    done
  done
done
