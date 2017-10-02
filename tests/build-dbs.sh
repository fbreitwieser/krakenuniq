#!/bin/bash

set -xeu

[[ "$#" -ne 1 ]] && DIR=`pwd` || DIR=$1

export PATH="$DIR/install:$PATH"
for K in 31 26 21; do
  mkdir -p $DIR/dbs/refseq-viral-k$K
  time krakenu-build --kmer-len $K --minimizer-len 12 --threads 4 --db $DIR/dbs/refseq-viral-k$K --build --taxids-for-genomes --taxids-for-sequences --library-dir=$DIR/data/library/viral --taxonomy-dir=$DIR/data/taxonomy 2>&1 | tee $DIR/dbs/refseq-viral-k$K/build.log

  mkdir -p $DIR/dbs/refseq-viral-k$K/taxonomy
  dump_taxdb $DIR/dbs/refseq-viral-k$K/taxDB $DIR/dbs/refseq-viral-k$K/taxonomy/names.dmp $DIR/dbs/refseq-viral-k$K/taxonomy/nodes.dmp

  if [[ `uname` != "Darwin" ]]; then
    mkdir -p $DIR/dbs/refseq-bacteria-k$K
    krakenu-build --kmer-len $K --threads 4 --db $DIR/dbs/refseq-bacteria-k$K --build --taxids-for-genomes --taxids-for-sequences --library-dir=$DIR/data/library/bacteria --library-dir=$DIR/data/library/archaea --taxonomy-dir=$DIR/data/taxonomy
    mkdir -p $DIR/dbs/refseq-oct2017-k$K
    krakenu-build --kmer-len $K --threads 4 --db $DIR/dbs/refseq-oct2017-k$K --build --taxids-for-genomes --library-dir=$DIR/data/library/viral-dusted --library-dir=$DIR/data/library/viral-neighbors-dusted --library-dir=$DIR/data/library/bacteria-dusted --library-dir=$DIR/data/library/archaea-dusted --library-dir=$DIR/data/libray/vertebrate_mammalia --taxonomy-dir=$DIR/data/taxonomy
  fi
done

