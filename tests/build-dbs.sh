#!/bin/bash

set -xeu

[[ "$#" -ne 1 ]] && DIR=`pwd` || DIR=$1

mkdir -p $DIR/dbs/refseq-viral-plus/library
[[ -L $DIR/dbs/refseq-viral-plus/library/viral ]] || ln -s $DIR/data/library/viral/ $DIR/dbs/refseq-viral-plus/library/
[[ -L $DIR/dbs/refseq-viral-plus/library/viral-neighbors ]] || ln -s $DIR/data/library/viral-neighbors/ $DIR/dbs/refseq-viral-plus/library/

export PATH="$DIR/install:$PATH"
for K in 21 26 31; do
  mkdir -p $DIR/dbs/refseq-viral-k$K
  krakenu-build --kmer-len $K --minimizer-len 12 --threads 4 --db $DIR/dbs/refseq-viral-k$K --build --taxids-for-genomes --taxids-for-sequences --library-dir=$DIR/data/library/viral --taxonomy-dir=$DIR/data/taxonomy

  if [[ `uname` != "Darwin" ]]; then
    krakenu-build --kmer-len $K --threads 4 --db $DIR/dbs/refseq-bacteria-k$K --build --taxids-for-genomes --taxids-for-sequences --library-dir=$DIR/data/library/bacteria --taxonomy-dir=$DIR/data/taxonomy

  fi
done

