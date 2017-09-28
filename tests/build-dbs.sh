#!/bin/bash

set -xeu

[[ "$#" -ne 1 ]] && DIR=`pwd` || DIR=$1

mkdir -p $DIR/dbs/refseq-viral-plus/library
[[ -L $DIR/dbs/refseq-viral-plus/library/viral ]] || ln -s $DIR/data/library/viral/ $DIR/dbs/refseq-viral-plus/library/
[[ -L $DIR/dbs/refseq-viral-plus/library/viral-neighbors ]] || ln -s $DIR/data/library/viral-neighbors/ $DIR/dbs/refseq-viral-plus/library/

export PATH="$DIR/install:$PATH"
krakenu-build --db $DIR/dbs/refseq-viral --build --taxids-for-genomes --taxids-for-sequences --library-dir=$DIR/data/library/viral --taxonomy-dir=$DIR/data/taxonomy

