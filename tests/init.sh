#!/bin/bash

DIR=$1
[[ "$DIR" == "" ]] && DIR=`pwd`

## Install KrakenU locally into install/
$(dirname $0)/../install_kraken.sh --install-jellyfish $DIR/install

## Download taxonomy and genomic data into data/
#$DIR/install/krakenu-download --db $DIR/data -R --include-viral-neighbors taxonomy refseq/archaea refseq/bacteria refseq/viral/Any

for i in viral viral-neighbors archaea bacteria; do 
  if [[ ! -f "$DIR/data/all-$i.fna" ]]; then 
    find $DIR/data/library/$i -name '*.fna' -exec cat {} \; > $DIR/data/all-$i.fna
  fi
  if [[ ! -f "$DIR/data/all-$i.map" ]]; then 
    find $DIR/data/library/$i -name '*.map' -exec cat {} \; > $DIR/data/all-$i.map
  fi
done
