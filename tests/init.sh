#!/bin/bash

set -xeu

[[ $# -eq 1 ]] && DIR=$1 || DIR=`pwd`

## Install KrakenU locally into install/
#$(dirname $0)/../install_kraken.sh --install-jellyfish $DIR/install

## Download taxonomy and genomic data into data/
#$DIR/install/krakenu-download --db $DIR/data -R --include-viral-neighbors taxonomy refseq/archaea refseq/bacteria refseq/viral/Any
#$DIR/install/krakenu-download --db $DIR/data --fna rna,genomic -R refseq/vertebrate_mammalian/Chromosome/taxid9606 
$DIR/install/krakenu-download --db $DIR/data -R contaminants

for i in viral viral-neighbors archaea bacteria; do 
  [[ -s "$DIR/data/all-$i.fna" ]] || find $DIR/data/library/$i -name '*.fna' -exec cat {} \; > $DIR/data/all-$i.fna
  [[ -s "$DIR/data/all-$i.map" ]] || find $DIR/data/library/$i -name '*.map' -exec cat {} \; > $DIR/data/all-$i.map
  DUSTED_F="$DIR/data/all-$i-dusted.fna"
  [[ -s $DUSTED_F ]] || dustmasker -infmt fasta -in $DIR/data/all-$i.fna -level 20 -outfmt fasta | sed '/^>/! s/[^AGCT]/N/g' > "$DUSTED_F"
  mkdir -p $DIR/data/library/$i-dusted
  [[ -f "$DIR/data/library/$i-dusted/all-$i-dusted.fna" ]] || ln "$DUSTED_F" "$DIR/data/library/$i-dusted/all-$i-dusted.fna"
  [[ -f "$DIR/data/library/$i-dusted/all-$i.map" ]] || ln "$DIR/data/all-$i.map" "$DIR/data/library/$i-dusted/all-$i.map"
done
