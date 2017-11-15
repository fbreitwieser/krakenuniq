#!/bin/bash

set -xeu

[[ $# -eq 1 ]] && DIR=$1 || DIR=`pwd`

## Install KrakenHLL locally into install/
#$(dirname $0)/../install_kraken.sh --install-jellyfish $DIR/install

## Download taxonomy and genomic data into data/
time krakenhll-download --db $DIR/data taxonomy refseq/archaea refseq/bacteria
time krakenhll-download --db $DIR/data --include-viral-neighbors refseq/viral/Any
time krakenhll-download --db $DIR/data refseq/fungi refseq/fungi/Chromosome refseq/protozoa refseq/protozoa/Chromosome
time krakenhll-download --db $DIR/data --fna rna,genomic refseq/vertebrate_mammalian/Chromosome/taxid9606 
time krakenhll-download --db $DIR/data contaminants

for i in fungi protozoa viral viral-neighbors archaea bacteria; do 
  [[ -s "$DIR/data/all-$i.fna" ]] || find $DIR/data/library/$i -name '*.fna' -print0 | xargs -0 -n 100 cat > $DIR/data/all-$i.fna
  [[ -s "$DIR/data/all-$i.map" ]] || find $DIR/data/library/$i -name '*.map' -print0 | xargs -0 -n 100 cat > $DIR/data/all-$i.map
  DUSTED_F="$DIR/data/all-$i-dusted.fna"
  [[ -s $DUSTED_F ]] || dustmasker -infmt fasta -in $DIR/data/all-$i.fna -level 20 -outfmt fasta | sed '/^>/! s/[^AGCT]/N/g' > "$DUSTED_F"
  mkdir -p $DIR/data/library/$i-dusted
  [[ -f "$DIR/data/library/$i-dusted/all-$i-dusted.fna" ]] || ln "$DUSTED_F" "$DIR/data/library/$i-dusted/all-$i-dusted.fna"
  [[ -f "$DIR/data/library/$i-dusted/all-$i.map" ]] || ln "$DIR/data/all-$i.map" "$DIR/data/library/$i-dusted/all-$i.map"
done
