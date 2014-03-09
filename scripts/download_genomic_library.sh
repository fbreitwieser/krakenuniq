#!/bin/bash

# Copyright 2013-2014, Derrick Wood <dwood@cs.umd.edu>
#
# This file is part of the Kraken taxonomic sequence classification system.
#
# Kraken is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Kraken is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Kraken.  If not, see <http://www.gnu.org/licenses/>.

# Download specific genomic libraries for use with Kraken.
# Supported choices are:
#   bacteria - NCBI RefSeq complete bacterial/archaeal genomes
#   viruses - NCBI RefSeq complete viral DNA and RNA genomes
#   human - NCBI RefSeq GRCh38 human reference genome

set -u  # Protect against uninitialized vars.
set -e  # Stop on error

LIBRARY_DIR="$KRAKEN_DB_NAME/library"
NCBI_SERVER="ftp.ncbi.nih.gov"
FTP_SERVER="ftp://$NCBI_SERVER"
RSYNC_SERVER="rsync://$NCBI_SERVER"
THIS_DIR=$PWD

case "$1" in
  "bacteria")
    mkdir -p $LIBRARY_DIR/Bacteria
    cd $LIBRARY_DIR/Bacteria
    if [ ! -e "lib.complete" ]
    then
      rm -f all.fna.tar.gz
      wget $FTP_SERVER/genomes/Bacteria/all.fna.tar.gz
      echo -n "Unpacking..."
      tar zxf all.fna.tar.gz
      rm all.fna.tar.gz
      echo " complete."
      touch "lib.complete"
    else
      echo "Skipping download of bacterial genomes, already downloaded here."
    fi
    ;;
  "viruses")
    mkdir -p $LIBRARY_DIR/Viruses
    cd $LIBRARY_DIR/Viruses
    if [ ! -e "lib.complete" ]
    then
      rm -f all.fna.tar.gz
      rm -f all.ffn.tar.gz
      wget $FTP_SERVER/genomes/Viruses/all.fna.tar.gz
      wget $FTP_SERVER/genomes/Viruses/all.ffn.tar.gz
      echo -n "Unpacking..."
      tar zxf all.fna.tar.gz
      tar zxf all.ffn.tar.gz
      rm all.fna.tar.gz
      rm all.ffn.tar.gz
      echo " complete."
      touch "lib.complete"
    else
      echo "Skipping download of viral genomes, already downloaded here."
    fi
    ;;
  "human")
    mkdir -p $LIBRARY_DIR/H_sapiens
    if [ ! -e "$LIBRARY_DIR/H_sapiens/lib.complete" ]
    then
      rsync -avz --include "H_sapiens/" --include "*/CHR_*/" --include "hs_ref_GRCh*.fa.gz" --exclude "*" \
        $RSYNC_SERVER/genomes/H_sapiens $LIBRARY_DIR
      for gzfile in $LIBRARY_DIR/H_sapiens/*/*.gz
      do
        dir=$(dirname $gzfile)
        cd $dir
        gunzip $(basename $gzfile)
        echo "Unzipped $(basename $gzfile)"
        cd $THIS_DIR
      done
      for fafile in $LIBRARY_DIR/H_sapiens/*/*.fa
      do
        mv $fafile ${fafile%.fa}.fna
        echo "Renamed $(basename $fafile)"
      done
      touch "$LIBRARY_DIR/H_sapiens/lib.complete"
    else
      echo "Skipping download of human genome, already downloaded here."
    fi
    ;;
  *)
    echo "Unsupported library.  Valid options are: "
    echo "  bacteria virus human"
    ;;
esac
