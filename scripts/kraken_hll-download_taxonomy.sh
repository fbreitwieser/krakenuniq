#!/bin/bash

# Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
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

# Download NCBI taxonomy information for Kraken.
# Designed to be called by kraken_build

set -u  # Protect against uninitialized vars.
set -e  # Stop on error

TAXONOMY_DIR="$KRAKEN_DB_NAME/taxonomy"
NCBI_SERVER="ftp.ncbi.nih.gov"
FTP_SERVER="ftp://$NCBI_SERVER"
THIS_DIR=$PWD

mkdir -p "$TAXONOMY_DIR"
cd "$TAXONOMY_DIR"

if [ ! -e "nucl_gb.accession2taxid.flag" ]
then
  wget $FTP_SERVER/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
  time gunzip -c nucl_gb.accession2taxid.gz | cut -f 2,3 | sort -k 1,1 > nucl_gb.accession2taxid.sorted
  touch nucl_gb.accession2taxid.flag
  echo "Downloaded and sorted GB to taxon map"
fi

if [ ! -e "taxdump.flag" ]
then
  wget $FTP_SERVER/pub/taxonomy/taxdump.tar.gz
  tar zxf taxdump.tar.gz
  touch taxdump.flag
  echo "Downloaded and uncompressed taxonomy tree data"
fi
