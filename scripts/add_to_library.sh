#!/bin/bash

# Copyright 2013, Derrick Wood <dwood@cs.umd.edu>
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

# Copy specified file into a Kraken library

set -u  # Protect against uninitialized vars.
set -e  # Stop on error

LIBRARY_DIR="$KRAKEN_DB_NAME/library"

if [ ! -e "$1" ]
then
  echo "Can't add \"$1\": file does not exist"
  exit 1
fi
if [ ! -f "$1" ]
then
  echo "Can't add \"$1\": not a regular file"
  exit 1
fi

if ! head -1 "$1" | perl -nle 'exit 1 unless /^>gi\|(\d+)\|/'
then
  echo "Can't add \"$1\": could not find GI number"
  exit 1
fi
seq_ct=$(grep -m2 '^>' "$1" | wc -l)
if (( seq_ct > 1 ))
then
  echo "Can't add \"$1\": multiple sequences found"
  exit 1
fi

mkdir -p "$LIBRARY_DIR/added"
ct=0
freefile=""
while [ -z "$freefile" ]
do
  freefile=$(seq -f '%015g' $ct $ct).fna
  if [ -e "$LIBRARY_DIR/added/$freefile" ]
  then
    ct=$(($ct + 1))
    freefile=""
  fi
done
cp "$1" "$LIBRARY_DIR/added/$freefile"
echo "Added \"$1\" to library ($KRAKEN_DB_NAME)"
