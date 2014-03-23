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

if ! (grep '^>' "$1" | perl -nle 'exit 1 if ! /^>gi\|\d+\|/')
then
  echo "Can't add \"$1\": sequence is missing GI number"
  exit 1
fi

seq_ct=$(grep -c -m2 '^>' "$1")

add_dir="$LIBRARY_DIR/added"
mkdir -p "$add_dir"

free_num=0
if [ -e "$add_dir/next.free" ]
then
  free_num=$(cat "$add_dir/next.free")
fi
filename=$(printf '%015g' "$free_num")

if (( seq_ct > 1 ))
then
  cp "$1" "$add_dir/$filename.ffn"
  fasta_split.pl "$add_dir/$filename.ffn"
else
  cp "$1" "$add_dir/$filename.fna"
fi

free_num=$(( free_num + 1 ))
echo "$free_num" > "$add_dir/next.free"

echo "Added \"$1\" to library ($KRAKEN_DB_NAME)"
