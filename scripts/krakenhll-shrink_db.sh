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

# Shrink a Kraken database
# Designed to be called by kraken_build

set -u  # Protect against uninitialized vars.
set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands

new_ct="$1"
new_db="$2"
offset="$3"

OLD_DB_DIR="$KRAKEN_DB_NAME"
NEW_DB_DIR="$new_db"

if [ -e "$NEW_DB_DIR" ]
then
  echo "$new_db already exists ($NEW_DB_DIR), aborting shrink operation."
  exit 1
else
  mkdir -p "$NEW_DB_DIR/taxonomy"
fi

cp "$OLD_DB_DIR/taxonomy/nodes.dmp" "$NEW_DB_DIR/taxonomy"
cp "$OLD_DB_DIR/taxonomy/names.dmp" "$NEW_DB_DIR/taxonomy"
db_shrink -n $new_ct -d "$OLD_DB_DIR/database.kdb" \
  -o "$NEW_DB_DIR/database.jdb.tmp" -O "$offset"
mv "$NEW_DB_DIR/database.jdb.tmp" "$NEW_DB_DIR/database.jdb"
echo "Reduced database created, now sorting..."
db_sort -M -t $KRAKEN_THREAD_CT -n $KRAKEN_MINIMIZER_LEN \
  -d "$NEW_DB_DIR/database.jdb" -o "$NEW_DB_DIR/database.kdb.tmp" \
  -i "$NEW_DB_DIR/database.idx"
mv "$NEW_DB_DIR/database.kdb.tmp" "$NEW_DB_DIR/database.kdb"
echo "Sort complete, database is ready."
