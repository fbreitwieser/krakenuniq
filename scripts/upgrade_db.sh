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

# Upgrade a pre-v0.10.0-beta Kraken DB to use scrambled minimizer order
# Designed to be called by kraken_build

set -u  # Protect against uninitialized vars.
set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands

function report_time_elapsed() {
  curr_time=$(date "+%s.%N")
  perl -e '$time = $ARGV[1] - $ARGV[0];' \
       -e '$sec = int($time); $nsec = $time - $sec;' \
       -e '$min = int($sec/60); $sec %= 60;' \
       -e '$hr = int($min/60); $min %= 60;' \
       -e 'print "${hr}h" if $hr;' \
       -e 'print "${min}m" if $min || $hr;' \
       -e 'printf "%.3fs", $sec + $nsec;' \
       $1 $curr_time
}

start_time=$(date "+%s.%N")

DATABASE_DIR="$KRAKEN_DB_NAME"

if [ ! -d "$DATABASE_DIR" ]
then
  echo "Can't find Kraken DB directory \"$KRAKEN_DB_NAME\""
  exit 1
fi
cd "$DATABASE_DIR"

MEMFLAG=""
if [ -n "$KRAKEN_WORK_ON_DISK" ]
then
  MEMFLAG="-M"
fi

if [ -e "old_database.kdb" ]
then
  echo "old_database.kdb found - it appears database is already upgraded."
  exit 1
fi

if [ ! -e "database.kdb" ]
then
  echo "Can't find database.kdb file!"
  exit 1
fi
if [ ! -e "database.idx" ]
then
  echo "Can't find database.idx file!"
  exit 1
fi

idx_size=$(stat -c '%s' database.idx)
# Calculate minimizer length based on existing index size
minimizer_len=$(perl -le 'print int(log(shift() / 8 - 2) / log(4))' $idx_size)

db_sort $MEMFLAG -t $KRAKEN_THREAD_CT -n $minimizer_len -d database.kdb \
  -o newdb.kdb -i newdb.idx
mv database.idx old_database.idx
mv database.kdb old_database.kdb
mv newdb.idx database.idx
mv newdb.kdb database.kdb

echo "Database upgrade complete. [$(report_time_elapsed $start_time)]"
echo "Old database files are at $KRAKEN_DB_NAME/old_database.{kdb,idx}"
