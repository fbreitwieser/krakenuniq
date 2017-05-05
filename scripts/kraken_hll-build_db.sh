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

# Build a Kraken database
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
script_dir=`dirname $0`

DATABASE_DIR="$KRAKEN_DB_NAME"
FIND_OPTS=-L
JELLYFISH_BIN=`$script_dir/kraken_hll-check_for_jellyfish.sh`

if [ ! -d "$DATABASE_DIR" ]
then
  echo "Can't find Kraken DB directory \"$KRAKEN_DB_NAME\""
  exit 1
fi
cd "$DATABASE_DIR"

MEMFLAG=""
if [ -z "$KRAKEN_WORK_ON_DISK" ]
then
  MEMFLAG="-M"
  echo "Kraken build set to minimize disk writes."
else
  echo "Kraken build set to minimize RAM usage."
fi

if [ "$KRAKEN_REBUILD_DATABASE" == "1" ]
then
  rm -f database.* *.map lca.complete library/seq-files.txt
fi

if [ !-f "library/seq-files.txt" ]; then
    echo "Finding all library files"
    find $FIND_OPTS library/ '(' -name '*.fna' -o -name '*.fa' -o -name '*.ffn' ')' > library/seq-files.txt
fi
N_FILES=`cat library/seq-files.txt | wc -l`
echo "Found $N_FILES sequence files (*.{fna,fa,ffn} in the library)"

if [ -e "database.jdb" ]
then
  echo "Skipping step 1, k-mer set already exists."
else
  echo "Creating k-mer set (step 1 of 6)..."
  start_time1=$(date "+%s.%N")

  echo "Using $JELLYFISH_BIN"
  [[ "$JELLYFISH_BIN" != "" ]] || exit 1
  # Estimate hash size as 1.15 * chars in library FASTA files
  if [ -z "$KRAKEN_HASH_SIZE" ]
  then
    KRAKEN_HASH_SIZE=$(find $FIND_OPTS library/ '(' -name '*.fna' -o -name '*.fa' -o -name '*.ffn' ')' -printf '%s\n' | perl -nle '$sum += $_; END {print int(1.15 * $sum)}')
    echo "Hash size not specified, using '$KRAKEN_HASH_SIZE'"
  fi

  cat library/seq-files.txt | tr '\n' '\0' | xargs -0 cat | \
    $JELLYFISH_BIN count -m $KRAKEN_KMER_LEN -s $KRAKEN_HASH_SIZE -C -t $KRAKEN_THREAD_CT \
      -o database /dev/fd/0

  # Merge only if necessary
  if [ -e "database_1" ]
  then
    $JELLYFISH_BIN merge -o database.jdb.tmp database_*
  else
    mv database_0 database.jdb.tmp
  fi

  # Once here, DB is finalized, can put file in place.
  mv database.jdb.tmp database.jdb

  echo "K-mer set created. [$(report_time_elapsed $start_time1)]"
fi

if [ -z "$KRAKEN_MAX_DB_SIZE" ]
then
  echo "Skipping step 2, no database reduction requested."
else
  if [ -e "database.jdb.big" ]
  then
    echo "Skipping step 2, database reduction already done."
  else
    start_time1=$(date "+%s.%N")
    kdb_size=$(stat -c '%s' database.jdb)
    idx_size=$(echo "8 * (4 ^ $KRAKEN_MINIMIZER_LEN + 2)" | bc)
    resize_needed=$(echo "scale = 10; ($kdb_size+$idx_size)/(2^30) > $KRAKEN_MAX_DB_SIZE" | bc)
    if (( resize_needed == 0 ))
    then
      echo "Skipping step 2, database reduction unnecessary."
    else
      echo "Reducing database size (step 2 of 6)..."
      max_kdb_size=$(echo "$KRAKEN_MAX_DB_SIZE*2^30 - $idx_size" | bc)
      idx_size_gb=$(printf %.2f $(echo "$idx_size/2^30" | bc) )
      if (( $(echo "$max_kdb_size < 0" | bc) == 1 ))
      then
        echo "Maximum database size too small - index alone needs $idx_size_gb GB.  Aborting reduction."
        exit 1
      fi
      # Key ct is 8 byte int stored 48 bytes from start of file
      key_ct=$(perl -MFcntl -le 'open F, "database.jdb"; seek F, 48, SEEK_SET; read F, $b, 8; $a = unpack("Q", $b); print $a')
      # key_bits is 8 bytes from start
      key_bits=$(perl -MFcntl -le 'open F, "database.jdb"; seek F, 8, SEEK_SET; read F, $b, 8; $a = unpack("Q", $b); print $a')
      # this is basically ceil(key_bits / 8) - why no ceiling function, bc?
      key_len=$(echo "($key_bits + 7) / 8" | bc)
      # val_len is 16 bytes from start
      val_len=$(perl -MFcntl -le 'open F, "database.jdb"; seek F, 16, SEEK_SET; read F, $b, 8; $a = unpack("Q", $b); print $a')
      record_len=$(( key_len + val_len ))
      new_ct=$(echo "$max_kdb_size / $record_len" | bc)
      echo "Shrinking DB to use only $new_ct of the $key_ct k-mers"
      db_shrink -d database.jdb -o database.jdb.small -n $new_ct
      mv database.jdb database.jdb.big.tmp
      mv database.jdb.small database.jdb
      mv database.jdb.big.tmp database.jdb.big
      echo "Database reduced. [$(report_time_elapsed $start_time1)]"
    fi
  fi
fi

if [ -e "database.kdb" ]
then
  echo "Skipping step 3, k-mer set already sorted."
else
  echo "Sorting k-mer set (step 3 of 6)..."
  start_time1=$(date "+%s.%N")
  db_sort -z $MEMFLAG -t $KRAKEN_THREAD_CT -n $KRAKEN_MINIMIZER_LEN \
    -d database.jdb -o database.kdb.tmp \
    -i database.idx

  # Once here, DB is sorted, can put file in proper place.
  mv database.kdb.tmp database.kdb

  echo "K-mer set sorted. [$(report_time_elapsed $start_time1)]"
fi

if [ -e "seqid2taxid.map" ]
then
  echo "Skipping step 4, seqID to taxID map already complete."
else
  echo "Creating seqID to taxID map (step 4 of 6).."
  start_time1=$(date "+%s.%N")
  cat library/seq-files.txt | tr '\n' '\0' | xargs -0 grep '^>' | sed 's/.//' | sed 's/ .*//' | sort > library/seq-headers.txt
  join -t $'\t' nucl_gb.accession2taxid.sorted library/seq-headers.txt > seqid2taxid.map.tmp
  mv seqid2taxid.map.tmp seqid2taxid.map
  line_ct=$(wc -l seqid2taxid.map | awk '{print $1}')

  echo "$line_ct sequences mapped to taxa. [$(report_time_elapsed $start_time1)]"
fi

if [ -e "lca.complete" ]
then
  echo "Skipping step 5, LCAs already set."
else
  echo "Setting LCAs in database (step 5 of 6)..."
  PARAM=""
  if [[ "$KRAKEN_ADD_TAXIDS_FOR_SEQ" == "1" ]]; then
	echo " Adding taxonomy IDs for sequences"
	PARAM=" -a"
  fi
  start_time1=$(date "+%s.%N")
  cat library/seq-files.txt | tr '\n' '\0' | xargs -0 cat | \
    set_lcas $MEMFLAG -x -d database.kdb -i database.idx -v \
    -b taxDB $PARAM -t $KRAKEN_THREAD_CT -m seqid2taxid.map -F /dev/fd/0
  touch "lca.complete"

  echo "Database LCAs set. [$(report_time_elapsed $start_time1)]"
fi

if [ -s "taxDB" ]
then
  echo "Skipping step 6, taxDB exists."
else
  echo "Creating taxDB (step 6 of 6)... "
  jellyfish1 dump database.kdb | grep '^>' | sed 's/.//' | sort | uniq -c | sort -rn | sed 's/^ *\([0-9]\+\) \+\([0-9]\+\)$/\2\t\1/' > database.taxon_count
  build_taxdb taxonomy/names.dmp taxonomy/nodes.dmp database.taxon_count > taxDB.tmp
  mv taxDB.tmp taxDB
fi


echo "Database construction complete. [Total: $(report_time_elapsed $start_time)]
You can delete all files but database.{kdb,idx} and taxDB now, if you want"
