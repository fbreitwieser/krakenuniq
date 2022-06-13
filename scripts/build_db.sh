#!/bin/bash
#vim: noai:ts=2:sw=2

# Original file Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
# Portions (c) 2017, Florian Breitwieser <fbreitwieser@jhu.edu> as part of KrakenUniq
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

export VERBOSE=1

LOG_FILE=database-build.log
exe () {
  echo -e "$(date "+%Y-%m-%d %H:%M:%S")\t$@" >> $LOG_FILE
  $@
}
                             

function cmd () {
  export start_time1=$(date "+%s.%N")
  if [[ $VERBOSE -eq 1 ]]; then
    echo "EXECUTING $@"
  fi
  $@
}


start_time=$(date "+%s.%N")
script_dir=`dirname $0`

DATABASE_DIR="$KRAKEN_DB_NAME"
FIND_OPTS=-L
[[ "$JELLYFISH_BIN" == "" ]] && JELLYFISH_BIN=`$script_dir/check_for_jellyfish.sh`
NCBI_SERVER="ftp.ncbi.nih.gov"
FTP_SERVER="ftp://$NCBI_SERVER"

if [ ! -d "$DATABASE_DIR" ]
then
  echo "Can't find KrakenUniq DB directory \"$KRAKEN_DB_NAME\""
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
  rm -f database.* *.map lca.complete library-files.txt uid_database.* taxDB
fi

LIBRARY_DIR="library/"
[[ "$KRAKEN_LIBRARY_DIRS" != "" ]] && LIBRARY_DIR="$KRAKEN_LIBRARY_DIRS"

TAXONOMY_DIR="taxonomy/"
[[ "$KRAKEN_TAXONOMY_DIR" != "" ]] && TAXONOMY_DIR="$KRAKEN_TAXONOMY_DIR"

if [ ! -s "library-files.txt" ]; then
    echo "Finding all library files"
    find $FIND_OPTS $LIBRARY_DIR '(' -iname '*.fna' -o -iname '*.fa' -o -iname '*.ffn' -o -iname '*.fasta' -o -iname '*.fsa' ')' > library-files.txt
fi

cat_library() {
  cat library-files.txt | tr '\n' '\0' | xargs -0 cat
}

cat_libraryp() {
  find $FIND_OPTS $LIBRARY_DIR/{$@,} '(' -iname '*.fna' -o -iname '*.fa' -o -iname '*.ffn' -o -iname '*.fasta' -o -iname '*.fsa' ')' | tr '\n' '\0' | xargs -0 cat
}

N_FILES=`cat library-files.txt | wc -l`
if [[ "$N_FILES" -eq 0 ]]; then
  echo "ERROR: No fna, fa, or ffn files found in $LIBRARY_DIR!";
  exit 1
fi
echo "Found $N_FILES sequence files (*.{fna,fa,ffn,fasta,fsa}) in the library directory."


if [ -e "database.jdb" ] || [ -e "database0.kdb" ]
then
  echo "Skipping step 1, k-mer set already exists."
else
  echo "Creating k-mer set (step 1 of 6)..."
  start_time1=$(date "+%s.%N")

  echo "Using $JELLYFISH_BIN"
  [[ "$JELLYFISH_BIN" != "" ]] || exit 1
  if [ -z "$KRAKEN_HASH_SIZE" ]
  then
    KRAKEN_HASH_SIZE=$( cat_library | count_unique -t $KRAKEN_THREAD_CT -k $KRAKEN_KMER_LEN)
    echo "Hash size not specified, using '$KRAKEN_HASH_SIZE'"
  fi

  exe eval $JELLYFISH_BIN count -m $KRAKEN_KMER_LEN -s $KRAKEN_HASH_SIZE -C -t $KRAKEN_THREAD_CT \
    -o database <( cat_library )

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
      exe eval db_shrink -d database.jdb -o database.jdb.small -n $new_ct
      mv database.jdb database.jdb.big.tmp
      mv database.jdb.small database.jdb
      mv database.jdb.big.tmp database.jdb.big
      echo "Database reduced. [$(report_time_elapsed $start_time1)]"
    fi
  fi
fi

SORTED_DB_NAME=database0.kdb
if [ -e "$SORTED_DB_NAME" ]
then
  echo "Skipping step 3, k-mer set already sorted."
else
  echo "Sorting k-mer set (step 3 of 6)..."
  start_time1=$(date "+%s.%N")
  exe eval db_sort -z $MEMFLAG -t $KRAKEN_THREAD_CT -n $KRAKEN_MINIMIZER_LEN \
    -d database.jdb -o $SORTED_DB_NAME.tmp \
    -i database.idx

  # Once here, DB is sorted, can put file in proper place.
  mv $SORTED_DB_NAME.tmp $SORTED_DB_NAME

  echo "K-mer set sorted. [$(report_time_elapsed $start_time1)]"
fi

if [ -s "seqid2taxid.map" ]
then
  echo "Skipping step 4, seqID to taxID map already complete."
else
  echo "Creating seqID to taxID map (step 4 of 6).."
  start_time1=$(date "+%s.%N")
  find -L $LIBRARY_DIR/ -name '*.map' -exec cat {} \; > seqid2taxid.map
  line_ct=$(wc -l seqid2taxid.map | awk '{print $1}')

  echo "$line_ct sequences mapped to taxa. [$(report_time_elapsed $start_time1)]"
fi


if [ -s "taxDB" ]
then
  echo "Skipping step 5, taxDB exists."
else
  echo "Creating taxDB (step 5 of 6)... "
  start_time1=$(date "+%s.%N")
  if [ ! -f $TAXONOMY_DIR/names.dmp ] || [ ! -f $TAXONOMY_DIR/nodes.dmp ]; then
    echo "$TAXONOMY_DIR/names.dmp or $TAXONOMY_DIR/nodes.dmp does not exist - downloading it ..."
    [ -d $TAXONOMY_DIR ] || mkdir $TAXONOMY_DIR
    cd $TAXONOMY_DIR
    exe eval wget $FTP_SERVER/pub/taxonomy/taxdump.tar.gz
    exe eval tar zxf taxdump.tar.gz
    cd ..
  fi
  build_taxdb $TAXONOMY_DIR/names.dmp $TAXONOMY_DIR/nodes.dmp | sort -t$'\t' -rnk6,6 -rnk5,5 > taxDB.tmp
  mv taxDB.tmp taxDB
  echo "taxDB construction finished. [$(report_time_elapsed $start_time1)]"
fi

if [ "$KRAKEN_LCA_DATABASE" != "0" ]; then
  if [ -s "database.kdb" ] && [ "$KRAKEN_RESET_TAXIDS" != "1" ]
  then
    echo "Skipping step 6, LCAs already set."
  else
    echo "Building  KrakenUniq LCA database (step 6 of 6)..."
    PARAM=""
    if [[ "$KRAKEN_ADD_TAXIDS_FOR_SEQ" == "1" ]]; then
  	echo " Adding taxonomy IDs for sequences"
  	PARAM=" -a"
    fi
    if [[ "$KRAKEN_ADD_TAXIDS_FOR_GENOME" == "1" ]]; then
  	echo " Adding taxonomy IDs for genomes"
  	PARAM="$PARAM -A"
    fi
	if [[ "$KRAKEN_MIN_CONTIG_SIZE" != "" ]]; then
	  echo "Excluding sequences smaller than $KRAKEN_MIN_CONTIG_SIZE"
	  PARAM1="-E $KRAKEN_MIN_CONTIG_SIZE"
    else PARAM1=""
	fi
    start_time1=$(date "+%s.%N")

    if [ "$KRAKEN_ADD_TAXIDS_FOR_SEQ" == "1" ] || [ "$KRAKEN_ADD_TAXIDS_FOR_GENOME" == "1" ]; then
        cp taxDB taxDB.orig
    fi

	[[ -z "${KRAKEN_LCA_ORDER}" ]] && DC="-c database.kdb.counts" || DC=""
    set_lcas $MEMFLAG -x -d $SORTED_DB_NAME -o database.kdb -i database.idx -v \
        -b taxDB $PARAM $PARAM1 -t $KRAKEN_THREAD_CT -m seqid2taxid.map $DC \
        -F <( cat_library ) -T > seqid2taxid-plus.map
    if [ "$KRAKEN_ADD_TAXIDS_FOR_SEQ" == "1" ] || [ "$KRAKEN_ADD_TAXIDS_FOR_GENOME" == "1" ]; then
      mv seqid2taxid.map seqid2taxid.map.orig
      mv seqid2taxid-plus.map seqid2taxid.map
    else
      rm seqid2taxid-plus.map
    fi
	if [[ ! -z "${KRAKEN_LCA_ORDER}" ]]; then
      echo "  Re-setting LCA's hierarchically"
	  IFS=';' read -ra MYDIRS <<< "$KRAKEN_LCA_ORDER"
      let COUNTER=1
      let TOTAL=${#MYDIRS[@]}

	  for DDIR in "${MYDIRS[@]}"; do
        echo " Setting LCAs for $DDIR (substep 6.$COUNTER of 6.$TOTAL) ..."
	    [[ $COUNTER -eq $TOTAL ]] && DC="-c database.kdb.counts" || DC=""
		  ## First reset all taxids that appear in the set to zero (flag -R)
        exe eval set_lcas $MEMFLAG -x -d $SORTED_DB_NAME -o database.kdb -i database.idx -v \
          -b taxDB $PARAM1 -t $KRAKEN_THREAD_CT -m seqid2taxid.map \
          -F <( cat_libraryp $DDIR ) -TR

		  ## Then just re-set them
        exe eval set_lcas $MEMFLAG -x -d $SORTED_DB_NAME -o database.kdb -i database.idx -v \
          -b taxDB $PARAM1 -t $KRAKEN_THREAD_CT -m seqid2taxid.map $DC \
          -F <( cat_libraryp $DDIR ) -T

        COUNTER=$((COUNTER+1))
	  done

	fi

    echo "LCA database created. [$(report_time_elapsed $start_time1)]"
  fi
  ## Make a classification report
  REPNAME=database
  if [[ ! -s $REPNAME.report.tsv ]]; then
    echo "Creating database summary report $REPNAME.report.tsv ..."
    krakenuniq --db . --report-file $REPNAME.report.tsv --threads $KRAKEN_THREAD_CT $(cat library-files.txt | tr '\n' ' ') > $REPNAME.kraken.tsv
  fi
fi


if [ "$KRAKEN_UID_DATABASE" != "0" ]; then
  if [ -s "uid_database.kdb" ]
  then
    echo "Skipping step 6.3, UID datanbase already generated."
  else
    echo "Building UID database (step 6.3 of 6)..."
    PARAM=""
    if [[ "$KRAKEN_LCA_DATABASE" == "0" ]]; then
      if [[ "$KRAKEN_ADD_TAXIDS_FOR_SEQ" == "1" ]]; then
  	echo " Adding taxonomy IDs for sequences"
  	PARAM=" -a"
      fi
      if [[ "$KRAKEN_ADD_TAXIDS_FOR_GENOME" == "1" ]]; then
   	echo " Adding taxonomy IDs for genomes"
    	PARAM="$PARAM -A"
      fi
    fi
    start_time1=$(date "+%s.%N")
      set_lcas $MEMFLAG -x -d $SORTED_DB_NAME -I uid_to_taxid.map -o uid_database.kdb -i database.idx -v \
        -b taxDB $PARAM -t $KRAKEN_THREAD_CT -m seqid2taxid.map -c uid_database.kdb.counts -F <( cat_library )
  
    echo "UID Database created. [$(report_time_elapsed $start_time1)]"
  fi

  ## Make a classification report
  REPNAME=uid_database
  if [[ ! -s $REPNAME.report.tsv ]]; then
    echo "Creating UID database summary report $REPNAME.report.tsv ..."
    krakenuniq --db . --report-file $REPNAME.report.tsv --threads $KRAKEN_THREAD_CT --uid-mapping $(cat library-files.txt | tr '\n' ' ') > $REPNAME.kraken.tsv
  fi
fi

echo "Database construction complete. [Total: $(report_time_elapsed $start_time)]
You can delete all files but database.{kdb,idx} and taxDB now, if you want"


