#!/bin/bash

set -eu

#[[ "$#" -ne 1 ]] && DIR=`pwd` || DIR=$1
DIR=`pwd`
[[ `uname` == "Darwin" ]] && THREADS=4 || THREADS=10
TAXDB=""

build_db() {
  local PROG=$1; shift
  local K=$1; shift
  local MIN=$1; shift
  local NAM=$1; shift

  set -eu

  local DB_NAM=refseq-$NAM-k$K
  DB_DIR=$DIR/dbs-$PROG/$DB_NAM

  if [[ "$PROG" == "kraken" ]]; then
    mkdir -p $DB_DIR
    CMD="krakenuniq-build --kmer-len $K --minimizer-len $MIN --threads $THREADS --db $DB_DIR --build --taxids-for-genomes --taxids-for-sequences --taxonomy-dir=$DIR/data/taxonomy"
    for L in $@; do
      CMD="$CMD  --library-dir=$DIR/data/library/$L"
    done
  elif [[ "$PROG" == "kallisto" ]]; then
    CMD="kallisto index -k $K -i $DB_DIR"
    for L in $@; do
      CMD="$CMD  $DIR/data/all-$L.fna"
    done
  else 
    echo "Unknown command $PROG"
    exit 1
  fi
  if [[ -f "$TAXDB"  &&  ! -f "$DB_DIR/taxDB" ]]; then
    mkdir -p $DB_DIR
    echo "#################### Copying taxDB $TAXDB"
    cp -v $TAXDB $DB_DIR
  fi

  #if [[ ! -f "$DB_DIR-is.busy" ]]; then
    echo "EXECUTING $CMD"
    touch $DB_DIR-is.busy
    echo "#### KrakenBuild at "`date` >> $DB_DIR-build.log
    $CMD 2>&1 | tee -a $DB_DIR-build.log
    if [[ $PROG == "kraken" && ! -f "$DB_DIR/taxonomy/nodes.dmp" ]]; then
      mkdir -p $DB_DIR/taxonomy
      echo "EXECUTING dump_taxdb $DB_DIR/taxDB $DB_DIR/taxonomy/names.dmp $DB_DIR/nodes.dmp"
      dump_taxdb $DB_DIR/taxDB $DB_DIR/taxonomy/names.dmp $DB_DIR/taxonomy/nodes.dmp
    fi
    export TAXDB="$DB_DIR/taxDB"
    rm $DB_DIR-is.busy
  #else 
  #  echo "$DB_DIR-is.busy exists, ignoring directory."
  #fi
}



VERBOSE=false
HELP=false
DRY_RUN=false
K=31
M=15
THREADS=10
PATH1="."

USAGE="
`basename $0` [options] {kraken,kaiju} {viral|all-viral|prok|std|euk|archaea}

Options:
  -k KMER_SIZE       default $K
  -m MINIMIZER_SIZE  default $M
  -t THREADS         default $THREADS
"

OPTS=`getopt -o vhnk:t:p:m: --long verbose,dry-run,help,threads:,path: -n 'parse-options' -- "$@"`
if [ $? != 0 ] ; then echo "Failed parsing options. Usage: $USAGE" >&2 ; exit 1 ; fi
eval set -- "$OPTS"

while true; do
  case "$1" in
    -v | --verbose ) VERBOSE=true; shift ;;
    -h | --help )    HELP=true; shift ;;
    -n | --dry-run ) DRY_RUN=true; shift ;;
    -k | --kmer-size ) K="$2"; shift; shift ;;
    -m | --minimizer-size ) M="$2"; shift; shift ;;
    -t | --threads ) THREADS="$2"; shift; shift ;;
    -p | --path ) PATH1="$2"; shift; shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done
shift $((OPTIND -1))

if [[ "$#" -le 1 ]]; then
  echo "$USAGE"
  exit 1
fi

[[ "$PATH" != "" ]] && export PATH="$PATH1:$PATH"

PROG=$1
shift
for VAR in $@; do
  case "$VAR" in
    human)     build_db $PROG $K $M human contaminants vertebrate_mammalian ;;
    viral)     build_db $PROG $K $M viral viral ;;
    all-viral) build_db $PROG $K $M all-viral viral viral-neighbors  ;;
    prok)      build_db $PROG $K $M prok archaea-dusted bacteria-dusted ;;
    archaea)   build_db $PROG $K $M archaea archaea ;;
    std)   build_db $PROG $K $M std archaea-dusted bacteria-dusted viral-dusted viral-neighbors-dusted \
                               vertebrate_mammalian contaminants ;;
    euk)       build_db $PROG $K $M euk fungi-dusted protozoa-dusted ;;
  *) echo -e "Unknown database $VAR!\n$USAGE"
     exit 1 ;;
  esac
done

