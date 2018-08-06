#!/bin/bash

# Portions (c) 2017-2018, Florian Breitwieser <fbreitwieser@jhu.edu>
# Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
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

set -e

DIR=$(dirname $0)
VERSION=`cat $(dirname $0)/VERSION`
INSTALL_JELLYFISH=0
MAKE_ARGS=
MAKE_CLEAN="clean"
ADD_DEBUG_INFO=0

USAGE="Usage: $(basename $0) [OPTIONS] INSTALL_DIR

OPTIONS:
    -j          Install jellyfish v1.1 in INSTALL_DIR, too.
    -l BIN_DIR  Link KrakenHLL executables to BIN_DIR, e.g /usr/local/bin or ~/bin.
    -c BIN      Use compiler BIN instead of g++.
	-g          Add debug info
    -h          This help message

On MacOS, if you experience the error \"clang: fatal error: unsupported option '-fopenmp'\" on OSX, try installing g++ with brew, and using the option \"-c g++-7\".
"


while getopts "Chjc:g" OPTION; do
    case $OPTION in
    c) MAKE_ARGS="CXX=\"$OPTARG\"" ;;
    C) MAKE_CLEAN="" ;;
    j) INSTALL_JELLYFISH=1 ;;
	g) ADD_DEBUG_INFO=1 ;;
    h) echo "$USAGE"; exit 0 ;;
    *) echo "Incorrect options provided. $USAGE"
       exit 1 ;;
    esac
done
shift $((OPTIND -1))

if [ -z "$1" ]
then
  echo "$USAGE"
  exit 0
fi

if [ "$1" = "KRAKEN_DIR" ]
then
  echo "Please replace \"KRAKEN_DIR\" with the name of the directory"
  echo "that you want to install Kraken in."
  exit 1
fi


# Perl cmd used to canonicalize dirname - "readlink -f" doesn't work
# on OS X.
export KRAKEN_DIR=$(perl -MCwd=abs_path -le 'print abs_path(shift)' "$1")

mkdir -p "$KRAKEN_DIR"
if [ "$INSTALL_JELLYFISH" == "1" ]; then
  WD=`pwd`
  cd $KRAKEN_DIR
  if [[ ! -d jellyfish ]]; then
    wget http://www.cbcb.umd.edu/software/jellyfish/jellyfish-1.1.11.tar.gz
    tar xf jellyfish-1.1.11.tar.gz
    mv jellyfish-1.1.11 jellyfish-install
  fi
  cd jellyfish-install
  [[ -f Makefile ]] || ./configure
  make
  #make install ## doest not work for me on OSX
  cp $KRAKEN_DIR/jellyfish-install/bin/jellyfish $KRAKEN_DIR
  #rm -r jellyfish-1.1.11.tar.gz jellyfish-1.1.11
  cd $WD
fi

[[ "$ADD_DEBUG_INFO" == 1 ]] && MAKE_ARGS="$MAKE_ARGS NDEBUG=-g"
echo make -C $DIR/src $MAKE_CLEAN install $MAKE_ARGS
make -C $DIR/src $MAKE_CLEAN  install $MAKE_ARGS || { echo "Error building KrakenHLL. See $(basename $0) -h for options." >&2; exit 1; }
for file in $DIR/scripts/*
do
  [[ -f $file ]] || continue;
  perl -pl -e 'BEGIN { while (@ARGV) { $_ = shift; ($k,$v) = split /=/, $_, 2; $H{$k} = $v } }'\
           -e 's/#####=(\w+)=#####/$H{$1}/g' \
           "KRAKEN_DIR=$KRAKEN_DIR" "VERSION=$VERSION" \
           < "$file" > "$KRAKEN_DIR/$(basename $file)"
  if [ -x "$file" ]
  then
    chmod +x "$KRAKEN_DIR/$(basename $file)"
  fi
done
mkdir -p $KRAKEN_DIR/File
cp -r $DIR/scripts/File/* $KRAKEN_DIR/File

echo -n "
Kraken installation complete.

To make things easier for you, you may want to copy/symlink the following
files into a directory in your PATH:

ln -s"
for file in $KRAKEN_DIR/krakenhll*
do
  [ -x "$file" ] && echo -n " $file"
done
echo " DEST_DIR"
exit 0
