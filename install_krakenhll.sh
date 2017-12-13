#!/bin/bash

# Portions (c) 2017, Florian Breitwieser <fbreitwieser@jhu.edu>
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

if [ "$1" == "--install-jellyfish" ]; then
 INSTALL_JELLYFISH=1;
 shift;
fi

if [ -z "$1" ] || [ -n "$2" ]
then
  echo "Usage: $(basename $0) [--install-jellyfish] KRAKEN_DIR

If --install-jellyfish is specified, the source code for version 1.1
is downloaded from http://www.cbcb.umd.edu/software/jellyfish and installed 
in KRAKEN_DIR. Note that this may overwrite other jellyfish installation in 
the same path."
  exit 64
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
    tar xvvf jellyfish-1.1.11.tar.gz
    mv jellyfish-1.1.11 jellyfish
  fi
  cd jellyfish
  [[ -f Makefile ]] || ./configure
  make
  #make install ## doest not work for me on OSX
  #cp $KRAKEN_DIR/jellyfish-install/bin/jellyfish $KRAKEN_DIR
  #rm -r jellyfish-1.1.11.tar.gz jellyfish-1.1.11
  cd $WD
fi

make -C src clean
make -C $DIR/src install
for file in $DIR/scripts/*
do
  perl -pl -e 'BEGIN { while (@ARGV) { $_ = shift; ($k,$v) = split /=/, $_, 2; $H{$k} = $v } }'\
           -e 's/#####=(\w+)=#####/$H{$1}/g' \
           "KRAKEN_DIR=$KRAKEN_DIR" "VERSION=$VERSION" \
           < "$file" > "$KRAKEN_DIR/$(basename $file)"
  if [ -x "$file" ]
  then
    chmod +x "$KRAKEN_DIR/$(basename $file)"
  fi
done

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
