#!/bin/bash

# Original file Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
# Portions (c) 2017, Florian Breitwieser <fbreitwieser@jhu.edu> as part of KrakenHLL
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

# Check that jellyfish is executable and is proper version
# Designed to be called by kraken-build

set -u  # Protect against uninitialized vars.
set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands

JELLYFISH1_BIN=""
for JF in $JELLYFISH_BIN $(dirname $0)/jellyfish/bin/jellyfish  /usr/local/opt/jellyfish-1.1/bin/jellyfish jellyfish1 jellyfish; do
  if test -f $JF || command -v $JF >/dev/null 2>&1; then
    JELLYFISH1_BIN=$JF;
    break
  fi
done

if [ "$JELLYFISH1_BIN" == "" ]; then
    echo "Did not find jellyfish!" 1>&2
    exit 1
fi

JELLYFISH_VERSION=$( $JELLYFISH1_BIN --version | awk '{print $2}')
if [[ $JELLYFISH_VERSION =~ ^1\. ]]
then
  echo "Found jellyfish v$JELLYFISH_VERSION" 1>&2
else
  echo "Found jellyfish v$JELLYFISH_VERSION" 1>&2
  echo "Kraken requires jellyfish version 1" 1>&2
  exit 1
fi
echo $JELLYFISH1_BIN
