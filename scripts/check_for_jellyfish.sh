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

# Check that jellyfish is executable and is proper version
# Designed to be called by kraken-build

set -u  # Protect against uninitialized vars.
set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands

JELLYFISH_VERSION=$(jellyfish --version | awk '{print $2}')
if [[ $JELLYFISH_VERSION =~ ^1\. ]]
then
  echo "Found jellyfish v$JELLYFISH_VERSION"
else
  echo "Found jellyfish v$JELLYFISH_VERSION"
  echo "Kraken requires jellyfish version 1"
  exit 1
fi
