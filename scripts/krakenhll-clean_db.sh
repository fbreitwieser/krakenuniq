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

# Clean unneeded files from a database

set -u  # Protect against uninitialized vars.
set -e  # Stop on error

cd "$KRAKEN_DB_NAME"

[ -e "database.kdb" ] || (echo "Incomplete database, clean aborted."; exit 1)
[ -e "database.idx" ] || (echo "Incomplete database, clean aborted."; exit 1)
[ -e "taxonomy/nodes.dmp" ] || (echo "Incomplete database, clean aborted."; exit 1)
[ -e "taxonomy/names.dmp" ] || (echo "Incomplete database, clean aborted."; exit 1)

rm -rf library
rm -f database.jdb* database_* *.map lca.complete 
mkdir newtaxo
mv taxonomy/{nodes,names}.dmp newtaxo
rm -rf taxonomy
mv newtaxo taxonomy
