#!/usr/bin/perl

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

# Reads multi-FASTA input and for each sequence ID reports a
# tab-delimited line:
#   <GI number> <sequence ID>
# 
#   or in the case of a sequence with Kraken taxid information:
#
#   TAXID <taxonomy ID> <sequence ID>
#
# Assumes all sequence IDs actually have GI numbers or Kraken
# taxid information.

use strict;
use warnings;
use File::Basename;

my $PROG = basename $0;

while (<>) {
  next unless /^>(\S+)/;
  my $seq_id = $1;
  if ($seq_id =~ /(^|\|)kraken:taxid\|(\d+)/) {
    print "TAXID\t$2\t$seq_id\n";
    next;
  }

  if ($seq_id !~ /(^|\|)gi\|(\d+)/) {
    die "$PROG: sequence ID $seq_id lacks GI number, aborting.\n";
  }
  print "$2\t$seq_id\n";
}
