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

# Checks each sequence header to ensure it has a GI number to
# enable taxonomic ID lookup later.  Also has some (very basic)
# FASTA-format checking.

use strict;
use warnings;
use File::Basename;

my $PROG = basename $0;

die "$PROG: must specify one filename!\n" if @ARGV != 1;

my $filename = shift;

open FASTA, "<", $filename
  or die "$PROG: can't open $filename: $!\n";
my $seq_ct = 0;
my $errors = 0;
while (<FASTA>) {
  next unless /^>/;
  $seq_ct++;
  if (! /^>(\S+)/) {
    $errors++;
    warn "file $filename, line $. lacks sequence ID\n";
  }
  if ($1 !~ /(^|\|)(gi|kraken:taxid)\|(\d+)/) {
    $errors++;
    warn "file $filename, line $.: sequence ID lacks GI number\n";
  }
}
close FASTA;

if ($errors) {
  exit 1;
}
