#!/usr/bin/perl

# Copyright 2013-2014, Derrick Wood <dwood@cs.umd.edu>
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

# Split a multi-fasta .ffn file into many single-fasta .fna files

use strict;
use warnings;
use File::Basename;

my $PROG = basename $0;

die "$PROG: must specify one filename!\n" if @ARGV != 1;

my $filename = shift;

my ($prefix, $extension);
if ($filename =~ /^(.*)\.ffn$/) {
  $prefix = $1;
}
else {
  die "$PROG: can't determine prefix!\n";
}
$extension = "fna";

my $ct = 0;
open MULTI, "<", $filename
  or die "$PROG: can't read $filename: $!\n";
while (<MULTI>) {
  if (/^>/) {
    close SINGLE if $ct;
    my $new_filename = "$prefix.$ct.$extension";
    open SINGLE, ">", $new_filename
      or die "$PROG: can't write to $new_filename: $!\n";
    $ct++;
  }
  print SINGLE;
}
close SINGLE;
close MULTI;
