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

# Create a file in a specified directory, then copy an
# existing file's contents into the new file.  Write name of
# new file to standard output.
#
# Thanks to everyone who wrote the mktemp program and couldn't be
# bothered to standardize the behavior.

use strict;
use warnings;
use File::Basename;
use File::Temp 'tempfile';
use Getopt::Std;

my $PROG = basename $0;
getopts('d:t:s:', \my %opts) or usage();
$opts{$_} or usage() for qw/d t s/;  # all switches mandatory
my ($directory, $template, $suffix) = @opts{qw/d t s/};
die "$PROG: '$directory' not a directory!\n" unless -d $directory;
die "$PROG: must specify a single filename\n" unless @ARGV == 1;

$suffix =~ s/^\.//;
my $old_filename = shift @ARGV;
open FILE, "<", $old_filename
  or die "$PROG: can't read $old_filename: $!\n";

my ($fh, $new_filename) = tempfile($template, DIR => $directory,
                                   UNLINK => 0, SUFFIX => ".$suffix");
# copy loop
while (<FILE>) {
  print {$fh} $_;
}
close FILE;
close $fh;

print "$new_filename\n";

sub usage {
  die "$PROG: <-d directory> <-t template> <-s suffix> <filename>\n";
}
