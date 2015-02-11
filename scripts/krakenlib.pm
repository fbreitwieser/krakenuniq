package krakenlib;

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

# Common subroutines for other Kraken scripts

use strict;
use warnings;

# Input: the argument for a --db option (possibly undefined)
# Returns: the DB to use, taking KRAKEN_DEFAULT_DB and KRAKEN_PATH
#   into account.
sub find_db {
  my $supplied_db_prefix = shift;
  my $db_prefix;
  if (! defined $supplied_db_prefix) {
    if (! exists $ENV{"KRAKEN_DEFAULT_DB"}) {
      die "Must specify DB with either --db or \$KRAKEN_DEFAULT_DB\n";
    }
    $supplied_db_prefix = $ENV{"KRAKEN_DEFAULT_DB"};
  }
  my @db_path = (".");
  if (exists $ENV{"KRAKEN_DB_PATH"}) {
    my $path_str = $ENV{"KRAKEN_DB_PATH"};
    # Allow zero-length path to be current dir
    $path_str =~ s/^:/.:/;
    $path_str =~ s/:$/:./;
    $path_str =~ s/::/:.:/;

    @db_path = split /:/, $path_str;
  }
  
  # Use supplied DB if abs. or rel. path is given
  if ($supplied_db_prefix =~ m|/|) {
    $db_prefix = $supplied_db_prefix;
  }
  else {
    # Check all dirs in KRAKEN_DB_PATH
    for my $dir (@db_path) {
      my $checked_db = "$dir/$supplied_db_prefix";
      if (-e $checked_db && -d _) {
        $db_prefix = $checked_db;
        last;
      }
    }
    if (! defined $db_prefix) {
      my $printed_path = exists $ENV{"KRAKEN_DB_PATH"} ? qq|"$ENV{'KRAKEN_DB_PATH'}"| : "undefined";
      die "unable to find $supplied_db_prefix in \$KRAKEN_DB_PATH ($printed_path)\n";
    }
  }

  for my $file (qw/database.kdb database.idx/) {
    if (! -e "$db_prefix/$file") {
      die "database (\"$db_prefix\") does not contain necessary file $file\n";
    }
  }

  return $db_prefix;
}

1;
