#!/usr/bin/env perl

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

# Merges two files specified on the command line, print to stdout
# Designed to be called by kraken

use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $PROG = basename $0;

my $fasta_input = 0;
my $fastq_input = 0;
my $check_names = 0;
my $fh1;
my $fh2;

GetOptions(
  "check-names" => \$check_names
);

if (@ARGV != 2) {
  die "$PROG: must have exactly two filename arguments\n";
}
for my $file (@ARGV) {
  if (! -e $file) {
    die "$PROG: $file does not exist\n";
  } 
  if (! (-f $file || -p $file)) {
    die "$PROG: $file is not a regular file\n";
  }
}

#auto-detect compression types
my $ftype1 = determine_file_type($ARGV[0]);
my $ftype2 = determine_file_type($ARGV[1]);

if ($ftype1 eq "gzip") {
  open $fh1, "-|", "gunzip", "-c", $ARGV[0]
    or die "$PROG: can't open gunzip pipe with $ARGV[0]: $!\n";
} elsif ($ftype1 eq "bzip") {
  open $fh1, "-|", "bunzip2", "-c", $ARGV[0]
    or die "$PROG: can't open bunzip2 pipe with $ARGV[0]: $!\n";
} elsif($ftype1 eq "fastq" || $ftype1 eq "fasta") {
  open $fh1, "<", $ARGV[0]
    or die "$PROG: can't open $ARGV[0]: $!\n";
}

#auto-detect file type
my $line=<$fh1>;
if ($line =~ /^>/){
  $fasta_input=1;
} elsif ($line =~ /^@/){
  $fastq_input=1;
} else {
  die "Unknown file format for $ARGV[0]";
}
close($fh1);
#we assume that both files are fastq or both files are fasta

if ($ftype1 eq "gzip") {
  open $fh1, "-|", "gunzip", "-c", $ARGV[0]
    or die "$PROG: can't open gunzip pipe with $ARGV[0]: $!\n";
} elsif ($ftype1 eq "bzip") {
  open $fh1, "-|", "bunzip2", "-c", $ARGV[0]
    or die "$PROG: can't open bunzip2 pipe with $ARGV[0]: $!\n";
} elsif($ftype1 eq "fastq" || $ftype1 eq "fasta") {
  open $fh1, "<", $ARGV[0]
    or die "$PROG: can't open $ARGV[0]: $!\n";
}

if ($ftype2 eq "gzip") {
  open $fh2, "-|", "gunzip", "-c", $ARGV[1]
    or die "$PROG: can't open gunzip pipe with $ARGV[1]: $!\n";
} elsif ($ftype2 eq "bzip") {
  open $fh2, "-|", "bunzip2", "-c", $ARGV[1]
    or die "$PROG: can't open bunzip2 pipe with $ARGV[1]: $!\n";
} elsif($ftype2 eq "fastq" || $ftype2 eq "fasta") {
  open $fh2, "<", $ARGV[1]
    or die "$PROG: can't open $ARGV[1]: $!\n";
}

# read/merge/print loop
# make sure names match before merging
my ($seq1, $seq2);
while (defined($seq1 = read_sequence($fh1))) {
  $seq2 = read_sequence($fh2);
  if (! defined $seq2) {
    print STDERR "$PROG: mismatched sequence counts - file 1 has more reads\n
  Outputting the further reads unpaired\n";
    print_sequence($seq1);
    while (defined($seq1 = read_sequence($fh1))) {
      print_sequence($seq1);
    }
  }
  if ($check_names && $seq1->{id} ne $seq2->{id}) {
    die "$PROG: mismatched mate pair names ('$seq1->{id}' & '$seq2->{id}')\n";
  }
  print_merged_sequence($seq1, $seq2);
}
if (defined($seq2 = read_sequence($fh2))) {
  print STDERR "$PROG: mismatched sequence counts - file 2 has more reads\n
  Outputting the further reads unpaired\n";
  print_sequence($seq2);
  while (defined($seq2 = read_sequence($fh2))) {
    print_sequence($seq2);
  }

}
close $fh1;
close $fh2;

{
  my %buffers;  # Needed due to fasta formatting
  sub read_sequence {
    my $fh = shift;
    my $id;
    my $seq = "";
    if (! exists $buffers{$fh}) {
      $buffers{$fh} = <$fh>;
    }
    if (! defined $buffers{$fh}) {  # No more file
      return undef;
    }
    if ($fasta_input) {
      if ($buffers{$fh} =~ /^>(\S+)/) {
        $id = $1;
      }
      else {
        die "$PROG: malformed fasta file (line = '$buffers{$fh}')\n";
      }
      delete $buffers{$fh};
      while (<$fh>) {
        if (/^>/) {
          $buffers{$fh} = $_;
          last;
        }
        else {
          chomp;
          $seq .= $_;
        }
      }
    }
    elsif ($fastq_input) {
      if ($buffers{$fh} =~ /^@(\S+)/) {
        $id = $1;
      }
      else {
        if ($buffers{$fh} =~ /^$/) {
          return undef;  # Some fastq files end with blank line, allow it
        }
        die "$PROG: malformed fastq file (line = '$buffers{$fh}')\n";
      }
      delete $buffers{$fh};
      chomp($seq = <$fh>);
      scalar <$fh>;  # quality header
      scalar <$fh>;  # quality values
    }
    else {
      # should never get here
      die "$PROG: I have no idea what kind of input I'm reading!!!\n";
    }

    $id =~ s/[\/_.][12]$//;  # strip /1 (or .1, _1) or /2 to help comparison
    return { id => $id, seq => $seq };
  }
}

sub print_merged_sequence {
  my ($seq1, $seq2) = @_;
  print ">" . $seq1->{id} . "\n";
  print $seq1->{seq} . "N" . $seq2->{seq} . "\n";
}

sub print_sequence {
  my ($seq1) = @_;
  print ">" . $seq1->{id} . "\n";
  print $seq1->{seq} . "\n";
}

sub determine_compression_type {
    my ($file_path) = @_;
    open(my $fh, '<:raw', $file_path) or die "Could not open file '$file_path': $!";
    my $header = '';
    read($fh, $header, 3);
    close($fh);
    # Changing this line to just check for first character @
    if ($header =~ /^@/) {
        return 'fastq';
    }elsif ($header =~ /^>/) {
        return 'fasta';
    }elsif ($header eq 'BZh') {
        return 'bzip';
    }elsif ($header eq "\x1f\x8b\x08") {
        return 'gzip';
    }else {
        die "Unknown input file type for $file_path";
    }
}
