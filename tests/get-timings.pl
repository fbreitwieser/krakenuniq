#!/usr/bin/env perl

use strict;
use warnings;

my $time = "wall clock";
my $mem = "Maximum resident set size";
my $sep = ";";

sub uutime {
		my $wt_h = `grep '$time' $_[0] | sed 's/.* //'`; chomp $wt_h;
		return $wt_h;
}

my $K="kraken";

print join($sep,qw/sample mbpm1 mbpm2 time1 time2 time3 mem1 mem2/),"\n";
for my $hll_t (@ARGV) {
		my $bn = `basename $hll_t .krakenhll.timing.log`; chomp $bn;
		die "$hll_t should be an hll timing log" if $bn == $hll_t;
		print "$bn$sep";
		(my $hll = $hll_t) =~ s/.timing//;
		(my $k = $hll) =~ s/krakenhll/$K/;
		(my $k_t = $hll_t) =~ s/krakenhll/$K/;
		(my $kr_t = $hll_t) =~ s/krakenhll/$K-report/;

		my $hm = `grep -o '[0-9\.]* Mbp/m' $hll | sed 's/ .*//'`; chomp $hm;
		print $hm,$sep;
		my $km = `grep -o '[0-9\.]* Mbp/m' $k | sed 's/ .*//'`; chomp $km;
		print $km,$sep;

		print uutime($hll_t),$sep;
		print uutime($k_t),$sep;
		print uutime($kr_t),$sep;

		my $mem_h = `grep '$mem' $hll_t | sed 's/.* //'`; chomp $mem_h;
		print $mem_h,$sep;
		my $mem_k = `grep '$mem' $k_t | sed 's/.* //'`; chomp $mem_k;
		print $mem_k,"\n";
}
