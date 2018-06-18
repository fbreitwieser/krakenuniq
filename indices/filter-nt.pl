#! /usr/bin/env perl
#
# Sort nt file sequences according to their taxonomy ID
# Uses the new mapping file format available at
# ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/
#
# Author fbreitwieser <fbreitwieser@sherman>
#
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;

my $nodes_file;
my $names_file;
my $taxonlist_file;
my $nt_file;
my $new_map_file;
my $ac_wo_mapping_file;
my $opt_help;

my $USAGE = 
"USAGE: ".basename($0)." [OPTIONS] -n nodes.dmp -N names.dmp -t taxonlist.csv -f nt.fna <mapping file> [<mapping file>*]\n

OPTIONS:
  -m str      Output mappings that are present in sequence file to file str
  -a str      Output ACs w/o mapping to file str
  -h          This message
";

GetOptions(
    "nodes=s" => \$nodes_file,
    "names=s" => \$names_file,
    "t|taxonlist=s" => \$taxonlist_file,
    "f|fasta=s" => \$nt_file,
    "m|map=s" => \$new_map_file,
    "a=s" => \$ac_wo_mapping_file,
    "h|help" => \$opt_help) or die "Error in command line arguments";

scalar(@ARGV) >= 1 or die $USAGE;
if (defined $opt_help) {
  print STDERR $USAGE;
  exit(0);
}

my @ac_taxid_files = @ARGV;

my %ac_to_pos;
my %taxid_to_ac;
my %ac_to_taxid;


my %parent_map;
my %child_map;
-f $nodes_file or die "Nodes files $nodes_file is not a file";
-f $names_file or die "Names files $names_file is not a file";
print STDERR "Reading taxonomy tree from $nodes_file ... ";
open(my $N, "<", $nodes_file) or die $!;
while (<$N>) {
  my ($taxid, $parent_taxid) = split /\t\|\t/;
  $parent_map{$taxid} = $parent_taxid;
  push @{$child_map{$parent_taxid}}, $taxid if $taxid ne $parent_taxid;
  die "parent??" unless defined $parent_taxid;
}
close($N);
print STDERR "Got ",(scalar keys %parent_map), " nodes.\n";

my %name_map;
print STDERR "Readin names from $names_file ... ";
open(my $M, "<", $names_file) or die $!;
while (<$M>) {
  next unless /scientific name\t|$/;
  my ($taxID, $name) = split /\t\|\t/;
  push @{$name_map{$name}}, $taxID;
}
close($M);
print STDERR "Done\n";

my %accepted_taxa;

sub recurse_tree_set_taxa {
  my ($taxid, $invert) = @_;
  if ($invert) {
    delete $accepted_taxa{$taxid};
  } else { 
    $accepted_taxa{$taxid} = 1;
  }
  return unless defined $child_map{$taxid};
  foreach my $child_taxid (@{$child_map{$taxid}}) {
    recurse_tree_set_taxa($child_taxid, $invert);
  }
}

my @accepted_taxon_list;
my @rejected_taxon_list;
print STDERR "Reading taxon list from $taxonlist_file ... ";
open (my $TL, "<", $taxonlist_file) or die $!;
while (<$TL>) {
  next if /^#/;
  chomp;
  my ($taxid, $name) = split /\t/;
  if ($taxid == 0) {
    $name =~ s/#.*$//;
    $name =~ s/\s+$//;
    die "Got no taxIDs for name '$name'!" unless defined $name_map{$name};
    push @rejected_taxon_list, @{$name_map{$name}};
  } elsif ($taxid == abs($taxid)) {
    push @accepted_taxon_list, $taxid;
    recurse_tree_set_taxa($taxid, 0);
  } else {
    push @rejected_taxon_list, $taxid;
  }
}
close($TL);
foreach my $taxid (@rejected_taxon_list) {
  recurse_tree_set_taxa($taxid, 1);
}
print STDERR "Got ",(scalar(@accepted_taxon_list))," taxon listed with ",(scalar keys %accepted_taxa)," children.\n";

print STDERR "Reading headers from $nt_file ... ";
open(my $NT, "<", $nt_file) or die $!;
while (<$NT>) {
  # get the headers with (!) the version number
  if (/(^>([^ ]*).*)/) {
    # record the position of this AC
    $ac_to_pos{$2} = [tell($NT),$1];
  }
}
print STDERR "found ", scalar(keys %ac_to_pos), " ACs\n";

my $n_accepted_acs = 0;
my $n_declined_acs = 0;

foreach my $ac_taxid_file (@ac_taxid_files) {
print STDERR "Reading ac to taxid mapping from $ac_taxid_file ...\n";
  my $FP1;
  if ($ac_taxid_file =~ /.gz$/) {
    open($FP1, "-|", "gunzip -c '$ac_taxid_file'") or die $!;
  } else {
    open($FP1, "<", $ac_taxid_file) or die $!;
  }

  # format: accession <tab> accession.version <tab> taxid <tab> gi
  # currently we look for a mapping with the version number
  while ( <$FP1> ) {
    my (undef, $ac, $taxid) = split;
    next unless defined $taxid;
    if (defined $accepted_taxa{$taxid} && defined( $ac_to_pos{ $ac } ) ) {
      push @{ $taxid_to_ac{ $taxid } }, $ac;
      $ac_to_taxid{ $ac } = $taxid;
    }
  }
  close $FP1;
}
print STDERR "Got taxonomy mappings for ", scalar(keys %ac_to_taxid), " ACs\n";
if (defined $ac_wo_mapping_file && scalar(keys %ac_to_taxid) < scalar(keys %ac_to_pos)) {
  print STDERR "Writing ACs without taxonomy mapping to $ac_wo_mapping_file\n";
  open(my $FP2, ">", $ac_wo_mapping_file) or die $!;
  foreach my $ac (keys %ac_to_pos) {
    next unless defined $ac_to_taxid{$ac};
    print $FP2 $ac, "\n";
  }
  close($FP2);
}

if (defined $new_map_file) {
print STDERR "Writing taxonomy ID mapping to $new_map_file\n";
open(my $FP3, ">", $new_map_file) or die $!;
foreach my $ac (keys %ac_to_taxid) {
  print $FP3 $ac,"\t",$ac_to_taxid{$ac},"\n";
}
close($FP3);
}


print STDERR "Outputting sorted FASTA ...\n";
foreach my $taxid (sort {$a <=> $b} keys %taxid_to_ac) {
  my @acs = @{$taxid_to_ac{$taxid}};
  my @sorted_acs = sort { $ac_to_pos{$a}->[0] <=> $ac_to_pos{$b}->[0] } @acs;
  foreach (@sorted_acs) {
    print $ac_to_pos{$_}->[1],"\n";
    seek($NT, $ac_to_pos{$_}->[0], 0);
    while (<$NT>) {
      last if (/^>/);
      print $_;
    }
  }
}
close $NT;
