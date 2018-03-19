#!/usr/bin/env perl
#vim: et:ts=2:sw=2

# krakenhll-download.pl - based on centrifuge-download
# (c) Florian Breitwieser, 2017-2018
# licensed under GPL-3

use strict;
use warnings;
use File::Basename;
use File::Fetch;
use File::Copy;
use File::Path qw/make_path remove_tree/;
use IO::Uncompress::Gunzip qw/gunzip $GunzipError/;
use autodie;
use Term::ANSIColor;
use Getopt::Long;
use LWP::UserAgent;
use List::Util qw/min/;


sub download_taxonomy(@);
sub download_contaminats(@);
sub download(@);
sub print_header_lines(@);
sub print_header_line_ac(@);
sub download_domain(@);
sub download_viral_neighbors(@);
sub download_ncbi_search(@);
sub get_sorted_maps(@);

my $FTP="ftp://ftp.ncbi.nih.gov";
my @ALL_GENOMES=qw/bacteria viral archaea fungi protozoa invertebrate plant vertebrate_mammalian vertebrate_other/;
my @ALL_DATABASES=qw/refseq genbank taxonomy contaminants/;
my @ALL_ASSEMBLY_LEVELS=qw/Complete\ Genome Chromosome Scaffold Contig/;
my @SMALL_GENOMES=qw/mitochondrion plasmid plastid/;

## Option parsing
my $DATABASE="refseq";
my $ASSEMBLY_LEVEL="Complete Genome";
my $REFSEQ_CATEGORY;
my $TAXID;

my $BASE_DIR;
my $DB_DIR;
my $N_PROC=5;
my $MAP_DIV="nucl_gb";
my $CHANGE_HEADER=0;
my $DO_DUST=0;
my $FILTER_UNPLACED=0;
my $SEARCH_TERM;
my $VERBOSE=0;
my $OVERWRITE_FILES=0;
my $DOMAINS;
my $DL_MOD_RSYNC;
my $n_children = 0;
my $INCLUDE_VIRAL_NEIGHBORS = 0;
my %taxid_name_map;
my @pids;

my $downloaded_viral_refseq=0;
my $FNA_FILES="genomic";

my $USAGE="\n".basename($0).
" [<options>] <database> <database>*

ARGUMENT
 <database>  Possible databases
              'contaminants'    downloads contaminant sequences from UniVec and EmVec
              'taxonomy'        downloads the NCBI taxonomy mappings
              'genbank'         downloads GenBank genomes (see parameters -a and -d)
              'refseq'          downloads RefSeq genomes (see parameters -a and -d)
			  'nucleotide'      specify --term to download nucleotide sequences retrieved from a NCBI search.
              'viral-neighbors' downloads sequences listed as neighbors to RefSeq genomes in the NCBI viral genome resource.
                                https://www.ncbi.nlm.nih.gov/genome/viruses/

    for refseq and genbank, the assemblies can be filtered:
              'refseq/DOMAIN'   downloads RefSeq genomes of DOMAIN (e.g. 'refseq/bacteria'). 
			                    possible values: @ALL_GENOMES.
              'refseq/DOMAIN/ASS_LEVEL'
                                specifiy assembly level (e.g. 'refseq/fungi/Chromosome' or 'refseq/viral/Any')
              'refseq/DOMAIN/ASS_LEVEL/column1=value1(/column2=value2)*' 
                                require specific column values in assembly_summary.txt.
                                 e.g. 'refseq/vertebrate_mammalian/Any/species_taxid=9606' <- download human assemblies
                                   or 'genbank/bacteria/Any/species_taxid=9606' <- downloads Genbank assemblies not in RefSeq

              
  By default, 'refseq' and 'genbank' download the domains and assembly levels specied
  with other parameters, however
                   - refseq and genbank can be proceeded by '/DOMAIN' or '/DOMAIN/ASS_LEVEL', e.g.
                     - refseq/archaea, refseq/viral/Any, or genbank/bacteria
                     - if ASS_LEVEL is not given, the default is used

COMMON OPTIONS
 -o <directory>     Folder to which the files are downloaded. Default: '.'
 --db <directory>   Alternative to -o: Download to <directory>/{library,taxonomy}.
 --threads <# of threads>  Number of processes when downloading (uses xargs). Default: '$N_PROC'
 --rsync, -R        Download using rsync.
 --overwrite        Redownload and overwrite files with the same name.
 --verbose          Be verbose.

WHEN USING DATABASE nucleotide:
 --term search_term    Download all sequences returned from a NCBI nucleotide search.
 --mappings divisions  Try mapping accession IDs using the mapping files for the specified divisions (comma-separated).
                       Default: $MAP_DIV. Possible values: nucl_est, nucl_gb, nucl_gss, nucl_wgs.
					   Downloaded from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/.

WHEN USING DATABASE refseq OR genbank:
 --fna <seq types>  Comma-separated list of sequence types, including genomic, rna, rna_from_genomic, cds_from_genomic. Default: $FNA_FILES.
                    See the assembly project FTP site for available sequences
 -u                 Filter unplaced sequences.
 -l                 Modify sequence header to include taxonomy ID for Kraken (i.e. add '>kraken:taxid|TAXID' to each sequence).
 --dust, -D         Mask low-complexity regions using dustmasker.
 --include-viral-neighbors  Include viral neighbors. Deprecated. Add 'viral-neighbors' to the arguments, instead.
";

# arguments: $OPTFIND (current index), $OPTARG (argument for option), $OPTERR (bash-specific)
Getopt::Long::Configure('no_auto_abbrev','pass_through');
GetOptions(
  "output|o=s"  =>\$BASE_DIR,
  "db=s" => \$DB_DIR,
  "threads|P=i" =>\$N_PROC,
  "domain|d=s"  => \$DOMAINS,
  "assembly-level|a=s" => \$ASSEMBLY_LEVEL,
  "category|c=s" => \$REFSEQ_CATEGORY,
  "taxonomy-id|t=s" => \$TAXID,
  "fna=s" => \$FNA_FILES,
  "rsync|R" => \$DL_MOD_RSYNC,
  "include-viral-neighbors" => \$INCLUDE_VIRAL_NEIGHBORS,
  "filter-unplaced|u" => \$FILTER_UNPLACED,
  "term=s" => \$SEARCH_TERM,
  "mappings=s" => \$MAP_DIV,
  "dust|D" => \$DO_DUST,
  "change-header|l" => \$CHANGE_HEADER,
  "force" => \$OVERWRITE_FILES,
  "verbose|v" => \$VERBOSE) or die "Error in command line arguments";

if (defined $BASE_DIR && defined $DB_DIR) {
  print "Define either --db or -o, not both!";
  exit 1;
}

my $ua = LWP::UserAgent->new( ssl_opts => { verify_hostname => 0 } );

# Use current directory as base directory
$BASE_DIR = "." unless defined $DB_DIR || defined $BASE_DIR;

# If DB directory is defined, use that as base directory
#  -- kept -o and --db options to allow the use of either Kraken and Centrifuge type command line
my $add_dir = defined $DB_DIR;
$BASE_DIR = $DB_DIR if defined $DB_DIR;
sub get_dir {
  my ($dir, $name) = @_;
  my $dir1 = $add_dir? "$dir/$name" : $dir;
  make_path $dir1;
  return $dir1;
}

my %select_taxonomy_ids;
if (defined $TAXID) {
  %select_taxonomy_ids = map { $_ => 1 } split(/,/, $TAXID);
}

if (!defined $ARGV[0] && !$INCLUDE_VIRAL_NEIGHBORS) {
  print STDERR $USAGE;
  exit 1;
}

if ($INCLUDE_VIRAL_NEIGHBORS) {
  my $nbr_lib_dir = $add_dir? "$BASE_DIR/library/viral/Neighbors" : "$BASE_DIR/viral/Neighbors";
  download_viral_neighbors_url($nbr_lib_dir);
}

foreach my $DATABASE (@ARGV) {
  if ( $DATABASE eq "taxonomy" ) { 
    download_taxonomy(get_dir($BASE_DIR,"taxonomy"));
  } elsif ( $DATABASE eq "contaminants" ) { 
    download_contaminats(get_dir($BASE_DIR,"library/contaminants"));
  } elsif ( $DATABASE =~ /^refseq/ || $DATABASE =~ /^genbank/ ) {
    my ($db, $domains, $assembly_levels, @additional_filters) = split(/\//, $DATABASE);
    $domains = $DOMAINS unless defined $domains;
    $assembly_levels = $ASSEMBLY_LEVEL unless defined $assembly_levels;

    foreach my $domain (split(/,/, $domains)) {
      my $lib_dir = $add_dir? "$BASE_DIR/library/$domain" : "$BASE_DIR/$domain";
      foreach my $assembly_level (split(/,/, $assembly_levels)) {
          download_domain($db, $lib_dir, $domain, $assembly_level, @additional_filters);
      }
    }
  } elsif ($DATABASE eq'viral-neighbors') {
    my $nbr_lib_dir = $add_dir? "$BASE_DIR/library/viral/Neighbors" : "$BASE_DIR/viral/Neighbors";
    download_viral_neighbors($nbr_lib_dir);
  } elsif ($DATABASE eq 'nucleotide') {
    if (!defined $SEARCH_TERM) {
	  print STDERR "Please define a search term with --term when using database nucleotide. \n".
	               " E.g. --term 'srcdb_refseq[prop] plasmid[title] complete sequence[title] txid2[organism]'\n.";
	  exit 1;
	}
    my $nbr_dir = $add_dir? "$BASE_DIR/library/nucleotide/" : "$BASE_DIR/nucleotide/";
    download_ncbi_search($nbr_dir, $SEARCH_TERM, get_sorted_maps(split /,/, $MAP_DIV));
  } else {
    print STDERR "Unknown database $DATABASE. \n";
    print STDERR $USAGE;
    exit 1;
  }
}





#########################################################
## Functions

sub download(@) {
  my ($url, $file, $opts) = @_;
  my $gunzipped_filename;

  if (defined $opts && defined $opts->{'gunzip'}) {
    ($gunzipped_filename = $file) =~ s/.gz$//; 
  }

  if (defined $opts && defined $opts->{'verbose'} && !$VERBOSE) {
    print STDERR "$file ";
  }

  if (!$OVERWRITE_FILES && (( defined $gunzipped_filename && -s $gunzipped_filename) || (!defined $gunzipped_filename && -s $file))) {
    print STDERR "Not fetching $url - file $file exists.\n" if $VERBOSE;
  	if (defined $opts && defined $opts->{'verbose'} && !$VERBOSE) {
      print colored("check\n", "green");
	}
    return 1;
  }
  if (defined $opts && defined $opts->{'verbose'} && !$VERBOSE) {
    print STDERR ": downloading ...";
  }

  if ($url =~ /^http/) {
    print STDERR "Fetching $url to $file ..." if $VERBOSE;
    if (!-d dirname($file)) {
      make_path(dirname($file));
    }
    my $response = $ua->get($url, ':content_file' => $file);
    if (!$response->is_success) {
      print STDERR "\nFAIL: Error downloading $url ($response->status_line) - tying curl\n";
	  system("curl '$url' -o $file") == 0 or die "Error curling $url.\n";
    } else {
      print STDERR "SUCCESS\n" if $VERBOSE;
    }
  } else {
    if ( $DL_MOD_RSYNC && $url =~ /^ftp/ ) {
     $url =~ s/^ftp/rsync/;
    }
    print STDERR "Fetching $url to $file ..." if $VERBOSE;

    my $ff = File::Fetch->new(uri=>"$url");
    my $where = $ff->fetch(to=> dirname($file)) or die $ff->error;
    move($where, $file);
  }

  if (defined $gunzipped_filename && $gunzipped_filename ne $file) {
    if (defined $opts && defined $opts->{'verbose'} && !$VERBOSE) {
	  print STDERR " gunzipping ...";
	}
    print STDERR " GUNZIPPING" if $VERBOSE;
    gunzip $file => $gunzipped_filename or die "gunzip failed: $GunzipError";
    unlink $file;
    $file = $gunzipped_filename;
  }
  if (-s $file || (defined $gunzipped_filename && -s $gunzipped_filename)) {
    print STDERR " SUCCESS\n" if $VERBOSE;
    if (defined $opts && defined $opts->{'verbose'} && !$VERBOSE) {
      print STDERR " done.\n";
	}
  } else {
    print STDERR "failed.\n";
  }

  #my $where = $ff->fetch(to=> dirname($file)) or die "\n$ff->error for $url!";
  return -s $file;
}

sub start_fork() {
  my $pid;
  return if $N_PROC <= 1;
  if ($n_children == $N_PROC) {
    $pid = wait();
    --$n_children;
  }
  if (defined($pid = fork())) {
    if ($pid) {
      ++$n_children;
      #print STDERR "Parent: forked child $pid\n";
      push @pids, $pid;
    } 
  } else {
    print STDERR "ERROR: Failed to fork\n";
  }
  return $pid;
}

sub wait_children() {
  foreach my $pid (@pids) {
    waitpid $pid, 0;
  }
  @pids = ();
  $n_children = 0;
}

sub end_fork() {
  exit() unless $N_PROC <= 1;
}

sub download_ncbi_search(@) {
  my ($nbr_dir, $term, @sorted_map_files) = @_;
  $term =~ s/ /+/g;

  initialize_name_map();

  my $esearch_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi";
  my $complete_url="$esearch_url?db=nuccore&usehistory=y&retmax=1&retmode=json&term=$term";
  my $esearch_file = "$nbr_dir/esearch_neighbor_res.json";
  download($complete_url, $esearch_file);

  my $n_res=`grep -m1 '"count":' $esearch_file | sed -e 's/.*: "//' -e 's/",.*//'`;
  chomp $n_res;
  print STDERR "Downloading $n_res sequences into $nbr_dir.\n";
  return if (!defined $n_res || $n_res eq 0);

  # Step 2: Download FASTAs, 10k at a time
  my $url_params=`cat $esearch_file | grep -e 'querykey' -e 'webenv' | sed -e 's/^ *"querykey": "/query_key=/' -e 's/^ *"webenv": "/webenv=/' -e 's/",//' | paste -sd\\&`;
  chomp $url_params;
  die "error getting viral neighbors with url $complete_url" if (!defined $url_params || $url_params eq "");
  my $retstart = 0;
  my $retmax = 10000;
  my @all_fas = ();
  while ($retstart < $n_res) {

    my $curr_retmax = min($n_res, $retstart+$retmax);
    print STDERR "\r  Downloading sequences ".($retstart+1)." to $curr_retmax of $n_res ..." unless $VERBOSE;

    #start_fork() and next;
	my $part_file = "$nbr_dir/dl_part_".($retstart+1)."-to-$curr_retmax.fna.tmp";
    download("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&$url_params&rettype=fasta&retstart=$retstart&retmax=$curr_retmax", $part_file);

    my $FA_HANDLE;
    my $fa_handle_open = 0;
    open (my $F, "<", $part_file) or die "Couln't open downloaded file";
    while (my $line = <$F>) {
      next if $line eq "\n";
      if ($line =~ /^>/) {
        (my $nbr_ac = $line) =~ s/.//;
        $nbr_ac =~ s/ .*//;
        chomp $nbr_ac;
   
		my $taxid;
		foreach (@sorted_map_files) {
          $taxid = `look $nbr_ac $_ | cut -f 2`;
		  chomp $taxid;
		  last if defined $taxid;
        }

        if (!defined $taxid || !defined $taxid_name_map{$taxid}) {
          my $res = `curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$nbr_ac&rettype=fasta&retmode=xml" | head -n 12  | egrep '<TSeq_taxid>|<TSeq_orgname>'  | sed -e 's#</.*>##' -e 's#.*<.*>##' | paste -sd\$'\\t'`;
          chomp $res;
          ($taxid) = split /\t/, $res;
        }

        my $name = $taxid_name_map{$taxid};
        if (!defined $taxid || !defined $name) {
          print STDERR "\nNo mapping for viral neighbor $nbr_ac!\n";
          next;
        }
        (my $name1 = $name) =~ s/[^a-zA-Z0-9_]/_/g;
        $name1 =~ s/__/_/g;
        system("mkdir -p $nbr_dir/$name1-tax$taxid");
        if ($fa_handle_open) {
          close($FA_HANDLE);
        } else {
          $fa_handle_open = 1;
        }
        my $file = "$nbr_dir/$name1-tax$taxid/$nbr_ac.fna";
        open($FA_HANDLE, ">", $file) or die "Couln't open file";
        print_header_line_ac($file, $nbr_ac, $taxid);
      }
      print $FA_HANDLE $line;
      #end_fork();
    }
    close($F);
    $retstart += $curr_retmax;
  }
  print STDERR "\n";
  wait_children();
#  $pm->wait_all_children();
}

sub get_sorted_maps(@) {
  my (@parts) = @_;
  my @res;
  my $dir = get_dir($BASE_DIR,"taxonomy");
  foreach (@parts) {
    my $url = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/$_.accession2taxid.gz";
    download($url, "$dir/$_.accession2taxid.gz", { verbose => 1 });

    my $sorted_map_f = "$dir/$_.accession2taxid.sorted";
    if (!-s $sorted_map_f) {
      print STDERR "Sorting mapping file ...\n";
      my $sort_cmd = system("sort --help | grep -q parallel") == 0? "sort --parallel $N_PROC" : "sort";
      system("gunzip -c $dir/$_.accession2taxid.gz | cut -f 2,3 | $sort_cmd -T $dir > $sorted_map_f") == 0 or die "Error sorting: $!";
    }
	push @res, $sorted_map_f;
  }
  return @res;
}

sub initialize_name_map() {
  return if %taxid_name_map;
  print STDERR "Reading names file ...\n";
  my $dir = get_dir($BASE_DIR,"taxonomy");
  my $names_file = "$dir/names.dmp";
  if (!-f $names_file) {
    download_taxonomy($dir);
  }
  open (my $N, "<", $names_file);
  while (<$N>) {
    next unless /scientific name/;
    my ($taxid, $name) = split /\t\|\t/;
    $taxid_name_map{$taxid} = $name;
  }
  close($N);
}

sub download_viral_neighbors(@) {
  my ($nbr_dir) = @_;
  print STDERR "Downloading viral neighbors.\n";
  my $term="viruses[organism] not cellular organisms[orgn] not wgs[prop] not gbdiv syn[prop] and (nuccore genome samespecies[filter] and (complete[title]) not unverified[title])";
  download_ncbi_search($nbr_dir, $term, get_sorted_maps("nucl_gb"));
}

sub print_header_lines(@) {
  my ($file, $taxid, $name) = @_;
  return if -s "$file.map";
  if (! -f $file) {
    print STDERR "ERROR: $file does not exist - should not happen!?\n";
  }
  print STDERR "Writing mapping file for $file [$taxid, $name]\n" if $VERBOSE;
  $taxid = "$taxid\t$name" if defined $name;
  #`grep '^>' "$file" | sed -e 's/.//' -e 's/\\( .*\\|\$\\)/\t$taxid/' > "$file.map"`;
  open (my $F, ">", "$file.map");
  open (my $G, "<", $file);
  while (<$G>) {
   next unless /^>([^ ]*)/;
    my $ac = $1;
    print $F "$ac\t$taxid\n";
  }
  close($G);
  close($F);
}


sub print_header_line_ac(@) {
  my ($file, $ac, $taxid, $name) = @_;
  return if -s "$file.map";
  print STDERR "Making map file for $file\n" if ($VERBOSE);
  open (my $F, ">", "$file.map");
  if (defined $name) {
    print $F "$ac\t$taxid\t$name\n";
  } else {
    print $F "$ac\t$taxid\n";
  }
  close($F);
}


sub download_contaminats(@) {
  my ($CONTAMINANT_DIR) = @_;
  print STDERR "Downloading contaminant databases ... \n";
  my $CONTAMINANT_TAXID=32630;
  make_path $CONTAMINANT_DIR;

  # download UniVec and EmVec database
  download("ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec","$CONTAMINANT_DIR/UniVec.fna");
  download("ftp://ftp.ebi.ac.uk/pub/databases/emvec/emvec.dat.gz","$CONTAMINANT_DIR/emvec.dat.gz", { gunzip => 1 });

  open(my $E1, "<", "$CONTAMINANT_DIR/emvec.dat");
  open(my $E2, ">", "$CONTAMINANT_DIR/EmVec.fna");

  my ($ac,$de);
  my $in_seq = 0;
  while(<$E1>) {
    if (/^AC\s+(.*)/) {
      $ac = $1;
      $ac =~ s/;$//;
    } elsif (/^DE\s+(.*)/) {
      $de = $1;
   } elsif (/^SQ/) {
      $in_seq = 1;
      print $E2 ">$ac $de\n";
    } elsif ($in_seq) {
      if (/^\s+[agct]/) {
        s/\s+[0-9]+$//;
       s/ //g;
       print $E2 $_;
      } else {
        $in_seq = 0;
      }
    }
  }
  close($E2);
  close($E1);
  unlink("$CONTAMINANT_DIR/emvec.dat");
 
  if ( $CHANGE_HEADER ) {
    system("sed -i 's/^>/>taxid|$CONTAMINANT_TAXID /' $CONTAMINANT_DIR/UniVec.fna");
    system("sed -i 's/^>/>taxid|$CONTAMINANT_TAXID /' $CONTAMINANT_DIR/EmVec.fna");
  } else {
    print_header_lines("$CONTAMINANT_DIR/UniVec.fna", $CONTAMINANT_TAXID, "UniVec");
    print_header_lines("$CONTAMINANT_DIR/EmVec.fna", $CONTAMINANT_TAXID, "EmVec");
  }
}

sub download_taxonomy(@) {
  my ($dir) = @_;
  print STDERR "Downloading NCBI taxonomy ... \n";
  make_path $dir;

  download("$FTP/pub/taxonomy/taxdump.tar.gz", "$dir/taxdump.tar.gz");
  system("tar -C $dir -zxvf $dir/taxdump.tar.gz nodes.dmp names.dmp 1>&2");
  system("date > $dir/timestamp");
}

sub download_domain(@) {
  my ($DATABASE, $domain_dir, $domain, $_assembly_level, @additional_filters) = @_;
  $_assembly_level =~ s/ /_/g;
  print STDERR "Downloading assembly summary file for $domain genomes, and filtering to assembly level $_assembly_level";
  print STDERR (@additional_filters? " and additional filters @additional_filters.\n" : ".\n");
  die unless defined $domain_dir && defined $domain;
  if (-d $domain_dir) {
    print STDERR "WARNING: $domain_dir already exists - potentially overwriting files.\n";
  } else {
    make_path $domain_dir;
  }
  my $ass_file = "$domain_dir/assembly_summary.txt";
  my $ass_file_filtered = "$domain_dir/assembly_summary_filtered.txt";
  my $n_genomes = 0;
  download("ftp://ftp.ncbi.nlm.nih.gov/genomes/$DATABASE/$domain/assembly_summary.txt", $ass_file) or die "Could not download assembly summary file!";

  my $is_viral_refseq =1 if $domain eq "viral" && $DATABASE eq "refseq";

  my %cols = (
    assembly_accession=>0,
    bioproject=>1,
    biosample=>2,
    wgs_master=>3,
    refseq_category=>4,
    taxid=>5,
    species_taxid=>6,
    organism_name=>7,
    infraspecific_name=>8,
    isolate=>9,
    version_status=>10,
    assembly_level=>11,
    release_type=>12,
    genome_rep=>13,
    seq_rel_date=>14,
    asm_name=>15,
    submitter=>16,
    gbrs_paired_asm=>17,
    paired_asm_comp=>18,
    ftp_path=>19,
    excluded_from_refseq=>20,
    relation_to_type_material=>21
  );

  my %seen_acs;
  my @genomes_to_dl;
  open(my $A1, "<", $ass_file);
  open(my $A2, ">", $ass_file_filtered);
  while (<$A1>) {
    next if /^#/;
    my @fields = split /\t/;

    $fields[$cols{"assembly_level"}] =~ s/ /_/g;
    next unless $fields[$cols{"version_status"}] eq "latest";
    next if ($_assembly_level ne "Any" && $fields[$cols{"assembly_level"}] ne $_assembly_level);
    next if (defined $REFSEQ_CATEGORY && $fields[$cols{"refseq_category"}] ne $REFSEQ_CATEGORY);

    ## Kick out duplicates
    next if (defined $seen_acs{$fields[$cols{"assembly_accession"}]});
    $seen_acs{$fields[$cols{"assembly_accession"}]} = 1;
    
    my $keep_it = 1;
    foreach (@additional_filters) {
        my ($k, $v) = split(/=/);
        die "$k is not an available column filter. Available: ".(join(", ", sort keys %cols)) if (!defined $cols{$k});
        $keep_it = 0 if $fields[$cols{$k}] ne $v; 
    }
    next unless $keep_it;

    print $A2 $_;
    ++ $n_genomes;
    push @genomes_to_dl, [
            $fields[$cols{"ftp_path"}], $fields[$cols{"taxid"}], $fields[$cols{"organism_name"}], 
            $fields[$cols{"infraspecific_name"}], $fields[$cols{"assembly_accession"}], $fields[$cols{"assembly_level"}]];
  }
  close $A2;
  close $A1;

  my $downloaded_files = 0;
  my $existing_files = 0;

  my $i = 0;
  foreach my $g (@genomes_to_dl) {
    my ($ftp_path, $taxid, $organism_name, $infraspecific_name, $assembly_accession, $assembly_level) = @$g;
    ++$i;
    #print STDERR "\r                                                                               " unless $VERBOSE;
    print STDERR "\r Downloading $domain genomes:  $i/$n_genomes ... " unless $VERBOSE;

    if (defined $infraspecific_name) {
        (my $i1 = $infraspecific_name) =~ s/strain=//;
        $organism_name .= " $infraspecific_name" unless $organism_name =~ /\Q$i1\E/ || $i1 eq "";
    }


    my $bname = basename($ftp_path);
    ( my $organism_name1 = $organism_name ) =~ s/[^a-zA-Z0-9_]/_/g;
    $organism_name1 = substr($organism_name1, 0, 100);
    $organism_name1 =~ s/__/_/g;
    $organism_name1 =~ s/_$//;
    my $bname1 = "${organism_name1}-tax${taxid}-${bname}";
    
    $assembly_level =~ s/ /_/g;
    my $download_dir = "$domain_dir/$assembly_level";
    my $nbr_download_dir = "$domain_dir/Neighbors";
    foreach my $ext (split(/,/, $FNA_FILES)) {
      my $full_ftp_path = "$ftp_path/${bname}_${ext}.fna.gz";
      my $bfname = $bname1."_".$ext;
      my $nbr_fname = $bfname."-neighbors.fna";
      my $fname = $bfname.".fna";
      my $fullfname1 = $DO_DUST? "$download_dir/${bfname}_dustmasked.fna" : "$download_dir/$fname";

      if (!$OVERWRITE_FILES && -s $fullfname1) {
        print STDERR "$download_dir/$fname exists - not downloading.. \n" if $VERBOSE;
        ++$existing_files;
      } else {
        ++$downloaded_files;
        start_fork() and next;
        download($full_ftp_path, "$download_dir/$fname.gz", { gunzip => 1 });
        end_fork();
      }

      if ($CHANGE_HEADER) {
        system("sed -i 's/^>/>kraken:taxid|$taxid /' '$download_dir/$fname'");
      }
      if ($FILTER_UNPLACED) {
        die("Not implemented");
      }

      if ($DO_DUST || !-s $fullfname1) {
        start_fork() and next;
        ## TODO: Consider hard-masking only low-complexity stretches with 10 or more bps
        system("dustmasker -infmt fasta -in '$download_dir/$fname' -level 20 -outfmt fasta | sed '/^>/! s/[^AGCT]/N/g' > '$fullfname1'");
        unlink("$download_dir/$fname");
        end_fork();
      }

      ## Output sequenceID to taxonomy ID map
      print_header_lines($fullfname1, $taxid, "$assembly_accession $organism_name");
    }
  }

  wait_children();

  print STDERR "  downloaded $downloaded_files files, $existing_files already existed.\n";
}
