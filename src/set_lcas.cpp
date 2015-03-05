/*
 * Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken taxonomic sequence classification system.
 *
 * Kraken is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Kraken is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kraken.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "kraken_headers.hpp"
#include "quickfile.hpp"
#include "krakendb.hpp"
#include "krakenutil.hpp"
#include "seqreader.hpp"

#define SKIP_LEN 50000

using namespace std;
using namespace kraken;

void parse_command_line(int argc, char **argv);
void usage(int exit_code=EX_USAGE);
void process_files();
void process_single_file();
void process_file(string filename, uint32_t taxid);
void set_lcas(uint32_t taxid, string &seq, size_t start, size_t finish);

int Num_threads = 1;
string DB_filename, Index_filename, Nodes_filename,
  File_to_taxon_map_filename,
  ID_to_taxon_map_filename, Multi_fasta_filename;
bool Allow_extra_kmers = false;
bool Operate_in_RAM = false;
bool One_FASTA_file = false;
map<uint32_t, uint32_t> Parent_map;
map<string, uint32_t> ID_to_taxon_map;
KrakenDB Database;

int main(int argc, char **argv) {
  #ifdef _OPENMP
  omp_set_num_threads(1);
  #endif

  parse_command_line(argc, argv);
  Parent_map = build_parent_map(Nodes_filename);

  QuickFile db_file(DB_filename, "rw");
  Database = KrakenDB(db_file.ptr());
  KmerScanner::set_k(Database.get_k());

  char *temp_ptr = NULL;
  size_t db_file_size = db_file.size();
  if (Operate_in_RAM) {
    db_file.close_file();
    temp_ptr = new char[ db_file_size ];
    ifstream ifs(DB_filename.c_str(), ifstream::binary);
    ifs.read(temp_ptr, db_file_size);
    ifs.close();
    Database = KrakenDB(temp_ptr);
  }

  QuickFile idx_file(Index_filename);
  KrakenDBIndex db_index(idx_file.ptr());
  Database.set_index(&db_index);

  if (One_FASTA_file)
    process_single_file();
  else
    process_files();

  if (Operate_in_RAM) {
    ofstream ofs(DB_filename.c_str(), ofstream::binary);
    ofs.write(temp_ptr, db_file_size);
    ofs.close();
    delete temp_ptr;
  }

  return 0;
}

void process_single_file() {
  ifstream map_file(ID_to_taxon_map_filename.c_str());
  if (map_file.rdstate() & ifstream::failbit) {
    err(EX_NOINPUT, "can't open %s", ID_to_taxon_map_filename.c_str());
  }
  string line;
  while (map_file.good()) {
    getline(map_file, line);
    if (line.empty())
      break;
    string seq_id;
    uint32_t taxid;
    istringstream iss(line);
    iss >> seq_id;
    iss >> taxid;
    ID_to_taxon_map[seq_id] = taxid;
  }

  FastaReader reader(Multi_fasta_filename);
  DNASequence dna;
  uint32_t seqs_processed = 0;

  while (reader.is_valid()) {
    dna = reader.next_sequence();
    if (! reader.is_valid())
      break;
    uint32_t taxid = ID_to_taxon_map[dna.id];
    if (taxid) {
      #pragma omp parallel for schedule(dynamic)
      for (size_t i = 0; i < dna.seq.size(); i += SKIP_LEN)
        set_lcas(taxid, dna.seq, i, i + SKIP_LEN + Database.get_k() - 1);
    }
    cerr << "\rProcessed " << ++seqs_processed << " sequences";
  }
  cerr << "\r                                                       ";
  cerr << "\rFinished processing " << seqs_processed << " sequences" << endl;
}

void process_files() {
  ifstream map_file(File_to_taxon_map_filename.c_str());
  if (map_file.rdstate() & ifstream::failbit) {
    err(EX_NOINPUT, "can't open %s", File_to_taxon_map_filename.c_str());
  }
  string line;
  uint32_t seqs_processed = 0;

  while (map_file.good()) {
    getline(map_file, line);
    if (line.empty())
      break;
    string filename;
    uint32_t taxid;
    istringstream iss(line);
    iss >> filename;
    iss >> taxid;
    process_file(filename, taxid);
    cerr << "\rProcessed " << ++seqs_processed << " sequences";
  }
  cerr << "\r                                                       ";
  cerr << "\rFinished processing " << seqs_processed << " sequences" << endl;
}

void process_file(string filename, uint32_t taxid) {
  FastaReader reader(filename);
  DNASequence dna;
  
  // For the purposes of this program, we assume these files are
  // single-fasta files.
  dna = reader.next_sequence();

  #pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < dna.seq.size(); i += SKIP_LEN)
    set_lcas(taxid, dna.seq, i, i + SKIP_LEN + Database.get_k() - 1);
}

void set_lcas(uint32_t taxid, string &seq, size_t start, size_t finish) {
  KmerScanner scanner(seq, start, finish);
  uint64_t *kmer_ptr;
  uint32_t *val_ptr;

  while ((kmer_ptr = scanner.next_kmer()) != NULL) {
    if (scanner.ambig_kmer())
      continue;
    val_ptr = Database.kmer_query(
                Database.canonical_representation(*kmer_ptr)
              );
    if (val_ptr == NULL) {
      if (! Allow_extra_kmers)
        errx(EX_DATAERR, "kmer found in sequence that is not in database");
      else
        continue;
    }
    *val_ptr = lca(Parent_map, taxid, *val_ptr);
  }
}

void parse_command_line(int argc, char **argv) {
  int opt;
  long long sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "f:d:i:t:n:m:F:xM")) != -1) {
    switch (opt) {
      case 'f' :
        File_to_taxon_map_filename = optarg;
        break;
      case 'd' :
        DB_filename = optarg;
        break;
      case 'i' :
        Index_filename = optarg;
        break;
      case 'F' :
        Multi_fasta_filename = optarg;
        break;
      case 'm' :
        ID_to_taxon_map_filename = optarg;
        break;
      case 't' :
        sig = atoll(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive thread count");
        #ifdef _OPENMP
        if (sig > omp_get_num_procs())
          errx(EX_USAGE, "thread count exceeds number of processors");
        Num_threads = sig;
        omp_set_num_threads(Num_threads);
        #endif
        break;
      case 'n' :
        Nodes_filename = optarg;
        break;
      case 'x' :
        Allow_extra_kmers = true;
        break;
      case 'M' :
        Operate_in_RAM = true;
        break;
      default:
        usage();
        break;
    }
  }

  if (DB_filename.empty() || Index_filename.empty() ||
      Nodes_filename.empty())
    usage();
  if (File_to_taxon_map_filename.empty() &&
      (Multi_fasta_filename.empty() || ID_to_taxon_map_filename.empty()))
    usage();

  if (! File_to_taxon_map_filename.empty())
    One_FASTA_file = false;
  else
    One_FASTA_file = true;
}

void usage(int exit_code) {
  cerr << "Usage: set_lcas [options]" << endl
       << endl
       << "Options: (*mandatory)" << endl
       << "* -d filename      Kraken DB filename" << endl
       << "* -i filename      Kraken DB index filename" << endl
       << "* -n filename      NCBI Taxonomy nodes file" << endl
       << "  -t #             Number of threads" << endl
       << "  -M               Copy DB to RAM during operation" << endl
       << "  -x               K-mers not found in DB do not cause errors" << endl
       << "  -f filename      File to taxon map" << endl
       << "  -F filename      Multi-FASTA file with sequence data" << endl
       << "  -m filename      Sequence ID to taxon map" << endl
       << "  -h               Print this message" << endl
       << endl
       << "-F and -m must be specified together.  If -f is given, "
       << "-F/-m are ignored." << endl;
  exit(exit_code);
}
