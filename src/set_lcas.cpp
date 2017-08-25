// vim: noai:ts=2:sw=2:expandtab:smarttab
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
#include "taxdb.h"
#include "readcounts.hpp"
#include <unordered_map>
#include <map>

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
string DB_filename, Index_filename,
  Output_DB_filename, TaxDB_filename,
  File_to_taxon_map_filename,
  ID_to_taxon_map_filename, Multi_fasta_filename;
bool force_taxid = false;
int New_taxid_start = 1000000000;

bool Allow_extra_kmers = false;
bool verbose = false;
bool Operate_in_RAM = false;
bool One_FASTA_file = false;
bool Add_taxIds_for_Sequences = false;
bool Use_uids_instead_of_taxids = false;
bool Output_UID_map_to_STDOUT = false;
bool Pretend = false;

string UID_map_filename;
ofstream UID_map_file;

uint32_t current_uid = 0;
uint32_t max_uid = -1;
unordered_map<uint32_t, uint32_t> Parent_map;
//unordered_multimap<uint32_t, uint32_t> Children_map;
//typedef std::_Rb_tree_iterator<std::pair<const std::set<unsigned int>, unsigned int> > map_it;
//typedef std::_Rb_tree_iterator<std::pair<const std::vector<unsigned int>, unsigned int> > map_it;
typedef const vector<uint32_t>* map_it;
vector< map_it > UID_to_taxids_vec;
map< vector<uint32_t>, uint32_t> Taxids_to_UID_map;

unordered_map<string, uint32_t> ID_to_taxon_map;
unordered_map<uint32_t, bool> SeqId_added;
KrakenDB Database;
TaxonomyDB<uint32_t, ReadCounts> taxdb;

int main(int argc, char **argv) {
  #ifdef _OPENMP
  omp_set_num_threads(1);
  #endif

  parse_command_line(argc, argv);

  if (!TaxDB_filename.empty() && !force_taxid) {
    taxdb = TaxonomyDB<uint32_t, ReadCounts>(TaxDB_filename);
    for (const auto & tax : taxdb.taxIDsAndEntries) {
      if (tax.first != 0)
        Parent_map[tax.first] = tax.second.parentTaxonomyID;
//      Children_map[tax.second.parentTaxonomyID].insert(tax.first);
    }
    Parent_map[1] = 0;
  } else {
    cerr << "TaxDB argument is required!" << endl;
    return 1;
  }

  if (Use_uids_instead_of_taxids) {
    UID_map_file.open(UID_map_filename, ios_base::out | ios_base::binary);

    if (!UID_map_file.is_open()) {
      cerr << "Something went wrong while creating the file." << endl;
      exit(1);
    }

  }

  QuickFile db_file(DB_filename, "rw");

  char *temp_ptr = NULL;
  size_t db_file_size = db_file.size();
  if (Operate_in_RAM) {
    cerr << "Getting " << DB_filename << " into memory ... ";
    db_file.close_file();
    temp_ptr = new char[ db_file_size ];
    ifstream ifs(DB_filename.c_str(), ifstream::binary);
    ifs.read(temp_ptr, db_file_size);
    ifs.close();
    Database = KrakenDB(temp_ptr);
    cerr << "done" << endl;
  } else {
    if (Output_DB_filename.size() > 0) {
      cerr << "You need to operate in RAM (flag -M) to use output to a different file (flag -o)" << endl;
      return 1;
    }
    //std::ifstream ifs("input.txt", std::ios::binary);
    //std::ofstream ofs("output.txt", std::ios::binary);
    //ofs << ifs.rdbuf();

    Database = KrakenDB(db_file.ptr());
  }

  KmerScanner::set_k(Database.get_k());

  QuickFile idx_file(Index_filename);
  KrakenDBIndex db_index(idx_file.ptr());
  Database.set_index(&db_index);

  if (One_FASTA_file)
    process_single_file();
  else
    process_files();

  if (Operate_in_RAM && !Pretend) {
    if (Output_DB_filename.size() > 0) {
      DB_filename = Output_DB_filename;
    }
    cerr << "Writing database from RAM back to " << DB_filename << " ..." << endl;
    ofstream ofs(DB_filename.c_str(), ofstream::binary);
    ofs.write(temp_ptr, db_file_size);
    ofs.close();
    delete temp_ptr;
  }

  UID_map_file.close();

  // Write new TaxDB file if new taxids were added
  if (Add_taxIds_for_Sequences && !TaxDB_filename.empty() && !Pretend) {
    cerr << "Writing new TaxDB ..." << endl;
    ofstream ofs(TaxDB_filename.c_str());
    taxdb.writeTaxonomyIndex(ofs);
    ofs.close();
  }

  return 0;
}

void process_single_file() {
  cerr << "Processing FASTA files" << endl;
  ifstream map_file(ID_to_taxon_map_filename.c_str());
  if (map_file.rdstate() & ifstream::failbit) {
    err(EX_NOINPUT, "can't open %s", ID_to_taxon_map_filename.c_str());
  }
  string line, seq_id;
  uint32_t parent_taxid, taxid;
  while (map_file.good()) {
    getline(map_file, line);
    if (line.empty())
      break;
    istringstream iss(line);
    iss >> seq_id;
    if (ID_to_taxon_map.find(seq_id) != ID_to_taxon_map.end()) 
        continue;

    if (Add_taxIds_for_Sequences) {
      iss >> parent_taxid;
      taxid = ++New_taxid_start;
      Parent_map[taxid] = parent_taxid;
      auto itEntry = taxdb.taxIDsAndEntries.insert({taxid, TaxonomyEntry<uint32_t, ReadCounts>(taxid, parent_taxid, "sequence")});
      if (!itEntry.second)
          cerr << "Taxonomy ID " << taxid << " already in Taxonomy DB? Shouldn't happen - run set_lcas without the -a option." << endl;
    } else {
      iss >> taxid;
    }
    ID_to_taxon_map[seq_id] = taxid;
  }

  FastaReader reader(Multi_fasta_filename);
  DNASequence dna;
  uint32_t seqs_processed = 0;
  uint32_t seqs_skipped = 0;
  uint32_t seqs_no_taxid = 0;

  while (reader.is_valid()) {
    dna = reader.next_sequence();
    if (! reader.is_valid())
      break;

    if ( dna.seq.empty() ) {
      ++seqs_skipped;
      continue;
    }

    // Get the taxid. If the header specifies kraken:taxid, use that
    uint32_t taxid;
    string prefix = "kraken:taxid|";
    if (dna.id.substr(0,prefix.size()) == prefix) {
        taxid = std::stol(dna.id.substr(prefix.size()));
        if (taxid == 0) {
          cerr << "Error: taxid is zero for the line '" << dna.id << "'?!" << endl;
        }
        const auto strBegin = dna.header_line.find_first_not_of("\t ");
        if (strBegin != std::string::npos)
            dna.header_line = dna.header_line.substr(strBegin);
    } else {
        taxid = ID_to_taxon_map[dna.id];
    }
    
    if (Add_taxIds_for_Sequences) {
      auto entryIt = taxdb.taxIDsAndEntries.find(taxid);
      if (entryIt == taxdb.taxIDsAndEntries.end()) {
        cerr << "Error! Didn't find " << taxid << " in TaxonomyDB!!" << endl;
      } else {
        entryIt->second.scientificName = dna.header_line;
      }
    }

    if (taxid) {
      #pragma omp parallel for schedule(dynamic)
      for (size_t i = 0; i < dna.seq.size(); i += SKIP_LEN)
        set_lcas(taxid, dna.seq, i, i + SKIP_LEN + Database.get_k() - 1);

        ++seqs_processed;
    } else {
      if (verbose) 
        cerr << "Skipping sequence with header [" << dna.header_line << "] - no taxid" << endl;
      ++seqs_no_taxid;
    }

    cerr << "\rProcessed " << seqs_processed << " sequences";
  }
  cerr << "\r                                                                            ";
  cerr << "\rFinished processing " << seqs_processed << " sequences (skipping "<< seqs_skipped <<" empty sequences, and " << seqs_no_taxid<<" sequences with no taxonomy mapping)" << endl;
}

void process_files() {
  cerr << "Processing files in " << File_to_taxon_map_filename.c_str() << endl;
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
      if (! Allow_extra_kmers) {
        errx(EX_DATAERR, "kmer found in sequence that is not in database");
      } 
      else if (verbose) {
        cerr << "kmer found in sequence w/ taxid " << taxid << " that is not in database" << endl;
      }
      continue;
    }
    if (Use_uids_instead_of_taxids) {
      uint32_t kmer_uid = *val_ptr;
      bool new_taxid = kmer_uid == 0;
      vector<uint32_t> taxid_set;
      if (new_taxid) {
        taxid_set.push_back(taxid);
      } else {
        if (kmer_uid > UID_to_taxids_vec.size()) {
          // This can happen when set_lcas is called on a database that is not all zeros
          cerr << "kmer_uid ("<< kmer_uid <<") greater than UID vector size ("<< UID_to_taxids_vec.size()<<")!!" << endl;
          exit(1);
        }
        taxid_set = *(UID_to_taxids_vec.at(kmer_uid-1));
        auto it = std::lower_bound( taxid_set.begin(), taxid_set.end(), taxid); // find proper position in descending order

        if (it == taxid_set.end() || *it != taxid) {
          // add the taxid to the set, in the right position
           taxid_set.insert( it, taxid ); // insert before iterator it
           new_taxid = true;
        }
      }

      if (new_taxid) {
        if (max_uid <= current_uid) {
          cerr << "Maxxed out on the UIDs!!" << endl;
          exit(1);
        }

        // get a new taxid for this set
        #pragma omp critical(new_uid)
        {
        auto insert_res = Taxids_to_UID_map.insert( { std::move(taxid_set), current_uid + 1 } );
        if (insert_res.second) {
          ++current_uid;

          // print result for map:
          if (Output_UID_map_to_STDOUT) {
            auto tid_it = insert_res.first->first.begin();
            cout << current_uid << '\t' << *tid_it++; 
            while (tid_it != insert_res.first->first.end()) { cout << ' ' << *tid_it++; }
            cout << '\n';
          }

          // FORMAT: TAXID<uint32_t> PARENT<uint32_t>
          // TODO: Consider using mmap here
          UID_map_file.write((char*)&taxid, sizeof(taxid));
          UID_map_file.write((char*)&kmer_uid, sizeof(kmer_uid));

          //UID_to_taxids_vec[current_uid] = taxid_set;
          UID_to_taxids_vec.push_back( &(insert_res.first->first) );
          *val_ptr = current_uid;
        } else {
         *val_ptr = insert_res.first->second;
        }
        }
      }
    } else if (!force_taxid) {
      *val_ptr = lca(Parent_map, taxid, *val_ptr);
    } else {
      // When force_taxid is set, do not compute lca, but assign the taxid
      // of the (last) sequence to k-mers
      *val_ptr = taxid;
    }
  }
}

void parse_command_line(int argc, char **argv) {
  int opt;
  long long sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "f:d:i:t:n:m:F:xMTvb:apI:o:S")) != -1) {
    switch (opt) {
      case 'f' :
        File_to_taxon_map_filename = optarg;
        break;
      case 'I' :
        Use_uids_instead_of_taxids = true;
        UID_map_filename = optarg;
        break;
      case 'S' :
        Output_UID_map_to_STDOUT = true;
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
      case 'T' :
        force_taxid = true;
        break;
      case 'v' :
        verbose = true;
        break;
      case 'x' :
        Allow_extra_kmers = true;
        break;
      case 'a' :
        Add_taxIds_for_Sequences = true;
        break;
      case 'b' :
        TaxDB_filename = optarg;
        break;
      case 'M' :
        Operate_in_RAM = true;
        break;
      case 'o' :
        Output_DB_filename = optarg;
        break;
      case 'p' :
        Pretend = true;
        break;
      default:
        usage();
        break;
    }
  }

  if (DB_filename.empty() || Index_filename.empty() ||
      TaxDB_filename.empty())
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
       << "* -b filename      Taxonomy DB file" << endl
       << "  -t #             Number of threads" << endl
       << "  -M               Copy DB to RAM during operation" << endl
       << "  -o filename      Output database to filename, instead of overwriting the input database" << endl
       << "  -x               K-mers not found in DB do not cause errors" << endl
       << "  -f filename      File to taxon map" << endl
       << "  -F filename      Multi-FASTA file with sequence data" << endl
       << "  -m filename      Sequence ID to taxon map" << endl
       << "  -a               Add taxonomy IDs (starting with "<<New_taxid_start<<") for sequences to Taxonomy DB" << endl
       << "  -T               Do not set LCA as taxid for kmers, but the taxid of the sequence" << endl
       << "  -I filename      Write UIDs into database, and output (binary) UID-to-taxid map to filename" << endl
       << "  -S               Write UID-to-taxid map to STDOUT" << endl
       << "  -p               Pretend - do not write database back to disk (when working in RAM)" << endl
       << "  -v               Verbose output" << endl
       << "  -h               Print this message" << endl
       << endl
       << "-F and -m must be specified together.  If -f is given, "
       << "-F/-m are ignored." << endl;
  exit(exit_code);
}

