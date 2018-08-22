/*
 * Original file Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
 * Portions (c) 2017-2018, Florian Breitwieser <fbreitwieser@jhu.edu> as part of KrakenUniq
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
#include "taxdb.hpp"
#include "uid_mapping.hpp"
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
void set_lcas(uint32_t taxid, string &seq, size_t start, size_t finish, bool is_contaminant_taxid = false);

int Num_threads = 1;
string DB_filename, Index_filename,
  Output_DB_filename, TaxDB_filename,
  Kmer_count_filename,
  File_to_taxon_map_filename,
  ID_to_taxon_map_filename, Multi_fasta_filename;
bool Force_contaminant_taxid = false;
uint32_t New_taxid_start = 1000000000;

bool Allow_extra_kmers = false;
bool verbose = false;
bool Operate_in_RAM = false;
bool One_FASTA_file = false;
bool Add_taxIds_for_Assembly = false;
bool Add_taxIds_for_Sequences = false;
bool Use_uids_instead_of_taxids = false;
bool Output_UID_map_to_STDOUT = false;
bool Pretend = false;
int Minimum_sequence_size = -1;

string UID_map_filename;
ofstream UID_map_file;

uint32_t current_uid = 0;
unordered_map<uint32_t, uint32_t> Parent_map;
//unordered_multimap<uint32_t, uint32_t> Children_map;
//typedef std::_Rb_tree_iterator<std::pair<const std::set<unsigned int>, unsigned int> > map_it;
//typedef std::_Rb_tree_iterator<std::pair<const std::vector<unsigned int>, unsigned int> > map_it;
vector< const TaxidSet*  > UID_to_taxids_vec;
map< TaxidSet, uint32_t> Taxids_to_UID_map;

unordered_map<string, uint32_t> ID_to_taxon_map;
unordered_map<uint32_t, bool> SeqId_added;
KrakenDB Database;
TaxonomyDB<uint32_t> taxdb;

const string prefix = "kraken:taxid|";

// do not add sequence taxIDs for host sequences (currently only human and mouse)
const uint32_t TID_HUMAN = 9606;
const uint32_t TID_MOUSE = 10090;

// k-mers appearing in contaminant sequences will keep the contaminant
//  sequence taxid, even if they also appear in a genome
const uint32_t TID_CONTAMINANT1 = 32630; // 'synthetic construct'
const uint32_t TID_CONTAMINANT2 = 81077; // 'artificial sequences'


int main(int argc, char **argv) {
  #ifdef _OPENMP
  omp_set_num_threads(1);
  #endif

  parse_command_line(argc, argv);

  if (!TaxDB_filename.empty()) {
    taxdb = TaxonomyDB<uint32_t>(TaxDB_filename);
    Parent_map = taxdb.getParentMap();
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

  if (!Operate_in_RAM && Output_DB_filename.size() > 0) {
      cerr << "You need to operate in RAM (flag -M) to use output to a different file (flag -o)" << endl;
      return 1;
  }

  QuickFile db_file(DB_filename, "rw");
  size_t db_file_size = db_file.size();
  vector<char> dat;
  if (Operate_in_RAM) {
    db_file.close_file();
    dat = slurp_file(DB_filename, db_file_size);
    Database = KrakenDB(dat.data());
  } else {
    if (Output_DB_filename.size() > 0) {
      //system("cp " + DB_filename + " " + Output_DB_filename);
    }
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

  if (!Kmer_count_filename.empty()) {
    ofstream ofs(Kmer_count_filename.c_str());
    cerr << "Writing kmer counts to " << Kmer_count_filename << "..." << endl;
    auto counts = Database.count_taxons();
    for (auto it = counts.begin(); it != counts.end(); ++it) {
      ofs << it->first << '\t' << it->second << '\n';
    }
    ofs.close();
  }

  if (Operate_in_RAM && !Pretend) {
    if (Output_DB_filename.size() > 0) {
      DB_filename = Output_DB_filename;
    }
    cerr << "Writing database from RAM back to " << DB_filename << " ..." << endl;
    ofstream ofs(DB_filename.c_str(), ofstream::binary);
    ofs.write(dat.data(), db_file_size);
    ofs.close();
    dat.clear();
  }

  UID_map_file.close();

  // Write new TaxDB file if new taxids were added
  if ((Add_taxIds_for_Sequences || Add_taxIds_for_Assembly) && !TaxDB_filename.empty() && !Pretend) {
    cerr << "Writing new TaxDB ..." << endl;
    ofstream ofs(TaxDB_filename.c_str());
    taxdb.writeTaxonomyIndex(ofs);
    ofs.close();
  }

  return 0;
}

inline 
uint32_t get_new_taxid(
    unordered_map<string, uint32_t>& name_to_taxid_map, 
    unordered_map<uint32_t,uint32_t>& Parent_map,
    string name, uint32_t parent_taxid, const string & rank_name) {

  auto it = name_to_taxid_map.find(name);
  if (it == name_to_taxid_map.end()) {
    uint32_t new_taxid = ++New_taxid_start;
    bool insert_res = taxdb.insert(new_taxid, parent_taxid, rank_name, name);
    //cerr << "Adding assembly: " << name << " with taxid " << new_taxid;
    if (!insert_res) {
      return 0;
    }
    // insert_res shows if insert failed, but we don't care
    Parent_map[new_taxid] = parent_taxid;
    name_to_taxid_map[name] = new_taxid;
    return new_taxid;
   } else {
    return it->second;
   }
}

unordered_map<string,uint32_t> read_seqid_to_taxid_map(string ID_to_taxon_map_filename, 
    TaxonomyDB<uint32_t>& taxdb, unordered_map<uint32_t,uint32_t>& Parent_map, 
    bool Add_taxIds_for_Assembly, bool Add_taxIds_for_Sequences) {

  cerr << "Reading sequence ID to taxonomy ID mapping ... ";

  unordered_map<string, uint32_t> ID_to_taxon_map;
  ifstream map_file(ID_to_taxon_map_filename.c_str());
  if (map_file.rdstate() & ifstream::failbit) {
    err(EX_NOINPUT, "can't open %s", ID_to_taxon_map_filename.c_str());
  }
  string line, seq_id, name;
  uint32_t taxid;

  if (Add_taxIds_for_Assembly || Add_taxIds_for_Sequences) {
    for (auto it = taxdb.entries.begin(); it != taxdb.entries.end(); ++it) {
      if (it->first >= New_taxid_start) {
        New_taxid_start = it->first+100;
      } 
    }
    cerr << "[starting new taxonomy IDs with " << (New_taxid_start+1) << ']';
  }

  // Used when adding new taxids for assembly or sequence
  unordered_map<string, uint32_t> name_to_taxid_map;

  while (map_file.good()) {
    getline(map_file, line);
    if (line.empty())
      break;
    istringstream iss(line);
    iss >> seq_id >> taxid;

    auto it = ID_to_taxon_map.find(seq_id);
    if (it != ID_to_taxon_map.end()) {
      // The sequence ID has been seen before, ignore
      continue;
    }

    uint32_t orig_taxid = taxid;

    if (Add_taxIds_for_Assembly && iss.good()) {
      iss.get();
      getline(iss, name);
      if (!name.empty())
        taxid = get_new_taxid(name_to_taxid_map, Parent_map, name, taxid, "assembly");
    }

    if (Add_taxIds_for_Sequences && orig_taxid != TID_HUMAN && orig_taxid != TID_MOUSE) {
      taxid = get_new_taxid(name_to_taxid_map, Parent_map, seq_id, taxid, "sequence");
    }
    if (Add_taxIds_for_Assembly || Add_taxIds_for_Sequences) {
      cout << seq_id << '\t' << taxid << '\n';
    }
    ID_to_taxon_map[seq_id] = taxid;
  }
  if (ID_to_taxon_map.size() == 0) {
    cerr << "Error: No ID mappings present!!" << endl;
  }
  cerr << " got " << ID_to_taxon_map.size() << " mappings." << endl;
  return ID_to_taxon_map;
}

void process_single_file() {
  cerr << "Processing FASTA files" << endl;
 
  ID_to_taxon_map = read_seqid_to_taxid_map(ID_to_taxon_map_filename, taxdb, Parent_map, Add_taxIds_for_Assembly, Add_taxIds_for_Sequences);

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
    auto it = ID_to_taxon_map.find(dna.id);
    if (it != ID_to_taxon_map.end()) {
      taxid = it->second;
    } else if (dna.id.size() >= prefix.size() && dna.id.substr(0,prefix.size()) == prefix) {
      // if the AC is not in the map, check if the fasta entry starts with '>kraken:taxid'
        taxid = std::stol(dna.id.substr(prefix.size()));
        if (taxid == 0) {
          cerr << "Error: taxonomy ID is zero for sequence '" << dna.id << "'?!" << endl;
        }
        const auto strBegin = dna.header_line.find_first_not_of("\t ");
        if (strBegin != std::string::npos)
            dna.header_line = dna.header_line.substr(strBegin);
    } else {
        cerr << "Error! Didn't find taxonomy ID mapping for sequence " <<  dna.id << "!!" << endl;
        ++seqs_skipped;
        continue;
    }

	if (Minimum_sequence_size > 0 && dna.seq.size() < Minimum_sequence_size) {
      cerr << "Skipping sequence " << dna.id << " as it's too short (" << dna.seq.size() << ")" << endl;
	  ++ seqs_skipped;
	  continue;
	}

    auto it_p = Parent_map.find(taxid);
    if (it_p == Parent_map.end()) {
      cerr << "Skipping sequence " << dna.id << " since taxonomy ID " << taxid << " is not in taxonomy database!" << endl;
      ++ seqs_skipped;
      continue;
    }
    
    bool is_contaminant_taxid = taxid == TID_CONTAMINANT1 || taxid == TID_CONTAMINANT2;
    if (Add_taxIds_for_Sequences && taxid != TID_HUMAN && it_p->second != TID_HUMAN && taxid != TID_MOUSE && it_p->second != TID_MOUSE) {
      // Update entry based on header line
      auto entryIt = taxdb.entries.find(taxid);
      if (entryIt == taxdb.entries.end()) {
        cerr << "Error! Didn't find taxid " << taxid << " in TaxonomyDB - can't update it!! ["<<dna.header_line<<"]" << endl;
      } else {
        entryIt->second.scientificName = dna.header_line;
      }
    }

    // TODO: Allow exclusion of certain taxids in the building process
    //if (Excluded_taxons.count(taxid) > 0) {
      // exclude taxid!
    //}

    if (taxid) {
      if (Parent_map.find(taxid) == Parent_map.end() || taxdb.entries.find(taxid) == taxdb.entries.end()) {
        cerr << "Ignoring sequence for taxID " << taxid << " - not in taxDB\n";
      } else {
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < dna.seq.size(); i += SKIP_LEN)
          set_lcas(taxid, dna.seq, i, i + SKIP_LEN + Database.get_k() - 1, is_contaminant_taxid);
         ++seqs_processed;
      }
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
    // TODO: Support a mapping file with only file names, not taxids
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

//void process_sequence(DNASequence dna) {
  // TODO: Refactor such that a list of files + taxid can be given.
  // Or maybe asembly_summary file?
//}

void set_lcas(uint32_t taxid, string &seq, size_t start, size_t finish, bool is_contaminant_taxid) {
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

    // TODO: Should I use pragma omp critical here?
    if (Use_uids_instead_of_taxids) {
      #pragma omp critical(new_uid)
      *val_ptr = uid_mapping(Taxids_to_UID_map, UID_to_taxids_vec, taxid, *val_ptr, current_uid, UID_map_file);
    } else {
      if (!Force_contaminant_taxid) {
        *val_ptr = lca(Parent_map, taxid, *val_ptr);
      } else {
        if (*val_ptr == TID_CONTAMINANT1 || *val_ptr == TID_CONTAMINANT2) {
          // keep value
        } else if (is_contaminant_taxid) {
          // When Force_contaminant_taxid is set, do not compute lca, but assign the taxid
          // of the (last) sequence to k-mers
          *val_ptr = taxid;
        } else {
          *val_ptr = lca(Parent_map, taxid, *val_ptr);
        }
      }
    }
  }
}

void parse_command_line(int argc, char **argv) {
  int opt;
  long long sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "f:d:i:t:n:m:F:xMTvb:aApI:o:Sc:E:")) != -1) {
    switch (opt) {
      case 'f' :
        File_to_taxon_map_filename = optarg;
        break;
      case 'I' :
        Use_uids_instead_of_taxids = true;
        UID_map_filename = optarg;
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
        Force_contaminant_taxid = true;
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
      case 'A' :
        Add_taxIds_for_Assembly = true;
        break;
      case 'b' :
        TaxDB_filename = optarg;
        break;
      case 'c' :
        Kmer_count_filename = optarg;
        break;
      case 'M' :
        Operate_in_RAM = true;
        break;
      case 'o' :
        Output_DB_filename = optarg;
        break;
      case 'E' :
        Minimum_sequence_size = atoi(optarg);
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
       << "  -a               Add taxonomy IDs (starting with "<<(New_taxid_start+1)<<") for assemblies (third column in seqid2taxid.map) to Taxonomy DB" << endl
       << "  -A               Add taxonomy IDs (starting with "<<(New_taxid_start+1)<<") for sequences to Taxonomy DB" << endl
       //<< "  -T               Do not set LCA as taxid for kmers, but the taxid of the sequence" << endl
       << "  -T               When a k-mer appears in a 'synthetic construct' sequence, force the taxID to be the 'synthetic construct' taxID, instead of the LCA." << endl
	   << "  -E #             Exclude sequences that are shorter than the threshold." << endl
       << "  -I filename      Write UIDs into database, and output (binary) UID-to-taxid map to filename" << endl
       << "  -p               Pretend - do not write database back to disk (when working in RAM)" << endl
       << "  -v               Verbose output" << endl
       << "  -h               Print this message" << endl
       << endl
       << "-F and -m must be specified together.  If -f is given, "
       << "-F/-m are ignored." << endl;
  exit(exit_code);
}

