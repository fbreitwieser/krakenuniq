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
#include "krakendb.hpp"
#include "krakenutil.hpp"
#include "quickfile.hpp"
#include "readcounts.hpp"
#include "seqreader.hpp"
#include "taxdb.hpp"
#include "uid_mapping.hpp"
#include <sstream>

const size_t DEF_WORK_UNIT_SIZE = 500000;
int New_taxid_start = 1000000000;

using namespace std;
using namespace kraken;

#define USE_KHSET_FOR_EXACT_COUNTING

#ifdef EXACT_COUNTING
  #ifdef USE_KHSET_FOR_EXACT_COUNTING
    #include "khset.h"
    using READCOUNTS = ReadCounts< kh::khset64_t >;
  #else
    #include <unordered_set>
    using READCOUNTS = ReadCounts< unordered_set<uint64_t> >;
  #endif
#else
  using READCOUNTS = ReadCounts<HyperLogLogPlusMinus<uint64_t> >;
#endif

void parse_command_line(int argc, char **argv);
void usage(int exit_code = EX_USAGE);
void process_file(char *filename, managed_ostream &kraken_output,
                  managed_ostream &classified_output,
                  managed_ostream &unclassified_output);
bool classify_sequence(DNASequence &dna, ostringstream &koss,
                       ostringstream &coss, ostringstream &uoss,
                       unordered_map<uint32_t, READCOUNTS>&);
inline void print_sequence(ostringstream* oss_ptr, const DNASequence& dna);
string hitlist_string(const vector<uint32_t> &taxa, const vector<char>& ambig_list);


set<uint32_t> get_ancestry(uint32_t taxon);
void report_stats(struct timeval time1, struct timeval time2);
double get_seconds(struct timeval time1, struct timeval time2);
unordered_map<uint32_t, READCOUNTS> taxon_counts; // stats per taxon

int Num_threads = 1;
vector<string> DB_filenames;
vector<string> Index_filenames;
vector<string> Onefile_DB_filenames;
bool Lock_DB = false;
bool Quick_mode = false;
bool Ordered_mode = false;
bool Fastq_input = false;
bool Print_classified = false;
bool Print_unclassified = false;
bool Print_kraken = true;
bool Print_kraken_report = false;
bool Populate_memory = false;
bool Only_classified_kraken_output = false;
bool Print_sequence = false;
bool Print_Progress = true;
bool full_report = false;

bool Map_UIDs = false;
string UID_to_TaxID_map_filename;
map<uint32_t, vector<uint32_t> > UID_to_taxids_map;
QuickFile UID_to_TaxID_map_file;

uint32_t Minimum_hit_count = 1;
unordered_map<uint32_t, uint32_t> Parent_map;
unordered_map<uint32_t, vector<uint32_t>> Uid_dict;
string TaxDB_file;

vector<string> Kraken_output_files;
vector<string> Report_output_files;
vector<string> Classified_output_files;
vector<string> Unclassified_output_files;
size_t Work_unit_size = DEF_WORK_UNIT_SIZE;
TaxonomyDB<uint32_t> taxdb;
static vector<KrakenDB *> KrakenDatabases(DB_filenames.size());

struct db_status {
  db_status() : current_bin_key(0), current_min_pos(1), current_max_pos(0) {}
  uint64_t current_bin_key;
  int64_t current_min_pos;
  int64_t current_max_pos;
};

unsigned long long total_classified = 0;
unsigned long long total_sequences = 0;
unsigned long long total_bases = 0;
uint32_t ambig_taxon = -1;

void load_kraken_db(KrakenDB &database, string DB_filename,
                  string Index_filename) {
  QuickFile db_file;
  db_file.open_file(DB_filename);
  if (Populate_memory) {
    db_file.load_file();
  }
  database = KrakenDB(db_file.ptr());
  QuickFile idx_file;
  idx_file.open_file(Index_filename);
  if (Populate_memory)
    idx_file.load_file();

  KrakenDBIndex db_index(idx_file.ptr());
  database.set_index(&db_index);
}

void load_onefile_krakendb(string filename, TaxonomyDB<uint32_t> &taxdb) {
  ifstream in(filename.c_str(), ios::in | ios::binary);
  string counts_size_s;
  string idx_size_s;
  string db_size_s;
  std::getline(in, counts_size_s);
  unsigned long long counts_size = std::stoull(counts_size_s);
  std::getline(in, idx_size_s);
  unsigned long long idx_size = std::stoull(idx_size_s);
  std::getline(in, db_size_s);
  unsigned long long db_size = std::stoull(db_size_s);

  std::string counts_s(counts_size, '\0');
  in.read(&counts_s[0], counts_size);
  std::stringstream ss(counts_s);
  taxdb.readGenomeSizes(ss);

  char *buffer;
  buffer = (char *)malloc(idx_size);
  in.read(buffer, idx_size);
  KrakenDBIndex *db_index = new KrakenDBIndex(buffer);

  buffer = (char *)malloc(db_size);
  in.read(buffer, db_size);
  KrakenDB *database = new KrakenDB(buffer);
  database->set_index(db_index);
  KrakenDatabases.push_back(database);
}

void ensure_counts_file(vector<string>& db_filenames) {
  for (size_t i = 0; i < db_filenames.size(); ++i) {
    const auto fname = db_filenames[i] + ".counts";
    ifstream ifs(fname);
    bool counts_file_gd = false;
    if (ifs.good()) {
      if (ifs.peek() == std::ifstream::traits_type::eof()) {
        cerr << "Kmer counts file is empty - trying to regenerate ..." <<
          endl;
      } else {
        ifs.close();
        counts_file_gd = true;
      }
    }
    if (!counts_file_gd) {
      ofstream ofs(fname);
      cerr << "Writing kmer counts to " << fname << "... [only once for this database, may take a while] " << endl;
      auto counts = KrakenDatabases[i]->count_taxons();
      for (auto it = counts.begin(); it != counts.end(); ++it) {
        ofs << it->first << '\t' << it->second << '\n';
      }
      ofs.close();
    }
    taxdb.readGenomeSizes(fname);
  }
}

int main(int argc, char **argv) {
  #ifdef _OPENMP
  omp_set_num_threads(1);
  #endif

  parse_command_line(argc, argv);

  if (Map_UIDs) {
    if (DB_filenames.size() > 1) {
      cerr << "Cannot use more than one database with UID mapping!" << endl;
      return 1;
    }

    cerr << "Reading UID mapping file " << UID_to_TaxID_map_filename << endl;
    UID_to_TaxID_map_file.open_file(UID_to_TaxID_map_filename);

    // Always Populate memory
    //if (Populate_memory) {
    UID_to_TaxID_map_file.load_file();
    //}
  }

  if (!TaxDB_file.empty()) {
    taxdb = TaxonomyDB<uint32_t>(TaxDB_file, false);
    Parent_map = taxdb.getParentMap();
  } else {
    cerr << "TaxDB argument is required!" << endl;
    return 1;
  }

  if (Populate_memory)
    cerr << "Loading database(s)... " << endl;

  static vector<QuickFile> idx_files (DB_filenames.size());
  static vector<QuickFile> db_files (DB_filenames.size());
  static vector<KrakenDBIndex> db_indices (DB_filenames.size());


  // TODO: Check DB_filenames and Index_filesnames have the same length
  for (size_t i=0; i < DB_filenames.size(); ++i) {
    cerr << " Database " << DB_filenames[i] << endl;
    db_files[i].open_file(DB_filenames[i]);
    if (Populate_memory)
      db_files[i].load_file();

    KrakenDatabases.push_back(new KrakenDB(db_files[i].ptr()));
    idx_files[i].open_file(Index_filenames[i]);
    if (Populate_memory)
      idx_files[i].load_file();
    db_indices[i] = KrakenDBIndex(idx_files[i].ptr());
    KrakenDatabases[i]->set_index(&db_indices[i]);
  }

  for (size_t i = 0; i < Onefile_DB_filenames.size(); ++i) {
    load_onefile_krakendb(Onefile_DB_filenames[i], taxdb);
  }

  // TODO: Check all databases have the same k
  uint8_t kmer_size = KrakenDatabases[0]->get_k();
  for (size_t i = 1; i < KrakenDatabases.size(); ++i) {
    uint8_t kmer_size_i = KrakenDatabases[i]->get_k();
    if (kmer_size_i != kmer_size) {
      fprintf(stderr, "Different k-mer sizes in databases 1 and %lu: %i vs %i!\n", i+1, (int)kmer_size, (int)kmer_size_i);
      exit(1);
    }
  };
  KmerScanner::set_k(kmer_size);

  if (Populate_memory)
    cerr << "\ncomplete." << endl;

  size_t n_inputs = argc - optind;
  for (int i = optind, j = 0; i < argc; i++, j++) {

    total_classified = 0;
    total_sequences = 0;
    total_bases = 0;

    managed_ostream classified_output(Classified_output_files[j],
                                      Print_classified);
    managed_ostream unclassified_output(Unclassified_output_files[j],
                                        Print_unclassified);
    cerr << "Writing kraken output to: " << Kraken_output_files[j] << endl;
    managed_ostream kraken_output(Kraken_output_files[j], Print_kraken);

    struct timeval tv1, tv2;
    gettimeofday(&tv1, NULL);

    process_file(argv[i], kraken_output, classified_output,
                 unclassified_output);

    gettimeofday(&tv2, NULL);
    report_stats(tv1, tv2);

    if (Report_output_files.size() == n_inputs &&
        !Report_output_files[j].empty() && Report_output_files[j] != "off") {
      gettimeofday(&tv1, NULL);
      std::cerr << "Writing report file to " << Report_output_files[j]
                << "  ..\n";
      managed_ostream report_output(Report_output_files[j], true, true);

      TaxReport<uint32_t, READCOUNTS> rep = TaxReport<uint32_t, READCOUNTS>(
                                                                            *report_output, taxdb, taxon_counts, false);
      if (HLL_PRECISION > 0) {
        if (full_report) {
          rep.setReportCols(vector<string>{
              "%", "reads", "taxReads", "kmers", "taxKmers", "kmersDB",
                "taxKmersDB", "dup", "cov", "taxID", "rank", "taxName"});
        } else {
          rep.setReportCols(vector<string>{"%", "reads", "taxReads", "kmers",
                "dup", "cov", "taxID", "rank",
                "taxName"});
        }
      } else {
        rep.setReportCols(vector<string>{"%", "reads", "taxReads", "taxID",
              "rank", "taxName"});
      }
      rep.printReport("kraken");
      gettimeofday(&tv2, NULL);
      fprintf(stderr, "Report finished in %.3f seconds.\n",
              get_seconds(tv1, tv2));
    }
  }
  cerr << "Finishing up ..." << endl;
  return 0;
}

double get_seconds(struct timeval time1, struct timeval time2) {
  time2.tv_usec -= time1.tv_usec;
  time2.tv_sec -= time1.tv_sec;
  if (time2.tv_usec < 0) {
    time2.tv_sec--;
    time2.tv_usec += 1000000;
  }
  double seconds = time2.tv_usec;
  seconds /= 1e6;
  seconds += time2.tv_sec;
  return(seconds);
}

void report_stats(struct timeval time1, struct timeval time2) {
  double seconds = get_seconds(time1, time2);

  cerr << "\r";
  fprintf(stderr, 
          "%llu sequences (%.2f Mbp) processed in %.3fs (%.1f Kseq/m, %.2f Mbp/m).\n",
          (unsigned long long) total_sequences, total_bases / 1.0e6, seconds,
          total_sequences / 1.0e3 / (seconds / 60),
          total_bases / 1.0e6 / (seconds / 60) );
  fprintf(stderr, "  %llu sequences classified (%.2f%%)\n",
          (unsigned long long) total_classified, total_classified * 100.0 / total_sequences);
  fprintf(stderr, "  %llu sequences unclassified (%.2f%%)\n",
          (unsigned long long) (total_sequences - total_classified),
          (total_sequences - total_classified) * 100.0 / total_sequences);
}

void process_file(char *filename, managed_ostream &kraken_output,
                  managed_ostream &classified_output,
                  managed_ostream &unclassified_output) {
  string file_str(filename);
  DNASequenceReader *reader;
  DNASequence dna;

  if (Fastq_input)
    reader = new FastqReader(file_str);
  else
    reader = new FastaReader(file_str);
  if (Ordered_mode) {
    vector<DNASequence> work_unit;

    while (reader->is_valid()) {
      work_unit.clear();
      size_t total_nt = 0;

      while (total_nt < Work_unit_size) {
        dna = reader->next_sequence();
        if (!reader->is_valid())
          break;
        work_unit.push_back(dna);
        total_nt += dna.seq.size();
      }
      if (total_nt == 0)
        break;
      #pragma omp parallel
      {
        #ifdef _OPENMP
        #pragma omp for ordered schedule(static)
        #endif
        for (size_t j = 0; j < work_unit.size(); j++) {
          ostringstream kraken_output_ss, classified_output_ss,
            unclassified_output_ss;
          unordered_map<uint32_t, READCOUNTS> my_taxon_counts;
          uint64_t my_total_classified = 0;
          kraken_output_ss.str("");
          classified_output_ss.str("");
          unclassified_output_ss.str("");
          my_total_classified += classify_sequence(
                                                   work_unit[j], kraken_output_ss, classified_output_ss,
                                                   unclassified_output_ss, my_taxon_counts);
          #ifdef _OPENMP
          #pragma omp ordered
          #endif
          {
            total_classified += my_total_classified;
            for (auto it = my_taxon_counts.begin(); it != my_taxon_counts.end();
                 ++it) {
              taxon_counts[it->first] += std::move(it->second);
            }

            if (Print_kraken)
              (*kraken_output) << kraken_output_ss.str();
            if (Print_classified)
              (*classified_output) << classified_output_ss.str();
            if (Print_unclassified)
              (*unclassified_output) << unclassified_output_ss.str();
          }
        }
      }
      total_sequences += work_unit.size();
      total_bases += total_nt;
      // if (Print_Progress && total_sequences % 100000 < work_unit.size())
      if (Print_Progress) {
        fprintf(stderr, "\r Processed %llu sequences (%.2f%% classified)",
                total_sequences, total_classified * 100.0 / total_sequences);
      }
    }
  } else {
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
      vector<DNASequence> work_unit;
      ostringstream kraken_output_ss, classified_output_ss,
        unclassified_output_ss;

      while (reader->is_valid()) {
        work_unit.clear();
        size_t total_nt = 0;

        #ifdef _OPENMP
        #pragma omp critical(get_input)
        #endif
        {
          while (total_nt < Work_unit_size) {
            dna = reader->next_sequence();
            if (!reader->is_valid())
              break;
            work_unit.push_back(dna);
            total_nt += dna.seq.size();
          }
        }
        if (total_nt == 0)
          break;

        unordered_map<uint32_t, READCOUNTS> my_taxon_counts;
        uint64_t my_total_classified = 0;
        kraken_output_ss.str("");
        classified_output_ss.str("");
        unclassified_output_ss.str("");
        for (size_t j = 0; j < work_unit.size(); j++) {
          my_total_classified += classify_sequence(
                                                   work_unit[j], kraken_output_ss, classified_output_ss,
                                                   unclassified_output_ss, my_taxon_counts);
        }

        #ifdef _OPENMP
        #pragma omp critical(write_output)
        #endif
        {
          total_classified += my_total_classified;
          for (auto it = my_taxon_counts.begin(); it != my_taxon_counts.end();
               ++it) {
            taxon_counts[it->first] += std::move(it->second);
          }

          if (Print_kraken)
            (*kraken_output) << kraken_output_ss.str();
          if (Print_classified)
            (*classified_output) << classified_output_ss.str();
          if (Print_unclassified)
            (*unclassified_output) << unclassified_output_ss.str();
          total_sequences += work_unit.size();
          total_bases += total_nt;
          // if (Print_Progress && total_sequences % 100000 < work_unit.size())
          if (Print_Progress) {
            fprintf(stderr, "\r Processed %llu sequences (%.2f%% classified)",
                    total_sequences,
                    total_classified * 100.0 / total_sequences);
          }
        }
      }
    } // end parallel section
  }

  delete reader;
}


inline void print_sequence(ostringstream* oss_ptr, const DNASequence& dna) {
      if (Fastq_input) {
        (*oss_ptr) << "@" << dna.header_line << endl
            << dna.seq << endl
            << "+" << endl
            << dna.quals << endl;
      }
      else {
        (*oss_ptr) << ">" << dna.header_line << endl
            << dna.seq << endl;
      }
}

/*
inline
void append_hitlist_string(string& hitlist_string, uint32_t& last_taxon, uint32_t& last_counter, uint32_t current_taxon) {
  if (last_taxon == current_taxon) {
    ++last_counter;
  } else {
    if (last_counter > 0) {
      if (last_taxon == ambig_taxon) {
        hitlist_string += "A:" + std::to_string(last_counter) + ' ';
      } else {
        hitlist_string += std::to_string(last_taxon) + ':' + std::to_string(last_counter) + ' ';
      }
    }
    last_counter = 1;
    last_taxon = current_taxon;
  }
}
*/

string hitlist_string(const vector<uint32_t> &taxa, const vector<uint8_t> &ambig)
{
  int64_t last_code;
  int code_count = 1;
  ostringstream hitlist;

  if (ambig[0])   { last_code = -1; }
  else            { last_code = taxa[0]; }

  for (size_t i = 1; i < taxa.size(); i++) {
    int64_t code;
    if (ambig[i]) { code = -1; }
    else          { code = taxa[i]; }

    if (code == last_code) {
      code_count++;
    }
    else {
      if (last_code >= 0) {
        hitlist << last_code << ":" << code_count << " ";
      }
      else {
        hitlist << "A:" << code_count << " ";
      }
      code_count = 1;
      last_code = code;
    }
  }
  if (last_code >= 0) {
    hitlist << last_code << ":" << code_count;
  }
  else {
    hitlist << "A:" << code_count;
  }
  return hitlist.str();
}

/*
string hitlist_string_depr(const vector<uint32_t> &taxa)
{
  uint32_t last_code = taxa[0];
  int code_count = 1;
  ostringstream hitlist;

  for (size_t i = 1; i < taxa.size(); i++) {
    uint32_t code = taxa[i];

    if (code == last_code) {
      code_count++;
    }
    else {
      if (last_code >= 0) {
        hitlist << last_code << ":" << code_count << " ";
      }
      else {
        hitlist << "A:" << code_count << " ";
      }
      code_count = 1;
      last_code = code;
    }
  }
  if (last_code == -1) {
    hitlist << "A:" << code_count;
  }
  else {
    hitlist << last_code << ":" << code_count;
  }
  return hitlist.str();
}
*/

bool classify_sequence(DNASequence &dna, ostringstream &koss,
                       ostringstream &coss, ostringstream &uoss,
                       unordered_map<uint32_t, READCOUNTS>& my_taxon_counts) {
  vector<uint32_t> taxa;
  vector<uint8_t> ambig_list;
  unordered_map<uint32_t, uint32_t> hit_counts;
  uint64_t *kmer_ptr;
  uint32_t taxon = 0;
  uint32_t hits = 0;  // only maintained if in quick mode

  //string hitlist_string;
  //uint32_t last_taxon;
  //uint32_t last_counter;

  vector<db_status> db_statuses(KrakenDatabases.size());

  if (dna.seq.size() >= KrakenDatabases[0]->get_k()) {
    size_t n_kmers = dna.seq.size()-KrakenDatabases[0]->get_k()+1;
    taxa.reserve(n_kmers);
    ambig_list.reserve(n_kmers);
    KmerScanner scanner(dna.seq);
    while ((kmer_ptr = scanner.next_kmer()) != NULL) {
      taxon = 0;
      if (scanner.ambig_kmer()) {
        //append_hitlist_string(hitlist_string, last_taxon, last_counter, ambig_taxon);
        ambig_list.push_back(1);
      }
      else {
        uint64_t cannonical_kmer = KrakenDatabases[0]->canonical_representation(*kmer_ptr);
        ambig_list.push_back(0);
        // go through multiple databases to map k-mer
        for (size_t i=0; i<KrakenDatabases.size(); ++i) {
          uint32_t* val_ptr = KrakenDatabases[i]->kmer_query(
            cannonical_kmer, &db_statuses[i].current_bin_key,
            &db_statuses[i].current_min_pos, &db_statuses[i].current_max_pos);
          if (val_ptr) {
            taxon = *val_ptr;
            break;
          }
        }

        // cerr << "taxon for " << *kmer_ptr << " is " << taxon << endl;
        my_taxon_counts[taxon].add_kmer(cannonical_kmer);

        if (taxon) {
          hit_counts[taxon]++;
          if (Quick_mode && ++hits >= Minimum_hit_count)
            break;
        }
      }
      taxa.push_back(taxon);
      //append_hitlist_string(hitlist_string, last_taxon, last_counter, taxon);
    }
  }

  uint32_t call = 0;
  if (Map_UIDs) {
    if (Quick_mode) {
      cerr << "Quick mode not available when mapping UIDs" << endl;
      exit(1);
    } else {
      call = resolve_uids3(hit_counts, Parent_map, Uid_dict,
        UID_to_TaxID_map_file.ptr(), UID_to_TaxID_map_file.size());
    }
  } else {
    if (Quick_mode)
      call = hits >= Minimum_hit_count ? taxon : 0;
    else
      call = resolve_tree(hit_counts, Parent_map);
  }

  my_taxon_counts[call].incrementReadCount();

  if (Print_unclassified && !call) 
    print_sequence(&uoss, dna);

  if (Print_classified && call)
    print_sequence(&coss, dna);


  if (! Print_kraken)
    return call;

  if (call) {
    koss << "C\t";
  }
  else {
    if (Only_classified_kraken_output)
      return false;
    koss << "U\t";
  }
  koss << dna.id << '\t' << call << '\t' << dna.seq.size() << '\t';

  if (Quick_mode) {
    koss << "Q:" << hits;
  }
  else {
    if (taxa.empty())
      koss << "0:0";
    else
      koss << hitlist_string(taxa, ambig_list);
    //if (hitlist_string.empty() && last_counter == 0)
    //  koss << "0:0";
    //else {
    //  koss << hitlist_string
    //       << (last_taxon == ambig_taxon? "A" :  std::to_string(last_taxon))
    //       << ':' << std::to_string(last_counter);
    //}
  }

  if (Print_sequence)
      koss << "\t" << dna.seq;

  koss << "\n";
  return call;
}

set<uint32_t> get_ancestry(uint32_t taxon) {
  set<uint32_t> path;

  while (taxon > 0) {
    path.insert(taxon);
    taxon = Parent_map[taxon];
  }
  return path;
}

void parse_command_line(int argc, char **argv) {
  int opt;
  long long sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "d:i:D:t:u:n:m:o:qOfcC:U:Ma:r:sI:p:")) !=
         -1) {
    switch (opt) {
      case 'd':
        DB_filenames.push_back(optarg);
        break;
      case 'i':
        Index_filenames.push_back(optarg);
        break;
      case 'D':
        Onefile_DB_filenames.push_back(optarg);
        break;
      case 'l':
        Lock_DB = true;
        break;
      case 't':
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
      case 'p':
        HLL_PRECISION = stoi(optarg);
        break;
      case 'q':
        Quick_mode = true;
        break;
      case 'm':
        sig = atoll(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive minimum hit count");
        Minimum_hit_count = sig;
        break;
      case 'f':
        Fastq_input = true;
        break;
      case 'c':
        Only_classified_kraken_output = true;
        break;
      case 'C':
        Print_classified = true;
        // Classified_output_file = optarg;
        Classified_output_files.push_back(optarg);
        break;
      case 'U':
        Print_unclassified = true;
        // Unclassified_output_file = optarg;
        Unclassified_output_files.push_back(optarg);
        break;
      case 'o':
        // Kraken_output_file = optarg;
        Kraken_output_files.push_back(optarg);
        break;
      case 'r':
        // Report_output_file = optarg;
        Report_output_files.push_back(optarg);
        break;
      case 's':
        Print_sequence = true;
        break;
      case 'a':
        TaxDB_file = optarg;
        break;
      case 'u':
        sig = atoll(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive work unit size");
        Work_unit_size = sig;
        break;
      case 'M':
        Populate_memory = true;
        break;
      case 'I':
        UID_to_TaxID_map_filename = optarg;
        Map_UIDs = true;
        break;
      default:
        usage();
        break;
    case 'O':
      Ordered_mode = true;
      break;
    }
  }

  if (Onefile_DB_filenames.empty()) {
    if (DB_filenames.empty()) {
      cerr << "Missing mandatory option -d" << endl;
      usage();
    }
    if (Index_filenames.empty()) {
      cerr << "Missing mandatory option -i" << endl;
      usage();
    }
    if (optind == argc && !Populate_memory) {
      cerr << "No sequence data files specified" << endl;
    }
  }
}

void usage(int exit_code) {
  cerr << "Usage: classify [options] <fasta/fastq file(s)>" << endl
       << endl
       << "Options: (*mandatory)" << endl
       << "* -d filename      Kraken DB filename" << endl
       << "* -i filename      Kraken DB index filename" << endl
       << "  -D filename      Kraken DB stored in a single file" << endl
       << "  -o filename      Output file for Kraken output" << endl
       << "  -r filename      Output file for Kraken report output" << endl
       << "  -a filename      TaxDB" << endl
       << "  -I filename      UID to TaxId map" << endl
       << "  -l               Memory lock DB files" << endl
       << "  -p #             Precision for unique k-mer counting, between 10 "
    "and 18"
       << endl
       << "  -t #             Number of threads" << endl
       << "  -u #             Thread work unit size (in bp)" << endl
       << "  -O               Order output matching input" << endl
       << "  -q               Quick operation" << endl
       << "  -m #             Minimum hit count (ignored w/o -q)" << endl
       << "  -C filename      Print classified sequences" << endl
       << "  -U filename      Print unclassified sequences" << endl
       << "  -f               Input is in FASTQ format" << endl
       << "  -c               Only include classified reads in output" << endl
       << "  -M               Preload database files" << endl
       << "  -s               Print read sequence in Kraken output" << endl
       << "  -h               Print this message" << endl
       << endl
       << "At least one FASTA or FASTQ file must be specified." << endl
       << "Kraken output is to standard output by default." << endl;
  exit(exit_code);
}
