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
#include "seqreader.hpp"
#include "readcounts.hpp"
#include "taxdb.hpp"
#include "gzstream.h"
#include "uid_mapping.hpp"
#include <sstream>
#include <inttypes.h>
#include <cassert>

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


#ifndef _OPENMP
  int omp_get_thread_num() { return 0; }
#endif

void parse_command_line(int argc, char **argv);
void usage(int exit_code=EX_USAGE);
void process_file(char *filename);
void process_file_with_db_chunk(char *filename);
void classify_sequence_with_db_chunk(std::pair<DNASequence, uint32_t> & seq, std::fstream & fp, const uint32_t db_chunk_id, const uint32_t db_id);
bool classify_sequence(DNASequence &dna, ostringstream &koss,
                       ostringstream &coss, ostringstream &uoss,
                       unordered_map<uint32_t, READCOUNTS>&);
inline void print_sequence(ostream* oss_ptr, const DNASequence& dna);
string hitlist_string(const vector<uint32_t> &taxa, const vector<char>& ambig_list);


set<uint32_t> get_ancestry(uint32_t taxon);
void report_stats(struct timeval time1, struct timeval time2);
double get_seconds(struct timeval time1, struct timeval time2);
unordered_map<uint32_t, READCOUNTS> taxon_counts; // stats per taxon

int Num_threads = 1;
vector<string> DB_filenames;
vector<string> Index_filenames;
bool Quick_mode = false;
bool Fastq_input = false;
bool Print_classified = false;
bool Print_unclassified = false;
bool Print_kraken = true;
bool Print_kraken_report = false;
bool Populate_memory = false;
uint64_t Populate_memory_size = 0;
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
unordered_map<uint32_t, vector<uint32_t> > Uid_dict;
string Classified_output_file, Unclassified_output_file, Kraken_output_file, Report_output_file, TaxDB_file;
ostream *Classified_output;
ostream *Unclassified_output;
ostream *Kraken_output;
ostream *Report_output;
vector<ofstream*> Open_fstreams;
vector<ogzstream*> Open_gzstreams;
size_t Work_unit_size = DEF_WORK_UNIT_SIZE;
TaxonomyDB<uint32_t> taxdb;
static vector<KrakenDB*> KrakenDatabases (DB_filenames.size());

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

inline bool ends_with(std::string const & value, std::string const & ending)
{
        if (ending.size() > value.size()) return false;
            return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

ostream* cout_or_file(string file, bool append = false) {
    if (file == "-")
      return &cout;

    if (ends_with(file, ".gz")) {
      ogzstream* ogzs = new ogzstream(file.c_str());
      ogzs->exceptions( ifstream::failbit | ifstream::badbit );
      Open_gzstreams.push_back(ogzs);
      return ogzs;
    } else {
      ofstream* ofs = append? new ofstream(file.c_str(), std::ofstream::app) : new ofstream(file.c_str());
      ofs->exceptions( ifstream::failbit | ifstream::badbit );
      Open_fstreams.push_back(ofs);
      return ofs;
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

  if (Populate_memory && Populate_memory_size == 0)
    cerr << "Loading database(s)... " << endl;

  static vector<QuickFile> idx_files (DB_filenames.size());
  static vector<QuickFile> db_files (DB_filenames.size());
  static vector<KrakenDBIndex> db_indices (DB_filenames.size());


  // TODO: Check DB_filenames and Index_filesnames have the same length
  for (size_t i=0; i < DB_filenames.size(); ++i) {
    cerr << " Database " << DB_filenames[i] << endl;
    db_files[i].open_file(DB_filenames[i]);

    KrakenDatabases.push_back(new KrakenDB(db_files[i].ptr()));
    idx_files[i].open_file(Index_filenames[i]);
    // NOTE: we switched the order, i.e., we are creating the objects before loading everything into main memory
    db_indices[i] = KrakenDBIndex(idx_files[i].ptr());
    KrakenDatabases[i]->set_index(&db_indices[i]);

    if (Populate_memory && Populate_memory_size == 0) // only when no chunk size is passed!
    {
      db_files[i].load_file();
      idx_files[i].load_file();
    }
    else if (Populate_memory && Populate_memory_size > 0)
    {
      KrakenDatabases[i]->prepare_chunking(Populate_memory_size);
    }
  }

  // Check all databases have the same k
  uint8_t kmer_size = KrakenDatabases[0]->get_k();
  for (size_t i = 1; i < KrakenDatabases.size(); ++i) {
    uint8_t kmer_size_i = KrakenDatabases[i]->get_k();
    if (kmer_size_i != kmer_size) {
      fprintf(stderr, "Different k-mer sizes in databases 1 and %lu: %i vs %i!\n", i+1, (int)kmer_size, (int)kmer_size_i);
      exit(1);
    }
  };
  KmerScanner::set_k(kmer_size);

  if (Populate_memory && Populate_memory_size == 0)
    cerr << "\ncomplete." << endl;


  if (!TaxDB_file.empty()) {
      taxdb = TaxonomyDB<uint32_t>(TaxDB_file, false);
      Parent_map = taxdb.getParentMap();
  } else {
      cerr << "TaxDB argument is required!" << endl;
      return 1;
  }

  if (Print_classified) {
    Classified_output = cout_or_file(Classified_output_file);
  }

  if (Print_unclassified) {
    Unclassified_output = cout_or_file(Unclassified_output_file);
  }

  if (! Kraken_output_file.empty()) {
    if (Kraken_output_file == "off" || Kraken_output_file == "-") {
      Print_kraken = false;
    //else if (Kraken_output_file == "-") {
    //  Kraken_output = &cout;
    } else {
      cerr << "Writing Kraken output to " << Kraken_output_file << endl;
      Kraken_output = cout_or_file(Kraken_output_file);
    }
  } else {
    Kraken_output = &cout;
  }

  //cerr << "Print_kraken: " << Print_kraken << "; Print_kraken_report: " << Print_kraken_report << "; k: " << uint32_t(KrakenDatabases[0]->get_k()) << endl;

  struct timeval tv1, tv2;
  gettimeofday(&tv1, NULL);
  for (int i = optind; i < argc; i++) {
    if (Populate_memory && Populate_memory_size > 0)
      process_file_with_db_chunk(argv[i]);
    else
      process_file(argv[i]);
  }
  gettimeofday(&tv2, NULL);

  report_stats(tv1, tv2);

  if (!Report_output_file.empty() && Report_output_file != "off") {
    gettimeofday(&tv1, NULL);
    std::cerr << "Writing report file to " << Report_output_file <<"  ..\n";
    for (size_t i = 0; i < DB_filenames.size(); ++i) {
      const auto fname = DB_filenames[i] + ".counts";
      ifstream ifs(fname);
      bool counts_file_gd = false;
      if (ifs.good()) {
        if (ifs.peek() == std::ifstream::traits_type::eof()) {
          cerr << "Kmer counts file is empty - trying to regenerate ..." << endl;
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
     Report_output = cout_or_file(Report_output_file, true);
  
    TaxReport<uint32_t,READCOUNTS> rep = TaxReport<uint32_t, READCOUNTS>(*Report_output, taxdb, taxon_counts, false);
    if (HLL_PRECISION > 0) {
      if (full_report) {
        rep.setReportCols(vector<string> { 
          "%",
          "reads", 
          "taxReads",
          "kmers",
          "taxKmers",
          "kmersDB",
          "taxKmersDB",
          "dup",
          "cov", 
          "taxID", 
          "rank", 
          "taxName"});
      } else {
        rep.setReportCols(vector<string> { 
          "%",
          "reads", 
          "taxReads",
          "kmers",
          "dup",
          "cov", 
          "taxID", 
          "rank", 
          "taxName"});
      }
    } else {
      rep.setReportCols(vector<string> { 
        "%",
        "reads", 
        "taxReads",
        "taxID", 
        "rank", 
        "taxName"});
    }
    rep.printReport("kraken");
    gettimeofday(&tv2, NULL);
    fprintf(stderr, "Report finished in %.3f seconds.\n", get_seconds(tv1,tv2));
  }
  cerr << "Finishing up ...";

  for (size_t i = 0; i < Open_fstreams.size(); ++i) {
    ofstream* ofs = Open_fstreams[i];
    ofs->close();
  }

  for (size_t i = 0; i < Open_gzstreams.size(); ++i) {
    ogzstream* ogzs = Open_gzstreams[i];
    ogzs->close();
  }

  for (size_t i=0; i < KrakenDatabases.size(); ++i) {
    delete KrakenDatabases[i];
  }

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

void merge_intermediate_results_by_workers(const bool first_intermediate_output, const std::string & tmp_file_name) {
  const std::string filename_merged_summary = tmp_file_name;

  FILE *fp_prev_merged_summary = NULL;
  std::string filename_prev_merged_summary;
  if (!first_intermediate_output)
  {
    filename_prev_merged_summary = filename_merged_summary + ".prev";
    rename(filename_merged_summary.c_str(), filename_prev_merged_summary.c_str());
    fp_prev_merged_summary = fopen(filename_prev_merged_summary.c_str(), "rb");
  }

  std::fstream fp_merged_summary(filename_merged_summary, std::fstream::out | std::fstream::binary | std::fstream::trunc);
  fp_merged_summary.exceptions(std::fstream::badbit);

  std::vector<FILE*> worker_files(Num_threads);
  std::vector<uint32_t> next_read_id(Num_threads);
  for (int worker_id = 0; worker_id < Num_threads; ++worker_id)
  {
    const std::string worker_filename = filename_merged_summary + "." + std::to_string(worker_id);
    worker_files[worker_id] = fopen(worker_filename.c_str(), "rb");

    if (fread(&next_read_id[worker_id], sizeof(next_read_id[worker_id]), 1, worker_files[worker_id]) != 1)
      next_read_id[worker_id] = 0; // 0 indicates that there is no read left in this file (reads are indexed starting at 1)
  }

  std::vector<uint32_t> taxa, taxa_prev_summary;
  uint32_t seq_idx = 1;
  while (true)
  {
    // determine worker to get this read info from
    int worker_to_retrieve_from = -1;
    for (int worker_id = 0; worker_id < Num_threads; ++worker_id) {
      if (seq_idx == next_read_id[worker_id]) {
        worker_to_retrieve_from = worker_id;
        break;
      }
    }
    if (worker_to_retrieve_from == -1)
      break; // finished processing all files

    // retrieve from that worker and load next sequence name
    uint32_t taxa_size;
    fread(&taxa_size, sizeof(uint32_t), 1, worker_files[worker_to_retrieve_from]);
    taxa.resize(taxa_size);
    fread(&taxa[0], sizeof(uint32_t), taxa_size, worker_files[worker_to_retrieve_from]);

    if (!first_intermediate_output)
    {
      uint32_t taxa_prev_size;
      fread(&taxa_prev_size, sizeof(uint32_t), 1, fp_prev_merged_summary);
      assert(taxa_size == taxa_prev_size);
      taxa_prev_summary.resize(taxa_prev_size);
      fread(&taxa_prev_summary[0], sizeof(uint32_t), taxa_prev_size, fp_prev_merged_summary);

      for (uint32_t i = 0; i < taxa_prev_summary.size(); ++i)
      {
        assert(!(taxa[i] && taxa_prev_summary[i])); // a kmer cannot be found in more than one database chunk
        if (taxa_prev_summary[i])
        {
          taxa[i] = taxa_prev_summary[i]; // if previous chunk/database got a match, keep the tax info
        }
      }
    }

    try {
      fp_merged_summary.write((char*) &taxa_size, 1 * sizeof(uint32_t));
      fp_merged_summary.write((char*) &taxa[0], taxa_size * sizeof(uint32_t));
    }
    catch (const std::fstream::failure& e) {
      printf("ERROR: Could not write to temporary file: %s\n", e.what());
      exit(1);
    }

    // load the next read id from that file
    if (fread(&next_read_id[worker_to_retrieve_from], sizeof(next_read_id[worker_to_retrieve_from]), 1, worker_files[worker_to_retrieve_from]) != 1)
      next_read_id[worker_to_retrieve_from] = 0; // 0 indicates that there is no read left in this file (reads are indexed starting at 1)

    ++seq_idx;
  }

  for (int worker_id = 0; worker_id < Num_threads; ++worker_id)
  {
    fclose(worker_files[worker_id]);
    const std::string worker_filename = filename_merged_summary + "." + std::to_string(worker_id);
    unlink(worker_filename.c_str());
  }

  if (!first_intermediate_output)
  {
    fclose(fp_prev_merged_summary);
    unlink(filename_prev_merged_summary.c_str());
  }

  fp_merged_summary.close();
}

void process_file(char *filename) {
  string file_str(filename);
  DNASequenceReader *reader;
  DNASequence dna;

  if (Fastq_input)
    reader = new FastqReader(file_str);
  else
    reader = new FastaReader(file_str);

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
    vector<DNASequence> work_unit;
    ostringstream kraken_output_ss, classified_output_ss, unclassified_output_ss;

    while (reader->is_valid()) {
      work_unit.clear();
      size_t total_nt = 0;

#ifdef _OPENMP
      #pragma omp critical(get_input)
#endif
      {
        while (total_nt < Work_unit_size) {
          dna = reader->next_sequence();
          if (! reader->is_valid())
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
        my_total_classified += 
            classify_sequence( work_unit[j], kraken_output_ss,
                           classified_output_ss, unclassified_output_ss,
                           my_taxon_counts);
      }
 
#ifdef _OPENMP
      #pragma omp critical(write_output)
#endif
      {
        total_classified += my_total_classified;
        for (auto it = my_taxon_counts.begin(); it != my_taxon_counts.end(); ++it) {
          taxon_counts[it->first] += std::move(it->second);
        }

        if (Print_kraken)
          (*Kraken_output) << kraken_output_ss.str();
        if (Print_classified)
          (*Classified_output) << classified_output_ss.str();
        if (Print_unclassified)
          (*Unclassified_output) << unclassified_output_ss.str();
        total_sequences += work_unit.size();
        total_bases += total_nt;
        //if (Print_Progress && total_sequences % 100000 < work_unit.size()) 
        if (Print_Progress) {  
          fprintf(stderr, "\r Processed %llu sequences (%.2f%% classified)",
                          total_sequences, total_classified * 100.0 / total_sequences);
        }
      }
    }
  }  // end parallel section

  delete reader;
}

void process_file_with_db_chunk(char *filename) {
  string file_str(filename);
  DNASequence dna;

  const std::string tmp_file_name = std::tmpnam(nullptr);

  // iterate over databases
  bool first_intermediate_output = true;
  for (size_t i=0; i<KrakenDatabases.size(); ++i)
  {
    for (uint32_t db_chunk_id = 0; db_chunk_id < KrakenDatabases[0]->chunks(); ++db_chunk_id)
    {
      total_sequences = 0;
      total_bases = 0;
      KrakenDatabases[i]->load_chunk(db_chunk_id);

      DNASequenceReader *reader;
      if (Fastq_input)
        reader = new FastqReader(file_str);
      else
        reader = new FastaReader(file_str);
      uint32_t seq_idx = 1;

#ifdef _OPENMP
      #pragma omp parallel
#endif
      {
        vector <std::pair<DNASequence, uint32_t> > work_unit;
        const int worker_id = omp_get_thread_num();
        const std::string worker_filename = tmp_file_name + "." + std::to_string(worker_id);

        std::fstream fp(worker_filename, std::fstream::out | std::fstream::binary | std::fstream::trunc);
        fp.exceptions(std::fstream::badbit);

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
              work_unit.emplace_back(dna, seq_idx);
              total_nt += dna.seq.size();
              ++seq_idx;
            }
          }
          if (total_nt == 0)
            break;

          for (size_t j = 0; j < work_unit.size(); j++) {
            classify_sequence_with_db_chunk(work_unit[j], fp, db_chunk_id, 0);
          }

#ifdef _OPENMP
          #pragma omp critical(progress)
#endif
          {
            total_sequences += work_unit.size();
            total_bases += total_nt;
            if (Print_Progress) {
              fprintf(stderr, "\r Processed %llu sequences (database chunk %" PRIu32" of %" PRIu32 ")",
                      total_sequences, db_chunk_id + 1, KrakenDatabases[i]->chunks());
            }
          }
        }
        fp.close();
      }  // end parallel section

      delete reader;
      merge_intermediate_results_by_workers(first_intermediate_output, tmp_file_name);
      first_intermediate_output = false;
    }
  }

  fprintf(stderr, "\r Processed %llu sequences\n", total_sequences);

  // merge intermediate results of all chunks and classify
  // TODO: parallelize this (need to buffer reads), we also do not care about the final output order
  total_sequences = 0;
  total_classified = 0;

  DNASequenceReader *reader;
  if (Fastq_input)
    reader = new FastqReader(file_str);
  else
    reader = new FastaReader(file_str);

  unordered_map<uint32_t, uint32_t> hit_counts;
  vector<uint32_t> taxa;
  vector<char> ambig_list;

  const std::string taxa_summary_filename = tmp_file_name;
  FILE* fp_taxa_summary = fopen(taxa_summary_filename.c_str(), "rb");

  while (reader->is_valid()) {
    hit_counts.clear();
    taxa.clear();
    ambig_list.clear();

    dna = reader->next_sequence();
    if (!reader->is_valid())
      break;

    // merge results from chunks into 'taxa'
    uint32_t hits = 0;

    // get number of elements (taxa.size())
    uint32_t taxa_size;
    fread(&taxa_size, sizeof(taxa_size), 1, fp_taxa_summary);

    taxa.resize(taxa_size);
    fread(&taxa[0], sizeof(uint32_t), taxa_size, fp_taxa_summary);

    for (uint32_t i = 0; i < taxa.size(); ++i)
    {
      if (taxa[i])
      {
        ++hit_counts[taxa[i]];
        // TODO: stop querying this read with other databases
        if (Quick_mode && ++hits >= Minimum_hit_count)
          goto quick_mode_call;
      }
    }

quick_mode_call:
    uint64_t *kmer_ptr;
    uint32_t taxon = 0;
    if (dna.seq.size() >= KrakenDatabases[0]->get_k()) {
      KmerScanner scanner(dna.seq);
      uint32_t taxa_idx = 0;
      while ((kmer_ptr = scanner.next_kmer()) != NULL) {
        if (scanner.ambig_kmer()) {
          ambig_list.push_back(1);
        } else {
          ambig_list.push_back(0);
          uint64_t cannonical_kmer = KrakenDatabases[0]->canonical_representation(*kmer_ptr);
          taxon = taxa[taxa_idx];
          taxon_counts[taxon].add_kmer(cannonical_kmer);
        }
        ++taxa_idx;
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

    total_classified += (call != 0);
    taxon_counts[call].incrementReadCount();

    if (Print_unclassified && !call)
      print_sequence(Unclassified_output, dna);

    if (Print_classified && call)
      print_sequence(Classified_output, dna);

    if (!Print_kraken) {
      goto next_read; // return call;
    }

    if (call) {
      (*Kraken_output) << "C\t";
    }
    else {
      if (Only_classified_kraken_output) {
        goto next_read; // return false;
      }
      (*Kraken_output) << "U\t";
    }
    (*Kraken_output) << dna.id << '\t' << call << '\t' << dna.seq.size() << '\t';

    if (Quick_mode) {
      (*Kraken_output) << "Q:" << hits;
    }
    else {
      if (taxa.empty())
        (*Kraken_output) << "0:0";
      else
        (*Kraken_output) << hitlist_string(taxa, ambig_list);
    }

    if (Print_sequence)
      (*Kraken_output) << "\t" << dna.seq;

    (*Kraken_output) << "\n";
    // return call;
next_read:
    ++total_sequences;
    fprintf(stderr, "\r Processed %llu sequences (%.2f%% classified)",
            total_sequences, total_classified * 100.0 / total_sequences);
    (void)1;
  }

  fclose(fp_taxa_summary);
  unlink(taxa_summary_filename.c_str());

  delete reader;
}


inline void print_sequence(ostream* oss_ptr, const DNASequence& dna) {
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

string hitlist_string(const vector<uint32_t> &taxa, const vector<char> &ambig)
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
  vector<char> ambig_list;
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

void classify_sequence_with_db_chunk(std::pair<DNASequence, uint32_t> & seq, std::fstream & fp, const uint32_t db_chunk_id, const uint32_t db_id) {
  vector<uint32_t> taxa;
  uint64_t *kmer_ptr;
  uint32_t taxon;

  auto & dna = seq.first;
  const uint32_t & seq_idx = seq.second;

  vector<db_status> db_statuses(KrakenDatabases.size());

  if (dna.seq.size() >= KrakenDatabases[0]->get_k()) {
    size_t n_kmers = dna.seq.size()-KrakenDatabases[0]->get_k()+1;
    taxa.reserve(n_kmers);
    KmerScanner scanner(dna.seq);
    while ((kmer_ptr = scanner.next_kmer()) != NULL) {
      taxon = 0;
      if (!scanner.ambig_kmer()) {
        uint64_t cannonical_kmer = KrakenDatabases[db_id]->canonical_representation(*kmer_ptr);
        const uint64_t minimizer = KrakenDatabases[db_id]->bin_key(cannonical_kmer); // TODO: inefficient because minimizer will also be computed in kmer_query()

        if (KrakenDatabases[db_id]->is_minimizer_in_chunk(minimizer, db_chunk_id)) {
          uint32_t* val_ptr = KrakenDatabases[db_id]->kmer_query_with_db_chunks(
                  cannonical_kmer, &db_statuses[db_id].current_bin_key,
                  &db_statuses[db_id].current_min_pos, &db_statuses[db_id].current_max_pos);
          if (val_ptr)
            taxon = *val_ptr;
        }
      }
      taxa.push_back(taxon);
    }
  }

  try {
    const uint32_t taxa_size = taxa.size();
    fp.write((char*) &seq_idx, 1 * sizeof(uint32_t)); // seq_idx
    fp.write((char*) &taxa_size, 1 * sizeof(uint32_t)); // number of elements
    fp.write((char*) &taxa[0], taxa_size * sizeof(uint32_t)); // elements
  }
  catch (const std::fstream::failure& e) {
    printf("ERROR: Could not write to temporary file: %s\n", e.what());
    exit(1);
  }
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
  while ((opt = getopt(argc, argv, "d:i:t:u:n:m:o:qfcC:U:Ma:r:sI:p:x:")) != -1) {
    switch (opt) {
      case 'd' :
        DB_filenames.push_back(optarg);
        break;
      case 'i' :
        Index_filenames.push_back(optarg);
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
      case 'p' :
        HLL_PRECISION = stoi(optarg);
        break;
      case 'q' :
        Quick_mode = true;
        break;
      case 'm' :
        sig = atoll(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive minimum hit count");
        Minimum_hit_count = sig;
        break;
      case 'f' :
        Fastq_input = true;
        break;
      case 'c' :
        Only_classified_kraken_output = true;
        break;
      case 'C' :
        Print_classified = true;
        Classified_output_file = optarg;
        break;
      case 'U' :
        Print_unclassified = true;
        Unclassified_output_file = optarg;
        break;
      case 'o' :
        Kraken_output_file = optarg;
        break;
      case 'r' :
        Report_output_file = optarg;
        break;
      case 's' :
        Print_sequence = true;
        break;
      case 'a' :
        TaxDB_file = optarg;
        break;
      case 'u' :
        sig = atoll(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive work unit size");
        Work_unit_size = sig;
        break;
      case 'M' :
        Populate_memory = true;
        break;
      case 'x' :
        Populate_memory = true;
        Populate_memory_size = parse_human_readable_size(optarg); // strtoull(optarg, NULL, 0);
        break;
      case 'I' :
        UID_to_TaxID_map_filename = optarg;
        Map_UIDs = true;
        break;
      default:
        usage();
        break;
    }
  }

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

void usage(int exit_code) {
  cerr << "Usage: classify [options] <fasta/fastq file(s)>" << endl
       << endl
       << "Options: (*mandatory)" << endl
       << "* -d filename      Kraken DB filename" << endl
       << "* -i filename      Kraken DB index filename" << endl
       << "  -o filename      Output file for Kraken output" << endl
       << "  -r filename      Output file for Kraken report output" << endl
       << "  -a filename      TaxDB" << endl
       << "  -I filename      UID to TaxId map" << endl
       << "  -p #             Precision for unique k-mer counting, between 10 and 18" << endl
       << "  -t #             Number of threads" << endl
       << "  -u #             Thread work unit size (in bp)" << endl
       << "  -q               Quick operation" << endl
       << "  -m #             Minimum hit count (ignored w/o -q)" << endl
       << "  -C filename      Print classified sequences" << endl
       << "  -U filename      Print unclassified sequences" << endl
       << "  -f               Input is in FASTQ format" << endl
       << "  -c               Only include classified reads in output" << endl
       << "  -M               Preload database files" << endl
       << "  -x size          Preload database files using x amount of RAM (e.g. 10G)" << endl
       << "  -s               Print read sequence in Kraken output" << endl
       << "  -h               Print this message" << endl
       << endl
       << "At least one FASTA or FASTQ file must be specified." << endl
       << "Kraken output is to standard output by default." << endl;
  exit(exit_code);
}
