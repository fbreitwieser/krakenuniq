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
#include "krakendb.hpp"
#include "krakenutil.hpp"
#include "quickfile.hpp"
#include "seqreader.hpp"
#include "hyperloglogplus.h"
#include "taxdb.h"
#include "gzstream.h"

const size_t DEF_WORK_UNIT_SIZE = 500000;
int New_taxid_start = 1000000000;

using namespace std;
using namespace kraken;

void parse_command_line(int argc, char **argv);
void usage(int exit_code=EX_USAGE);
void process_file(char *filename);
bool classify_sequence(DNASequence &dna, ostringstream &koss,
                       ostringstream &coss, ostringstream &uoss,
                       unordered_map<uint32_t, ReadCounts>&);
string hitlist_string(vector<uint32_t> &taxa, vector<uint8_t> &ambig);
set<uint32_t> get_ancestry(uint32_t taxon);
void report_stats(struct timeval time1, struct timeval time2);
unordered_map<uint32_t, ReadCounts> taxon_counts; // stats per taxon

int Num_threads = 1;
vector<string> DB_filenames;
vector<string> Index_filenames;
bool Quick_mode = false;
bool Fastq_input = false;
bool Print_classified = false;
bool Print_unclassified = false;
bool Print_kraken = true;
bool Print_kraken_report = true;
bool Populate_memory = false;
bool Only_classified_kraken_output = false;
bool Print_sequence = false;
bool Print_Progress = true;
uint32_t Minimum_hit_count = 1;
unordered_map<uint32_t, uint32_t> Parent_map;
vector<KrakenDB*> KrakenDatabases;
string Classified_output_file, Unclassified_output_file, Kraken_output_file, Report_output_file, TaxDB_file;
ostream *Classified_output;
ostream *Unclassified_output;
ostream *Kraken_output;
ostream *Report_output;
vector<ofstream*> Open_fstreams;
vector<ogzstream*> Open_gzstreams;
size_t Work_unit_size = DEF_WORK_UNIT_SIZE;
TaxonomyDB<uint32_t> taxdb;

uint64_t total_classified = 0;
uint64_t total_sequences = 0;
uint64_t total_bases = 0;

inline bool ends_with(std::string const & value, std::string const & ending)
{
        if (ending.size() > value.size()) return false;
            return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

ostream* cout_or_file(string file) {
    if (file == "-")
      return &cout;

    if (ends_with(file, ".gz")) {
      ogzstream* ogzs = new ogzstream(file.c_str());
      Open_gzstreams.push_back(ogzs);
      return ogzs;
    } else {
      ofstream* ofs = new ofstream(file.c_str());
      Open_fstreams.push_back(ofs);
      return ofs;
    }
}



void loadKrakenDB(KrakenDB& database, string DB_filename, string Index_filename) {
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

int main(int argc, char **argv) {
  #ifdef _OPENMP
  omp_set_num_threads(1);
  #endif

  parse_command_line(argc, argv);
  
  if (!TaxDB_file.empty()) {
	  taxdb = TaxonomyDB<uint32_t>(TaxDB_file);
      for (const auto & tax : taxdb.taxIDsAndEntries) {
          if (tax.first != 0)
          Parent_map[tax.first] = tax.second.parentTaxonomyID;
      }
      Parent_map[1] = 0;
  } else {
      cerr << "TaxDB argument is required!" << endl;
      return 1;
  }

  if (Populate_memory)
    cerr << "Loading database(s)... " << endl;

  // TODO: Check DB_filenames and Index_filesnames have the same length
  for (size_t i=0; i < DB_filenames.size(); ++i) {
    //cerr << "\t " << DB_filenames[i] << endl;
    static QuickFile db_file;
    db_file.open_file(DB_filenames[i]);
    if (Populate_memory)
      db_file.load_file();
    static KrakenDB Database = KrakenDB(db_file.ptr());
    KmerScanner::set_k(Database.get_k());
  
    static QuickFile idx_file;
    idx_file.open_file(Index_filenames[i]);
    if (Populate_memory)
      idx_file.load_file();
    static KrakenDBIndex db_index(idx_file.ptr());
    Database.set_index(&db_index);
    
  
    KrakenDatabases.push_back(&Database);
  }

  // TODO: Check all databases have the same k
  KmerScanner::set_k(KrakenDatabases[0]->get_k());

  if (Populate_memory)
    cerr << "\ncomplete." << endl;

  if (Print_classified) {
    Classified_output = cout_or_file(Classified_output_file);
  }

  if (Print_unclassified) {
    Unclassified_output = cout_or_file(Unclassified_output_file);
  }

  if (! Kraken_output_file.empty()) {
    if (Kraken_output_file == "off")
      Print_kraken = false;
    else {
      cerr << "Writing Kraken output to " << Kraken_output_file << endl;
      Kraken_output = cout_or_file(Kraken_output_file);
    }
  } else {
    Kraken_output = &cout;
  }

  if (!Report_output_file.empty()) {
     Print_kraken_report = true;
      cerr << "Writing Kraken report output to " << Report_output_file << endl;
     Report_output = cout_or_file(Report_output_file);
  }

  cerr << "Print_kraken: " << Print_kraken << "; Print_kraken_report: " << Print_kraken_report << "; k: " << uint32_t(KrakenDatabases[0]->get_k()) << endl;

  struct timeval tv1, tv2;
  gettimeofday(&tv1, NULL);
  for (int i = optind; i < argc; i++)
    process_file(argv[i]);
  gettimeofday(&tv2, NULL);

  std::cerr << "Finishing up ..\n";

  if (Print_kraken_report) {
	taxdb.fillCounts(taxon_counts);
	TaxReport<uint32_t> rep = TaxReport<uint32_t>(*Report_output, taxdb, false);
	rep.printReport("kraken","blu");
  }

  for (ofstream* ofs : Open_fstreams) {
    ofs->close();
  }
  for (ogzstream* ogzs : Open_gzstreams) {
    ogzs->close();
  }


  report_stats(tv1, tv2);

  return 0;
}

void report_stats(struct timeval time1, struct timeval time2) {
  time2.tv_usec -= time1.tv_usec;
  time2.tv_sec -= time1.tv_sec;
  if (time2.tv_usec < 0) {
    time2.tv_sec--;
    time2.tv_usec += 1000000;
  }
  double seconds = time2.tv_usec;
  seconds /= 1e6;
  seconds += time2.tv_sec;

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

void process_file(char *filename) {
  string file_str(filename);
  DNASequenceReader *reader;
  DNASequence dna;

  if (Fastq_input)
    reader = new FastqReader(file_str);
  else
    reader = new FastaReader(file_str);

  #pragma omp parallel
  {
    vector<DNASequence> work_unit;
    ostringstream kraken_output_ss, classified_output_ss, unclassified_output_ss;

    while (reader->is_valid()) {
      work_unit.clear();
      size_t total_nt = 0;
      #pragma omp critical(get_input)
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
      
      unordered_map<uint32_t, ReadCounts> my_taxon_counts;
      uint64_t my_total_classified = 0;
      kraken_output_ss.str("");
      classified_output_ss.str("");
      unclassified_output_ss.str("");
      for (size_t j = 0; j < work_unit.size(); j++)
        my_total_classified += 
            classify_sequence( work_unit[j], kraken_output_ss,
                           classified_output_ss, unclassified_output_ss,
                           my_taxon_counts);

      #pragma omp critical(write_output)
      {
        total_classified += my_total_classified;
        for (auto &it : my_taxon_counts) {
            taxon_counts[it.first] += it.second;
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
        if (Print_Progress && total_sequences % 100000 < work_unit.size()) 
          cerr << "\rProcessed " << total_sequences << " sequences (" << total_classified << " classified) ...";
      }
    }
  }  // end parallel section

  delete reader;
}

uint32_t get_taxon_for_kmer(KrakenDB& database, uint64_t* kmer_ptr, uint64_t& current_bin_key,
		int64_t& current_min_pos, int64_t& current_max_pos) {
	uint32_t* val_ptr = database.kmer_query(
			database.canonical_representation(*kmer_ptr), &current_bin_key,
			&current_min_pos, &current_max_pos);
	uint32_t taxon = val_ptr ? *val_ptr : 0;
	return taxon;
}

bool classify_sequence(DNASequence &dna, ostringstream &koss,
                       ostringstream &coss, ostringstream &uoss,
                       unordered_map<uint32_t, ReadCounts>& my_taxon_counts) {
  vector<uint32_t> taxa;
  vector<uint8_t> ambig_list;
  unordered_map<uint32_t, uint32_t> hit_counts;
  uint64_t *kmer_ptr;
  uint32_t taxon = 0;
  uint32_t hits = 0;  // only maintained if in quick mode

  uint64_t current_bin_key; int64_t current_min_pos = 1;  int64_t current_max_pos = 0;

  if (dna.seq.size() >= KrakenDatabases[0]->get_k()) {
    KmerScanner scanner(dna.seq);
    while ((kmer_ptr = scanner.next_kmer()) != NULL) {
      taxon = 0;
      if (scanner.ambig_kmer()) {
        ambig_list.push_back(1);
      }
      else {
        ambig_list.push_back(0);

        for (auto& db : KrakenDatabases) {
            taxon = get_taxon_for_kmer(*db, kmer_ptr, current_bin_key, current_min_pos, current_max_pos);
            if (taxon) break;
        }

        my_taxon_counts[taxon].add_kmer(*kmer_ptr);

        if (taxon) {
          hit_counts[taxon]++;
          if (Quick_mode && ++hits >= Minimum_hit_count)
            break;
        }
      }
      taxa.push_back(taxon);
    }
  }

  uint32_t call = 0;
  if (Quick_mode)
    call = hits >= Minimum_hit_count ? taxon : 0;
  else
    call = resolve_tree(hit_counts, Parent_map);

  ++(my_taxon_counts[call].n_reads);

  if (Print_unclassified || Print_classified) {
    ostringstream *oss_ptr = call ? &coss : &uoss;
    bool print = call ? Print_classified : Print_unclassified;
    if (print) {
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
  }

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
  koss << dna.id << "\t" << call << "\t" << dna.seq.size() << "\t";

  if (Quick_mode) {
    koss << "Q:" << hits;
  }
  else {
    if (taxa.empty())
      koss << "0:0";
    else
      koss << hitlist_string(taxa, ambig_list);
  }

  if (Print_sequence)
      koss << "\t" << dna.seq;

  koss << "\n";
  return call;
}

string hitlist_string(vector<uint32_t> &taxa, vector<uint8_t> &ambig)
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
  while ((opt = getopt(argc, argv, "d:i:t:u:n:m:o:qfcC:U:Ma:r:s")) != -1) {
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
  if (optind == argc) {
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
       << "  -t #             Number of threads" << endl
       << "  -u #             Thread work unit size (in bp)" << endl
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
