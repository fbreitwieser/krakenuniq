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

using namespace std;
using namespace kraken;

string Input_DB_filename, Output_DB_filename, Index_filename;
uint8_t Bin_key_nt = 15;
int Num_threads = 1;
bool Zero_vals = false;
bool Operate_in_RAM = false;
// Global until I can find a way to pass this to the sorting function
size_t Key_len = 8;

static int pair_cmp(const void *a, const void *b);
static void parse_command_line(int argc, char **argv);
static void bin_and_sort_data(KrakenDB &kdb, char *data, KrakenDBIndex &idx);
static void usage(int exit_code=EX_USAGE);

int main(int argc, char **argv) {
  #ifdef _OPENMP
  omp_set_num_threads(1);
  #endif

  parse_command_line(argc, argv);

  QuickFile input_db_file(Input_DB_filename);
  KrakenDB *input_db = new KrakenDB(input_db_file.ptr());
  Key_len = input_db->get_key_len();
  uint64_t val_len = input_db->get_val_len();
  uint64_t key_ct = input_db->get_key_ct();
  input_db->make_index(Index_filename, Bin_key_nt);
  QuickFile index_file(Index_filename);
  KrakenDBIndex db_index(index_file.ptr());

  uint64_t skip_len = input_db->header_size();
  char *header = new char[ skip_len ];
  memcpy(header, input_db_file.ptr(), skip_len);

  delete input_db;
  // No longer legit for traversal, but allows access to KDB functions
  input_db = new KrakenDB(header);
  input_db_file.close_file();  // Stop using memory-mapped file

  char *data = new char[ key_ct * (Key_len + val_len) ];
  // Populate data w/ pairs from DB and sort bins in parallel
  bin_and_sort_data(*input_db, data, db_index);

  ofstream output_file(Output_DB_filename.c_str(), std::ofstream::binary);
  output_file.write(header, skip_len);
  output_file.write(data, key_ct * (Key_len + val_len));
  output_file.close();
  
  return 0;
}

static void bin_and_sort_data(KrakenDB &kdb, char *data, KrakenDBIndex &idx) {
  uint8_t nt = idx.indexed_nt();
  uint64_t *offsets = idx.get_array();
  uint64_t entries = 1ull << (nt * 2);
  uint64_t key_len = kdb.get_key_len();
  uint64_t val_len = kdb.get_val_len();
  uint64_t pair_size = key_len + val_len;
  char pair[pair_size];

  ifstream input_file(Input_DB_filename.c_str(), std::ifstream::binary);
  input_file.seekg(kdb.header_size(), ios_base::beg);

  // Create a copy of the offsets array for use as insertion positions
  vector<uint64_t> pos(offsets, offsets + entries);
  for (uint64_t i = 0; i < kdb.get_key_ct(); i++) {
    input_file.read(pair, pair_size);
    uint64_t kmer = 0;
    memcpy(&kmer, pair, key_len);
    uint64_t b_key = kdb.bin_key(kmer, nt);
    char *pair_pos = data + pair_size * pos[b_key]++;
    // Copy pair into correct bin (but not final position)
    memcpy(pair_pos, pair, pair_size);
    if (Zero_vals)
      memset(pair_pos + key_len, 0, val_len);
  }
  input_file.close();

  // Sort all bins
  #pragma omp parallel for schedule(dynamic)
  for (uint64_t i = 0; i < entries; i++) {
    qsort(data + offsets[i] * pair_size,
          offsets[i+1] - offsets[i], pair_size,
          pair_cmp);
  }
}

static int pair_cmp(const void *a, const void *b) {
  uint64_t aval = 0, bval = 0;
  memcpy(&aval, a, Key_len);
  memcpy(&bval, b, Key_len);
  if (aval < bval)
    return -1;
  else if (aval == bval)
    return 0;
  else
    return 1;
}

void parse_command_line(int argc, char **argv) {
  int opt;
  long long sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "n:d:o:i:t:zM")) != -1) {
    switch (opt) {
      case 'n' :
        sig = atoll(optarg);
        if (sig < 1 || sig > 31)
          errx(EX_USAGE, "bin key length out of range");
        Bin_key_nt = (uint8_t) sig;
        break;
      case 'd' :
        Input_DB_filename = optarg;
        break;
      case 'o' :
        Output_DB_filename = optarg;
        break;
      case 'i' :
        Index_filename = optarg;
        break;
      case 'M' :
        Operate_in_RAM = true;
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
      case 'z' :
        Zero_vals = true;
        break;
      default:
        usage();
        break;
    }
  }

  if (Input_DB_filename.empty() || Output_DB_filename.empty()
      || Index_filename.empty())
    usage();
}

void usage(int exit_code) {
  cerr << "Usage: db_sort [-z] [-M] [-t threads] [-n nt] <-d input db> <-o output db> <-i output idx>\n";
  exit(exit_code);
}
