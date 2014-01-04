/*
 * Copyright 2013-2014, Derrick Wood <dwood@cs.umd.edu>
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
#ifdef _OPENMP
omp_lock_t *Locks;
#endif

static void parse_command_line(int argc, char **argv);
static void bin_and_sort_data(KrakenDB &in, KrakenDB &out);
static void usage(int exit_code=EX_USAGE);

int main(int argc, char **argv) {
  #ifdef _OPENMP
  omp_set_num_threads(1);
  #endif

  parse_command_line(argc, argv);

  QuickFile input_db_file(Input_DB_filename);
  KrakenDB input_db(input_db_file.ptr());

  input_db.make_index(Index_filename, Bin_key_nt);
  QuickFile index_file(Index_filename);
  KrakenDBIndex db_index(index_file.ptr());
  input_db.set_index(&db_index);

  // Copy header from input DB to output DB, reopen output DB in R/W
  // open/reopen necessary due to "feature" of QuickFile/KrakenDB code
  //   (could be skipped if we didn't want to use KrakenDB features in basd()
  QuickFile output_file;
  char *temp_file;
  if (Operate_in_RAM) {
    temp_file = new char[ input_db_file.size() ];
  }
  else {
    output_file.open_file(Output_DB_filename, "w", input_db_file.size());
    temp_file = output_file.ptr();
  }
  memcpy(temp_file, input_db_file.ptr(), input_db.header_size());
  if (! Operate_in_RAM) {
    output_file.close_file();
    output_file.open_file(Output_DB_filename, "rw");
    temp_file = output_file.ptr();
  }
  KrakenDB output_db(temp_file);

  bin_and_sort_data(input_db, output_db);
  
  if (Operate_in_RAM) {
    size_t filesize = input_db_file.size();
    input_db_file.close_file();
    output_file.open_file(Output_DB_filename, "w", filesize);
    memcpy(output_file.ptr(), temp_file, filesize);
  }
  return 0;
}

static void bin_and_sort_data(KrakenDB &in, KrakenDB &out) {
  char *in_ptr = in.get_pair_ptr();
  char *out_ptr = out.get_pair_ptr();
  uint8_t nt = in.get_index()->indexed_nt();
  uint64_t *offsets = in.get_index()->get_array();
  uint64_t entries = 1ull << (nt * 2);

  #ifdef _OPENMP
  Locks = new omp_lock_t[entries];
  for (uint64_t i = 0; i < entries; i++)
    omp_init_lock(&Locks[i]);
  #endif

  // Create a copy of the offsets array for use as insertion positions
  vector<uint64_t> pos(offsets, offsets + entries);
  #pragma omp parallel for schedule(dynamic,400)
  for (uint64_t i = 0; i < in.get_key_ct(); i++) {
    uint64_t bin_key = in.bin_key(* (uint64_t *)(in_ptr + i * in.pair_size()));
    #ifdef _OPENMP
    omp_set_lock(&Locks[bin_key]);
    #endif
    #pragma omp flush(pos)
    char *pair_pos = out_ptr + in.pair_size() * pos[bin_key];
    // Copy pair into correct bin (but not final position)
    memcpy(pair_pos, in_ptr + i * in.pair_size(), in.pair_size());
    if (Zero_vals)
      memset(pair_pos + in.get_key_len(), 0, in.get_val_len());
    #pragma omp flush(pos)
    pos[bin_key]++;
    #ifdef _OPENMP
    omp_unset_lock(&Locks[bin_key]);
    #endif
  }

  #ifdef _OPENMP
  for (uint64_t i = 0; i < entries; i++)
    omp_destroy_lock(&Locks[i]);
  delete Locks;
  #endif

  // Sort all bins
  #pragma omp parallel for schedule(dynamic)
  for (uint64_t i = 0; i < entries; i++) {
    qsort(out_ptr + offsets[i] * out.pair_size(),
          offsets[i+1] - offsets[i], out.pair_size(),
          KrakenDB::pair_cmp);
  }
}

void parse_command_line(int argc, char **argv) {
  int opt;
  int sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "n:d:o:i:t:zM")) != -1) {
    switch (opt) {
      case 'n' :
        sig = atoi(optarg);
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
        sig = atoi(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "can't use nonpositive thread count");
        #ifdef _OPENMP
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
