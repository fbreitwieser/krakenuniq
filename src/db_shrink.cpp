/*
 * Copyright 2013, Derrick Wood <dwood@cs.umd.edu>
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

string Input_DB_filename, Output_DB_filename;
int Shrink_percentage = 0;
bool Operate_in_RAM = false;

static void parse_command_line(int argc, char **argv);
static void usage(int exit_code=EX_USAGE);

int main(int argc, char **argv) {
  #ifdef _OPENMP
  omp_set_num_threads(1);
  #endif

  parse_command_line(argc, argv);

  QuickFile input_db_file(Input_DB_filename);
  KrakenDB input_db(input_db_file.ptr());

  size_t new_file_size = input_db.header_size();
  uint64_t key_ct = input_db.get_key_ct();
  uint64_t new_key_ct = key_ct * Shrink_percentage / 100;
  uint64_t pair_size = input_db.pair_size();
  new_file_size += new_key_ct * pair_size;

  QuickFile output_file;
  char *temp_file = NULL;
  char *outptr;
  if (Operate_in_RAM) {
    temp_file = new char[new_file_size];
  }
  else {
    output_file.open_file(Output_DB_filename, "rw", new_file_size);
    temp_file = output_file.ptr();
  }
  outptr = temp_file;
  // Copy input header to output file
  memcpy(outptr, input_db.get_ptr(), input_db.header_size());
  // Change key count
  memcpy(outptr + 48, &new_key_ct, 8);
  outptr += input_db.header_size();

  char *inptr = input_db.get_pair_ptr();
  // For each 100 consecutive pairs, select the last P pairs of the block
  for (uint64_t i = 0; i < new_key_ct; i += Shrink_percentage) {
    inptr += (100 - Shrink_percentage) * pair_size;
    outptr += Shrink_percentage * pair_size;
    memcpy(outptr, inptr, Shrink_percentage * pair_size);
    inptr += Shrink_percentage * pair_size;
  }

  input_db_file.close_file();

  if (Operate_in_RAM) {
    // Actual file hasn't been opened yet if -M used
    output_file.open_file(Output_DB_filename, "w", new_file_size);
    memcpy(output_file.ptr(), temp_file, new_file_size);
    delete temp_file;
  }

  output_file.close_file();

  return 0;
}

void parse_command_line(int argc, char **argv) {
  int opt;
  int sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "d:o:p:M")) != -1) {
    switch (opt) {
      case 'p' :
        sig = atoi(optarg);
        if (sig < 1 || sig >= 100)
          errx(EX_USAGE, "shrink percentage out of range");
        Shrink_percentage = sig;
        break;
      case 'd' :
        Input_DB_filename = optarg;
        break;
      case 'o' :
        Output_DB_filename = optarg;
        break;
      case 'M' :
        Operate_in_RAM = true;
        break;
      default:
        usage();
        break;
    }
  }

  if (Input_DB_filename.empty() || Output_DB_filename.empty() || ! Shrink_percentage)
    usage();
}

void usage(int exit_code) {
  cerr << "Usage: db_shrink [-M] <-d input db> <-o output db> <-p percentage>\n";
  exit(exit_code);
}
