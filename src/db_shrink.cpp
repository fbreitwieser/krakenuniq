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

string Input_DB_filename, Output_DB_filename;
uint64_t Output_count = 0;
size_t Offset = 1;
bool Operate_in_RAM = false;

static void parse_command_line(int argc, char **argv);
static void usage(int exit_code=EX_USAGE);

int main(int argc, char **argv) {
  parse_command_line(argc, argv);
  
  uint64_t key_bits, val_len, key_count, key_len;
  uint64_t pair_size;

  // NOTE: Skipping the KrakenDB use here does make the code a
  //   bit redundant and messy, but it keeps the memory footprint
  //   low, which is important for this particular program.

  // Check file magic and get metadata to calculate header size
  char *buffer = new char[8];
  ifstream input_file(Input_DB_filename.c_str(), std::ifstream::binary);
  input_file.read(buffer, 8);
  if (strncmp("JFLISTDN", buffer, 8) != 0) {
    errx(EX_DATAERR, "input file not Jellyfish v1 database");
  }
  input_file.read(buffer, 8);
  memcpy(&key_bits, buffer, 8);
  size_t header_size = 72 + 2 * (4 + 8 * key_bits);
  delete buffer;

  // Read in header and get remaining metadata
  buffer = new char[header_size];
  input_file.seekg(0, ios_base::beg);
  input_file.read(buffer, header_size);
  memcpy(&val_len, buffer + 16, 8);
  memcpy(&key_count, buffer + 48, 8);
  key_len = key_bits / 8 + !! (key_bits % 8);
  pair_size = key_len + val_len;

  if (Output_count > key_count) {
    errx(EX_DATAERR, "Requested new key count %llu larger than old key count %llu, aborting...",
                      (long long unsigned int) Output_count,
                      (long long unsigned int) key_count);
  }

  // Change key count
  memcpy(buffer + 48, &Output_count, 8);

  // Copy input header to output file
  ofstream output_file(Output_DB_filename.c_str(), std::ofstream::binary);
  output_file.write(buffer, header_size);

  delete buffer;

  // Prep buffer for scan/select loop
  // We select one pair (the last) per "block"
  size_t block_size = key_count / Output_count;
  if (block_size < Offset) {
    errx(EX_DATAERR, "offset %llu larger than block size %llu, aborting.",
         (long long unsigned int) Offset,
         (long long unsigned int) block_size);
  }
  // Some blocks have an extra element
  size_t odd_block_count = key_count % Output_count;
  buffer = new char[pair_size * (block_size + 1)];
  size_t remaining_input_pairs = key_count;
  size_t current_output_count = 0;

  while (remaining_input_pairs) {
    size_t pairs_to_read;
    if (odd_block_count == 0) {
      pairs_to_read = block_size;
    }
    else {
      pairs_to_read = block_size + 1;
      odd_block_count--;
    }

    input_file.read(buffer, pairs_to_read * pair_size);
    remaining_input_pairs -= pairs_to_read;

    output_file.write(buffer + pair_size * (pairs_to_read - Offset), pair_size);
    current_output_count++;
    if (current_output_count % 10000 == 0) {
      cerr << "\rWritten " << current_output_count << "/" << Output_count << " k-mers to new file";
    }
  }
  cerr << "\rWrote " << current_output_count << "/" << Output_count << " k-mers to new file  \n";

  input_file.close();
  output_file.close();

  return 0;
}

void parse_command_line(int argc, char **argv) {
  int opt;
  long long sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "d:o:n:O:")) != -1) {
    switch (opt) {
      case 'n' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "output count must be positive");
        Output_count = sig;
        break;
      case 'd' :
        Input_DB_filename = optarg;
        break;
      case 'o' :
        Output_DB_filename = optarg;
        break;
      case 'O' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "offset count cannot be negative");
        Offset = sig;
        break;
      default:
        usage();
        break;
    }
  }

  if (Input_DB_filename.empty() || Output_DB_filename.empty() || ! Output_count)
    usage();
}

void usage(int exit_code) {
  cerr << "Usage: db_shrink [-O offset] <-d input db> <-o output db> <-n output count>\n";
  exit(exit_code);
}
