/*
 * Copyright 2017, Florian Breitwieser
 *
 * This file is part of the KrakenUniq taxonomic sequence classification system.
 *
 * KrakenUniq is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KrakenUniq is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kraken.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "kraken_headers.hpp"
#include "krakendb.hpp"
#include "quickfile.hpp"
#include <iostream>
#include <fstream>

using namespace std;
using namespace kraken;

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "USAGE:\n" 
      << "dump_db_kmers DATABASE\n"
      << "\n"
      << "Dumps database k-mers as 64-bit numbers \n";
    return 1;
  }

  char *db_name = argv[1];
  QuickFile db_file;
  db_file.open_file(db_name);
  //db_file.load_file();
  //cerr << "Fully loaded\n";
  KrakenDB db(db_file.ptr());

  char* ptr = db.get_ptr();
  //char* pair_ptr = db.get_pair_ptr();
  uint64_t key_len = db.get_key_len();     // how many bytes does each key occupy?
  //uint64_t val_len = db.get_val_len();     // how many bytes does each value occupy?
  uint64_t key_ct = db.get_key_ct();      // how many key/value pairs are there?
  uint64_t pair_sz = db.pair_size();       // how many bytes does each pair occupy?

  if (ptr == NULL) { 
    std::cerr << "Kraken database pointer is NULL [pair_sz: " << pair_sz << ", key_ct: "<<key_ct<<", key_len: "<< key_len<<"]!" << std::endl;
    exit(1);
  }
  for (uint64_t i = 0; i < key_ct; i++) {
    uint64_t* kmer = (uint64_t *) (ptr + pair_sz * i);
    cout << *kmer << '\n';
  }
}

