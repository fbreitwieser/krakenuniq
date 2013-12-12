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

using namespace std;
using namespace kraken;

map<uint64_t, string> requests;
uint64_t request_count = 0;

void fill_request_map(char *filename);
void report_taxo_numbers(char *filename);

int main(int argc, char **argv) {
  if (argc < 3) {
    cerr << "Usage: make_file_to_taxon_map <gi to taxid map> <gi to file list>"
         << endl;
    return 1;
  }
  char *map_filename = argv[1];
  char *list_filename = argv[2];
  fill_request_map(list_filename);
  report_taxo_numbers(map_filename);

  return 0;
}

void report_taxo_numbers(char *filename) {
  string file_str = filename;
  QuickFile file(file_str);
  char *fptr, *fptr_start;
  fptr_start = fptr = file.ptr();
  size_t file_size = file.size();

  while (request_count > 0 && (size_t)(fptr - fptr_start) < file_size) {
    char *nl_ptr = strchr(fptr, '\n');
    uint64_t gi = atoll(fptr);
    
    if (requests.count(gi) > 0) {
      char *tab_ptr = strchr(fptr, '\t');
      cout << requests[gi] << "\t";
      cout.write(tab_ptr + 1, nl_ptr - tab_ptr);
      request_count--;
    }

    fptr = nl_ptr + 1;
  }

  file.close_file();
}

void fill_request_map(char *filename) {
  string file_str = filename;
  QuickFile file(file_str);
  char *fptr, *fptr_start;
  fptr_start = fptr = file.ptr();
  size_t file_size = file.size();

  while ((size_t)(fptr - fptr_start) < file_size) {
    char *nl_ptr = strchr(fptr, '\n');
    char *sep_ptr = strchr(fptr, '|');
    uint64_t gi = atoll(fptr);
    requests[gi] = string(sep_ptr + 1, nl_ptr - sep_ptr - 1);
    request_count++;
    fptr = nl_ptr + 1;
  }

  file.close_file();
}
