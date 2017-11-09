/*
 * Copyright 2017, Florian Breitwieser
 *
 * This file is part of the KrakenHLL taxonomic sequence classification system.
 *
 * KrakenHLL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KrakenHLL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kraken.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "taxdb.h"
#include "quickfile.hpp"
#include <iostream>
#include <fstream>
#include <unordered_map>

using namespace std;

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "Usage: dump_taxdb taxDB names.dmp nodes.dmp\n";
    return 1;
  }

	cerr << "Reading taxonomy database from " << argv[1] << ", writing nodes dump to " << argv[3] << " and names dump to " << argv[2] << "." << endl;
  TaxonomyDB<uint32_t> taxdb {(string)argv[1]};
  ofstream names_file(argv[2]);
  names_file.exceptions(ifstream::failbit | ifstream::badbit);
  ofstream nodes_file(argv[3]);
  nodes_file.exceptions(ifstream::failbit | ifstream::badbit);

  for (auto it = taxdb.entries.begin(); it != taxdb.entries.end(); ++it) {
    const auto &taxon = *it;
    std::string scientificName;
    uint32_t parentTaxonomyID = taxon.second.parent == NULL? taxon.first : taxon.second.parent->taxonomyID;
    nodes_file << taxon.second.taxonomyID 
      << "\t|\t" << parentTaxonomyID
      << "\t|\t" << taxon.second.rank
      << endl; // there are further columns, but Kraken does not care about them
    
    names_file << taxon.second.taxonomyID 
      << "\t|\t" << taxon.second.scientificName
      << "\t|\t" 
      << "\t|\t" << "scientific name" << endl;
  }
  names_file.close();
  nodes_file.close();
}
