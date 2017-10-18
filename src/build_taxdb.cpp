/*
 * Copyright 2017, Florian Breitwieser
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

#include "taxdb.h"
#include "quickfile.hpp"
#include <iostream>
#include <fstream>
#include <unordered_map>

using namespace std;

int main(int argc, char **argv) {
	if (argc < 2 || argc > 4) {
      std::cerr << "USAGE:\n" 
        << "With two or three arguments, echo taxDB based on NCBI taxonomy dump:\n"
        << "build_taxdb names.dmp nodes.dmp [taxon-counts]\n"
        << "\n"
        << "With one argument, read in taxDB and echo it again for consistency checks:\n"
        << "build_taxdb taxDB\n";
      return 1;
    }
        
    TaxonomyDB<uint32_t, uint32_t> taxdb;
    if (argc == 2) {
    taxdb = TaxonomyDB<uint32_t, uint32_t> ((string)argv[1]);
    } else {
    taxdb = TaxonomyDB<uint32_t, uint32_t> ((string)argv[1], (string)argv[2]);
    }
    if (argc == 4) {
        ifstream ifs(argv[3]);
        uint32_t taxon; uint64_t count;
        while (ifs >> taxon >> count) {
            taxdb.setGenomeSize(taxon, count);
        }
        taxdb.genomeSizes_are_set = true;
    }
    taxdb.writeTaxonomyIndex(std::cout);
}
