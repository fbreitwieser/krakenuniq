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

#include "hyperloglogplus.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "USAGE:\n" 
      << "count_unique PRECISION SPARSE TEST_MODE\n"
      << "\n"
      << "Valid precision values: 10-18. SPARSE can be 0 or 1. If TEST_MODE is 1, then a HLL estimate is given with each number. \n"
      << "Returns the cardinality of the input stream (has to be uint64_t)\n";
    return 1;
  }

  size_t p = stoi(argv[1]);
  bool sparse = bool(stoi(argv[2]));
  bool test_mode = bool(stoi(argv[3]));
  HyperLogLogPlusMinus<uint64_t> hll(p, sparse); // unique k-mer count per taxon
  uint64_t nr;
  uint64_t ctr = 0;
  if (test_mode) {
    cout << "observed\testimated\n";
  }
  while (cin >> nr) {
    hll.add(nr);
    if (test_mode) {
      cout << ++ctr << '\t' << hll.cardinality() << '\n';
    }
  }
  if (!test_mode) {
    cout << hll.cardinality() << endl;
  }
  
}

