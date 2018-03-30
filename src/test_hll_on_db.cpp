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
 * along with KrakenHLL.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "hyperloglogplus.hpp"
#include "kraken_headers.hpp"
#include "quickfile.hpp"
#include "krakendb.hpp"
#include <iostream>
#include <fstream>
#include <random>

using namespace std;
using namespace kraken;

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "USAGE: test_hll_on_db DATABASE SPARSE NR_KMERS\n"
       "\n"
       "Tests precision values 10-18 on the k-mers in the database. SPARSE can be 0 or 1. \n"
       " Each k-mer in the database is sampled with a probability of NR_KMERS / (# of k-mers in the database).\n";
    return 1;
  }


  char *db_name = argv[1];
  QuickFile db_file;
  db_file.open_file(db_name);
  db_file.load_file();
  cerr << "Fully loaded\n";
  KrakenDB db(db_file.ptr());

  //size_t p = stoi(argv[2]);
  bool sparse = bool(stoi(argv[2]));
  size_t nr = stoi(argv[3]);

  HyperLogLogPlusMinus<uint64_t> hll10(10, sparse); // unique k-mer count per taxon
  HyperLogLogPlusMinus<uint64_t> hll11(11, sparse); // unique k-mer count per taxon
  HyperLogLogPlusMinus<uint64_t> hll12(12, sparse); // unique k-mer count per taxon
  HyperLogLogPlusMinus<uint64_t> hll13(13, sparse); // unique k-mer count per taxon
  HyperLogLogPlusMinus<uint64_t> hll14(14, sparse); // unique k-mer count per taxon
  HyperLogLogPlusMinus<uint64_t> hll15(15, sparse); // unique k-mer count per taxon
  HyperLogLogPlusMinus<uint64_t> hll16(16, sparse); // unique k-mer count per taxon
  HyperLogLogPlusMinus<uint64_t> hll17(17, sparse); // unique k-mer count per taxon
  HyperLogLogPlusMinus<uint64_t> hll18(18, sparse); // unique k-mer count per taxon

  hll10.use_n_observed = false;
  hll11.use_n_observed = false;
  hll12.use_n_observed = false;
  hll13.use_n_observed = false;
  hll14.use_n_observed = false;
  hll15.use_n_observed = false;
  hll16.use_n_observed = false;
  hll17.use_n_observed = false;
  hll18.use_n_observed = false;

  char* ptr = db.get_ptr();
  //char* pair_ptr = db.get_pair_ptr();
  uint64_t key_len = db.get_key_len();     // how many bytes does each key occupy?
  //uint64_t val_len = db.get_val_len();     // how many bytes does each value occupy?
  uint64_t key_ct = db.get_key_ct();      // how many key/value pairs are there?
  uint64_t pair_sz = db.pair_size();       // how many bytes does each pair occupy?

  if (nr > key_ct) {
    cerr << nr << " is greater than " << key_ct << "!!!" << endl;
    exit(1);
  }

  if (ptr == NULL) { 
    std::cerr << "Kraken database pointer is NULL [pair_sz: " << pair_sz << ", key_ct: "<<key_ct<<", key_len: "<< key_len<<"]!" << std::endl;
    exit(1);
  }
  double prob = double(nr)/double(key_ct);
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0.0, 1.0);

  //double lg_range = log10(nr)-2.;

  cout << "precision\ttrue_count\testimate\tertl_estimate\n";
  size_t ctr = 0;
  for (uint64_t i = 0; i < key_ct; i++) {
    if (dis(gen) < prob) {
      uint64_t* kmer = (uint64_t *) (ptr + pair_sz * i);
      //uint32_t* taxon = (uint32_t *) (ptr + pair_sz * i + key_len);
      //if (taxon == NULL) {
      //  std::cerr << "taxon is NULL (i is " << i << " and key_ct is " << key_ct << ")" << std::endl;
      hll10.add(*kmer);
      hll11.add(*kmer);
      hll12.add(*kmer);
      hll13.add(*kmer);
      hll14.add(*kmer);
      hll15.add(*kmer);
      hll16.add(*kmer);
      hll17.add(*kmer);
      hll18.add(*kmer);
      ++ctr;

	  double log_ctr = log10(ctr);

      uint64_t last_dec = pow(10,floor(log10(ctr)));
	  bool do_print = floor(log10(ctr)) == log10(ctr) || (ctr > 100 && (ctr*10) % last_dec == 0);
	  //bool do_print = dis(gen)
      if (do_print) {
        cout << 10 << '\t' << ctr << '\t' << hll10.cardinality() << '\t' << hll10.ertlCardinality() << '\n';
        cout << 11 << '\t' << ctr << '\t' << hll11.cardinality() << '\t' << hll11.ertlCardinality() << '\n';
        cout << 12 << '\t' << ctr << '\t' << hll12.cardinality() << '\t' << hll12.ertlCardinality() << '\n';
        cout << 13 << '\t' << ctr << '\t' << hll13.cardinality() << '\t' << hll13.ertlCardinality() << '\n';
        cout << 14 << '\t' << ctr << '\t' << hll14.cardinality() << '\t' << hll14.ertlCardinality() << '\n';
        cout << 15 << '\t' << ctr << '\t' << hll15.cardinality() << '\t' << hll15.ertlCardinality() << '\n';
        cout << 16 << '\t' << ctr << '\t' << hll16.cardinality() << '\t' << hll16.ertlCardinality() << '\n';
        cout << 17 << '\t' << ctr << '\t' << hll17.cardinality() << '\t' << hll17.ertlCardinality() << '\n';
        cout << 18 << '\t' << ctr << '\t' << hll18.cardinality() << '\t' << hll18.ertlCardinality() << '\n';
      }
    }
  }
}

