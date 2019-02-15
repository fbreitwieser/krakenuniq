/*
 * Copyright 2018, Florian P Breitwieser
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
#include "krakenutil.hpp"
#include "seqreader.hpp"
#include "hyperloglogplus.hpp"


#define SKIP_LEN 50000

using namespace std;
using namespace kraken;

void parse_command_line(int argc, char **argv);
void usage(int exit_code=EX_USAGE);
uint64_t count_unique(int, bool);

int Num_threads = 1;
int k = 31;
int hll_precision = 14;

int main(int argc, char **argv) {
  #ifdef _OPENMP
  omp_set_num_threads(1);
  #endif

  parse_command_line(argc, argv);

  KmerScanner::set_k(k);
  cout << count_unique(hll_precision, true) << endl;

  return 0;
}

uint64_t count_unique(int p, bool sparse) {
  FastaReader reader("/dev/fd/0");
  DNASequence dna;
  HyperLogLogPlusMinus<uint64_t> counter(p, false);
  
  while (true) {
    dna = reader.next_sequence();
    if (! reader.is_valid())
      break;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (size_t i = 0; i < dna.seq.size(); i += SKIP_LEN) {
      HyperLogLogPlusMinus<uint64_t> mycounter(p, false);
      KmerScanner scanner(dna.seq, i, i + SKIP_LEN + k - 1);
      uint64_t *kmer_ptr;

      while ((kmer_ptr = scanner.next_kmer()) != NULL) {
        if (scanner.ambig_kmer())
          continue;
        mycounter.insert(*kmer_ptr);
      }
#ifdef _OPENMP
      #pragma omp critical(set_access)
#endif
      counter += mycounter;
    }
  }
  return (uint64_t) counter.cardinality();
}

void parse_command_line(int argc, char **argv) {
  int opt;
  long long sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "t:k:m:p:")) != -1) {
    switch (opt) {
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
      case 'k' :
        sig = atoll(optarg);
        if (sig <= 0)
          errx(EX_USAGE, "k can't be <= 0");
        if (sig > 31)
          errx(EX_USAGE, "k can't be > 31");
        k = sig;
        break;
      case 'p' :
        hll_precision = atoi(optarg);
        break;
      default:
        usage();
        break;
    }
  }

  if (k == 0)
    usage(EX_USAGE);
}

void usage(int exit_code) {
  cerr << "Usage: count-unique [options]" << endl
       << "  Estimate the number of k-mer in stdin using HLL." << endl
       << "Options: " << endl
       << "  -k #          Length of k-mers ["<<k<<"]" << endl
       << "  -t #          Number of threads ["<<Num_threads<<"]" << endl
       << "  -p INT        HLL precision ["<<hll_precision<<"]" << endl
       << "  -h            Print this message" << endl;
  exit(exit_code);
}
