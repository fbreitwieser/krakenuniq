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

#include "hyperloglogplus.hpp"
#include "khset64.h"
#include <iostream>
#include <fstream>
#include <random>
#include <limits>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

using namespace std;

int usage(int exit_code) {
  std::cerr << 
"count_unique: Get cardinality of 64-bit input stream\n"
"\n"
"Usage: count_unique OPTIONS\n"
"\n"
"OPTIONS:\n"
"  -p PRECISION   Precision in range of 10 to 18 (required)\n"
"  -r INT         Create INT random numbers, instead of counting from STDIN\n"
"  -s             Use sparse representation for smaller cardinalities\n"
"  -t             Test mode - print cardinalities regularily\n"
"  -y             Show relative error along with cardinality estimates\n"
"  -e             Use improved cardinality estimator by Otmar Ertl, too\n"
"  -U             Use exact cardinality counting implemented w/ unordered_set (not working w/ test mode)\n"
"  -K             Use exact cardinality counting implemented w/ khash (not working w/ test mode)\n";
    return exit_code;
  }

double rel_error(uint64_t est, uint64_t truth) {
  if (est >= truth) 
    return double(est - truth)/double(truth);
  else
    return -double(truth - est)/double(truth);
}

void print_card(HyperLogLogPlusMinus<uint64_t>& hll, uint64_t ctr, bool show_rel_error, bool heule_too, bool flajolet_too, bool ertl_too) { 
 
    uint64_t esth, este, estf;

    if (heule_too) {
      esth = hll.heuleCardinality();
      cout << ctr << '\t' << esth;
    }
    if (flajolet_too) {
      estf = hll.flajoletCardinality(false);
      cout << '\t' << estf;
    }
    if (ertl_too) {
      este = hll.ertlCardinality();
      cout << '\t' << este;
    }
  
  if (show_rel_error) {
    double esth_err, este_err;
    if (heule_too) {
      esth_err = rel_error(esth, ctr);
      cout << '\t' << esth_err;
    }
    if (flajolet_too) {
      cout << '\t' << rel_error(estf, ctr);
    }
    if (ertl_too) {
      este_err = rel_error(este, ctr);
      cout << '\t' << este_err;
    }
    if (ertl_too && heule_too) {
      if (abs(este_err) == abs(esth_err)) {
        cout << "\tequal";
      } else if (abs(este_err) < abs(esth_err)) {
        cout << "\tErtl won!";
      } else {
        cout << "\tHeule won!";
      }
    }
  }
  cout << '\n';
}


void add_to_hll(HyperLogLogPlusMinus<uint64_t>& hll, uint64_t nr, uint64_t& ctr, bool test_mode, bool show_rel_error, bool heule_too, bool flajolet_too, bool ertl_too) {
  hll.insert(nr);
  ++ctr;
  if (test_mode) {
    // prints numbers at powers of 10 and 10 points in-between
    uint64_t last_dec = pow(10,floor(log10(ctr)));
    if (floor(log10(ctr)) == log10(ctr) || 
        (ctr > 100 && (ctr*10) % last_dec == 0)) {
      print_card(hll, ctr, show_rel_error, heule_too, flajolet_too, ertl_too);
    }
  }
}


  template <typename T>
  unordered_set<T>& operator+=(unordered_set<T>& left, const unordered_set<T>& right) {
    left.insert(right.begin(), right.end());
    return left;
  }


int main(int argc, char **argv) {

  size_t p = 10;
  bool sparse = false;
  bool test_mode = false;
  bool heule_too = true;
  bool ertl_too = false;
  bool exact_counting_unordered_set = false;
  bool exact_counting_khash = false;
  bool flajolet_too = false;
  bool show_rel_error = false;
  bool use_stdin = true;
  size_t n_rand = 1;
  size_t n_redo = 2;

  int c;

  while ((c = getopt (argc, argv, "shtep:r:yx:fUK")) != -1)
    switch (c) {
      case 's': sparse = true; break;
      case 't': test_mode = true; break;
      case 'e': ertl_too = true; break;
      case 'U': exact_counting_unordered_set = true; break;
      case 'K': exact_counting_khash = true; break;
      case 'f': flajolet_too = true; break;
      case 'y': show_rel_error = true; break;
      case 'p': p = stoi(optarg); break;
      case 'r': use_stdin = false; 
                n_rand = stoll(optarg); 
                break;
      case 'x': n_redo = stoi(optarg); break;
      case 'h': return usage(0); break;
      case '?':
        if (optopt == 'p' || optopt == 'r')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }

  bool exact_counting = exact_counting_khash || exact_counting_unordered_set;
  HyperLogLogPlusMinus<uint64_t> hll(p, sparse); // unique k-mer count per taxon
  //HyperLogLogPlusMinus<uint64_t> hll(p, sparse, wang_mixer); // unique k-mer count per taxon

  if (test_mode && ! exact_counting) {
    cout << "observed\testimate_heule";
    if (flajolet_too) {
      cout << "\testimate_flajolet";
    }
    if (ertl_too)
      cout << "\testimate_ertl";
    if (show_rel_error) {
      cout << "\trel_error_heule";
      if (flajolet_too)
        cout << "\trel_error_flajolet";

      if (ertl_too)
        cout << "\trel_error_ertl";

      cout << "\twho_won";

    }
    cout << '\n';
  } 
  uint64_t ctr = 0;
  
  unordered_set<uint64_t> exact_counter;
  khset64_t exact_counter_khash;

  if (use_stdin) {
    uint64_t nr;
    while (cin >> nr) {
      if (exact_counting_unordered_set) {
	      exact_counter.insert(nr);
      } else if (exact_counter_khash) {
        exact_counter_khash.insert(nr);
      } else {
        add_to_hll(hll, nr, ctr, test_mode, show_rel_error, heule_too, flajolet_too, ertl_too);
      }
    }
    if (!test_mode) {
      if (exact_counting_unordered_set) {
        cout << exact_counter.size() << "\n";
      } else if (exact_counter_khash) {
        cout << exact_counter_khash.size() << "\n";
      } else {
        print_card(hll, ctr, show_rel_error, heule_too, flajolet_too, ertl_too);
      }
    }
  } else {
    // get random seed from random_device RNG
    std::random_device rd;
    // use 64-bit Mersenne Twister 19937 as RNG, seed with rd()
    std::mt19937_64 rng(rd()); 
    
    // Define output distribution (default range for unsigned: 0 to MAX)
    std::uniform_int_distribution<uint64_t> distr;
    for (size_t j = 0; j < n_redo; ++j) {
      unordered_set<uint64_t> exact_counter1;
      khset64_t exact_counter_khash1;
      for(size_t i = 0; i < n_rand; i++) {
        if (exact_counting_unordered_set) {
	        exact_counter1.insert(distr(rng));
        } else if (exact_counting_khash) {
          exact_counter_khash1.insert(distr(rng));
        } else {
          add_to_hll(hll, distr(rng), ctr, test_mode, show_rel_error, heule_too, flajolet_too, ertl_too);
        }
      }
      if (!test_mode) {
        if (exact_counting_unordered_set) {
          cout << exact_counter1.size() << "\n";
        } else if (exact_counting_khash) {
          cout << exact_counter_khash1.size() << "\n";
        } else {
          print_card(hll, ctr, show_rel_error, heule_too, flajolet_too, ertl_too);
        }
      }
	  //exact_counter += exact_counter1;
	  exact_counter_khash += std::move(exact_counter_khash1);
      cout << exact_counter_khash.size() << "\n";
      //cout << exact_counter.size() << "\n";
      hll.reset();
      exact_counter1.clear();
      ctr = 0;
    }
  }
  
}

