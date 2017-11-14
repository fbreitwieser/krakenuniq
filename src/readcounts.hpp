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

#ifndef READCOUNTS_HPP
#define READCOUNTS_HPP

#include "kraken_headers.hpp"
#include "hyperloglogplus.h"

namespace kraken {
  static size_t HLL_PRECISION = 12;
  struct ReadCounts {
    uint64_t n_reads;
    uint64_t n_kmers;
    HyperLogLogPlusMinus<uint64_t> kmers; // unique k-mer count per taxon

    ReadCounts() : n_reads(0), n_kmers(0), kmers(HyperLogLogPlusMinus<uint64_t>(HLL_PRECISION)) { }

    ReadCounts(size_t precision) : kmers(HyperLogLogPlusMinus<uint64_t>(precision)) {
    }
    
    void add_kmer(uint64_t kmer) {
      ++ n_kmers;
      kmers.add(kmer);
    }
    
    ReadCounts& operator+=(const ReadCounts& b) {
      n_reads += b.n_reads;
      n_kmers += b.n_kmers;
      kmers += b.kmers;
      return *this;
    }

    bool operator<(const ReadCounts& rc) {
      if (n_reads < rc.n_reads) {
        return true;
      }
      if (n_reads == rc.n_reads && n_kmers < rc.n_kmers) {
        return true;
      }
      return false;
    }
  };
  
  uint64_t reads(const ReadCounts& read_count) {
    return(read_count.n_reads);
  }
}
#endif
