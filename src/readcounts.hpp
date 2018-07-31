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
#include "hyperloglogplus.hpp"

namespace kraken {
  static size_t HLL_PRECISION = 14;

  template <typename CONTAINER>
  class ReadCounts {

  public:
    uint64_t readCount() const { return n_reads; }
    void incrementReadCount() { ++n_reads; }
    uint64_t kmerCount() const { return n_kmers; }
    uint64_t uniqueKmerCount() const; // to be implemented for each CONTAINER

    ReadCounts() : n_reads(0), n_kmers(0) {
    }

    ReadCounts(const ReadCounts& other) : n_reads(other.n_reads), n_kmers(other.n_kmers), kmers(other.kmers) {
    }

    ReadCounts& operator=(const ReadCounts& other) {
      n_reads = other.n_reads;
      n_kmers =other.n_kmers;
      kmers = other.kmers;
      return *this;
    }

    ReadCounts& operator=(ReadCounts&& other) {
      n_reads = other.n_reads;
      n_kmers =other.n_kmers;
      kmers = std::move(other.kmers);
      return *this;
    }
  
    void add_kmer(uint64_t kmer) {
      ++n_kmers;
      kmers.insert(kmer);
    }

    ReadCounts& operator+=(const ReadCounts& other) {
      n_reads += other.n_reads;
      n_kmers += other.n_kmers;
      kmers += other.kmers;
      return *this;
    }
/*
    ReadCounts& operator+=(ReadCounts&& other) {
      n_reads += other.n_reads;
      n_kmers += other.n_kmers;
      kmers += std::move(other.kmers);
      return *this;
    }
*/

    bool operator<(const ReadCounts& other) {
      if (n_reads < other.n_reads) {
        return true;
      }
      if (n_reads == other.n_reads && n_kmers < other.n_kmers) {
        return true;
      }
      return false;
    }


  private:
    uint64_t n_reads;
    uint64_t n_kmers; 
    CONTAINER kmers; // unique k-mer count per taxon

  };

  // Overload operator += for set, so that it can be used for merging
  template <typename T>
  unordered_set<T>& operator+=(unordered_set<T>& left, const unordered_set<T>& right) {
    left.insert(right.begin(), right.end());
    return left;
  }
  
  template<>
  uint64_t ReadCounts< HyperLogLogPlusMinus<uint64_t> >::uniqueKmerCount() const {
    return(kmers.cardinality());
  }

  template<>
  uint64_t ReadCounts< unordered_set<uint64_t> >::uniqueKmerCount() const {
    return(kmers.size());
  }
}
#endif
