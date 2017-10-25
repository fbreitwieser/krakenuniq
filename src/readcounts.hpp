
#ifndef READCOUNTS_HPP
#define READCOUNTS_HPP

#include "kraken_headers.hpp"
#include "hyperloglogplus.h"

namespace kraken {
  struct ReadCounts {
    uint64_t n_reads = 0;
    uint64_t n_kmers = 0;
    HyperLogLogPlusMinus<uint64_t> kmers; // unique k-mer count per taxon

    ReadCounts() { }

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
  };
  
  uint64_t reads(const ReadCounts& read_count) {
    return(read_count.n_reads);
  }
}
#endif
