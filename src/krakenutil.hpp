/*
 * Original file Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
 * Portions (c) 2017-2018, Florian Breitwieser <fbreitwieser@jhu.edu> as part of KrakenUniq
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

#ifndef KRAKENUTIL_HPP
#define KRAKENUTIL_HPP

#include "kraken_headers.hpp"
#include <unordered_map>

namespace kraken {
  // Build a map of node to parent from an NCBI taxonomy nodes.dmp file
  std::unordered_map<uint32_t, uint32_t> build_parent_map(std::string filename);

    // Return lowest common ancestor of a and b
    // LCA(0,x) = LCA(x,0) = x
    // Default ancestor is 1 (root of tree)
  uint32_t lca(const std::unordered_map<uint32_t, uint32_t> &parent_map, uint32_t a, uint32_t b);
  inline bool ends_with(std::string const &value, std::string const &ending);

  class managed_ostream {
  public:
    managed_ostream &operator=(const managed_ostream &) = delete;
    managed_ostream(const managed_ostream &) = delete;

    managed_ostream();
    managed_ostream(std::ostream *managed_ostream);
    managed_ostream(const std::string &file, bool use = true, bool append = false);
    ~managed_ostream();
    std::ostream &operator*();
    std::ostream *operator->();

  private:
    std::ostream *s;
  };

  // Resolve classification tree
  uint32_t resolve_tree(const std::unordered_map<uint32_t, uint32_t> &hit_counts,
                        const std::unordered_map<uint32_t, uint32_t> &parent_map);

  class KmerScanner {
    public:

    KmerScanner(std::string &seq, size_t start=0, size_t finish=~0);
    uint64_t *next_kmer();  // NULL when seq exhausted
    bool ambig_kmer();  // does last returned kmer have non-ACGT?


    static uint8_t get_k();
    // MUST be called before first invocation of KmerScanner()
    static void set_k(uint8_t n);

    private:
    std::string *str;
    size_t curr_pos, pos1, pos2;
    uint64_t kmer;  // the kmer, address is returned (don't share b/t thr.)
    uint32_t ambig; // is there an ambiguous nucleotide in the kmer?
    int64_t loaded_nt;

    static uint8_t k;  // init. to 0 b/c static
    static uint64_t kmer_mask;
    static uint32_t mini_kmer_mask;
  };
}

#endif // ndef KRAKENUTIL_HPP
