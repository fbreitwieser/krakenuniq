/*
 * Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
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

#include "assert_helpers.h"
#include "kraken_headers.hpp"
#include "krakenutil.hpp"

using namespace std;

namespace kraken {

  // Build a node->parent unordered_map from NCBI Taxonomy nodes.dmp file
  unordered_map<uint32_t, uint32_t> build_parent_map(string filename) {
    unordered_map<uint32_t, uint32_t> pmap;
    uint32_t node_id, parent_id;
    string line;
    ifstream ifs(filename.c_str());
    if (ifs.rdstate() & ifstream::failbit) {
      err(EX_NOINPUT, "error opening %s", filename.c_str());
    }

    while (ifs.good()) {
      getline(ifs, line);
      if (line.empty())
        break;
      sscanf(line.c_str(), "%d\t|\t%d", &node_id, &parent_id);
      pmap[node_id] = parent_id;
    }
    pmap[1] = 0;
    return pmap;
  }

  // Return lowest common ancestor of a and b
  // LCA(0,x) = LCA(x,0) = x
  // Default ancestor is 1 (root of tree)
  uint32_t lca(const unordered_map<uint32_t, uint32_t> &parent_map,
    uint32_t a, uint32_t b)
  {
    if (a == 0 || b == 0)
      return a ? a : b;

    // create a path from a to the root
    set<uint32_t> a_path;
    while (a > 0) {
      a_path.insert(a);
      assert(parent_map.find(a) != parent_map.end());
      a = parent_map.at(a);
    }

    // search for b in the path from a to the root
    while (b > 0) {
      if (a_path.count(b) > 0)
        return b;
      assert(parent_map.find(b) != parent_map.end());
      b = parent_map.at(b);
    }
    return 1;
  }

  // Tree resolution: take all hit taxa (plus ancestors), then
  // return leaf of highest weighted leaf-to-root path.
  uint32_t resolve_tree(const unordered_map<uint32_t, uint32_t> &hit_counts,
                        const unordered_map<uint32_t, uint32_t> &parent_map)
  {
    set<uint32_t> max_taxa;
    uint32_t max_taxon = 0, max_score = 0;

    // Sum each taxon's LTR path
    for (auto it = hit_counts.begin();
         it != hit_counts.end(); ++it) {
      uint32_t taxon = it->first;
      uint32_t node = taxon;
      uint32_t score = 0;
      while (node > 0) {
        auto it2 = hit_counts.find(node);
        if (it2 != hit_counts.end()) {
          score += it2->second;
        }
        node = parent_map.at(node);

      }

      if (score > max_score) {
        max_taxa.clear();
        max_score = score;
        max_taxon = taxon;
      }
      else if (score == max_score) {
        if (max_taxa.empty())
          max_taxa.insert(max_taxon);
        max_taxa.insert(taxon);
      }
    }

    // If two LTR paths are tied for max, return LCA of all
    if (! max_taxa.empty()) {
      set<uint32_t>::iterator sit = max_taxa.begin();
      max_taxon = *sit;
      for (sit++; sit != max_taxa.end(); sit++)
        max_taxon = lca(parent_map, max_taxon, *sit);
    }

    return max_taxon;
  }


  // Tree resolution: take all hit taxa (plus ancestors), then
  // return leaf of highest weighted leaf-to-root path.
  uint32_t resolve_uids(
      const unordered_map<uint32_t, uint32_t> &uid_hit_counts,
      const unordered_map<uint32_t, uint32_t> &parent_map,
      const vector< vector<uint32_t> > &UID_to_taxids_vec) {
    unordered_map<uint32_t, uint32_t> taxid_counts;
    unordered_map<uint32_t, double> frac_taxid_counts;

    if (uid_hit_counts.size() == 0) {
      return(0);
    }

    for (auto it = uid_hit_counts.begin(); it != uid_hit_counts.end(); ++it) {
      uint32_t uid = it->first;
      double frac_count = ((double)it->second / (double)UID_to_taxids_vec[uid-1].size());
      for (auto taxid : UID_to_taxids_vec[uid-1]) {
        taxid_counts[taxid] += it->second;
        frac_taxid_counts[taxid] += frac_count;
      }
    }
    vector<uint32_t> max_taxids;
    uint32_t max_count = 0;
    double max_frac_count = 0;
    for (auto it : taxid_counts) {
      if (it.second == max_count) {
        if (frac_taxid_counts[it.first] == max_frac_count) {
          max_taxids.push_back(it.first);
        } else if (frac_taxid_counts[it.first] > max_frac_count) {
          max_frac_count = frac_taxid_counts[it.first];
          max_taxids = { it.first };
        }
      } else if (it.second > max_count) {
        max_taxids = { it.first };
        max_count = it.second;
        max_frac_count = frac_taxid_counts[it.first];
      }
    }

    uint32_t max_taxon = max_taxids[0];
    auto sit = max_taxids.begin();
    for (++sit; sit != max_taxids.end(); ++sit) {
      max_taxon = lca(parent_map, max_taxon, *sit);

    }

    // return the taxid that appeared most often
    return max_taxon;
  }

  // Tree resolution: take all hit taxa (plus ancestors), then
  // return leaf of highest weighted leaf-to-root path.
  uint32_t resolve_uids2(
      const unordered_map<uint32_t, uint32_t> &uid_hit_counts,
      const unordered_map<uint32_t, uint32_t> &parent_map,
      char* fptr) {
    unordered_map<uint32_t, uint32_t> taxid_counts;
    unordered_map<uint32_t, double> frac_taxid_counts;

    if (uid_hit_counts.size() == 0) {
      return(0);
    }

    size_t int_size = sizeof(int);
    size_t block_size = sizeof(int)*2;
    for (auto it = uid_hit_counts.begin(); it != uid_hit_counts.end(); ++it) {
      uint32_t uid = it->first;
      if (uid == 0) {
	continue;
      }
      uint32_t taxid;
      // TODO: Just get a uint64_t and shift the bits, probably faster
      vector<uint32_t> taxids;
      do {
        taxid = *(uint32_t*)(fptr+(uid-1)*block_size);
        uid = *(uint32_t*)(fptr+(uid-1)*block_size + int_size);
  
        taxid_counts[taxid] += it->second;
	taxids.push_back(taxid);
      } while (uid != 0);

      double frac_count = (double)it->second / (double)taxids.size();
      for (uint32_t taxid : taxids) {
        frac_taxid_counts[taxid] += frac_count;
      }
    }

    if (taxid_counts.size() == 0) {
      return(0);
    }
    vector<uint32_t> max_taxids;
    uint32_t max_count = 0;
    double max_frac_count = 0;
    for (auto it : taxid_counts) {
      if (it.second == max_count) {
        if (frac_taxid_counts[it.first] == max_frac_count) {
          max_taxids.push_back(it.first);
        } else if (frac_taxid_counts[it.first] > max_frac_count) {
          max_frac_count = frac_taxid_counts[it.first];
          max_taxids = { it.first };
        }
      } else if (it.second > max_count) {
        max_taxids = { it.first };
        max_count = it.second;
        max_frac_count = frac_taxid_counts[it.first];
      }
    }

    uint32_t max_taxon = max_taxids[0];
    auto sit = max_taxids.begin();
    for (++sit; sit != max_taxids.end(); ++sit) {
      max_taxon = lca(parent_map, max_taxon, *sit);

    }

    // return the taxid that appeared most often
    return max_taxon;
  }




  uint8_t KmerScanner::k = 0;
  uint64_t KmerScanner::kmer_mask = 0;
  uint32_t KmerScanner::mini_kmer_mask = 0;

  // Create a scanner for the string over the interval [start, finish)
  KmerScanner::KmerScanner(string &seq, size_t start, size_t finish) {
    if (! k)
      errx(EX_SOFTWARE, "KmerScanner created w/o setting k");
    if (finish > seq.size())
      finish = seq.size();

    kmer = 0;
    ambig = 0;
    str = &seq;
    curr_pos = start;
    pos1 = start;
    pos2 = finish;
    loaded_nt = 0;
    if (pos2 - pos1 + 1 < k)
      curr_pos = pos2;
  }

  uint8_t KmerScanner::get_k() { return k; }

  void KmerScanner::set_k(uint8_t n) {
    if (k)  // Only allow one setting per execution
      return;
    k = n;
    kmer_mask = ~0;
    kmer_mask >>= sizeof(kmer_mask) * 8 - (k * 2);
    mini_kmer_mask = ~0;
    mini_kmer_mask >>= sizeof(mini_kmer_mask) * 8 - k;
  }

  uint64_t *KmerScanner::next_kmer() {
    if (curr_pos >= pos2)
      return NULL;
    if (loaded_nt)  
      loaded_nt--;
    while (loaded_nt < k) {
      loaded_nt++;
      kmer <<= 2;
      ambig <<= 1;
      switch ((*str)[curr_pos++]) {
        case 'A': case 'a':
          break;
        case 'C': case 'c':
          kmer |= 1;
          break;
        case 'G': case 'g':
          kmer |= 2;
          break;
        case 'T': case 't':
          kmer |= 3;
          break;
        default:
          ambig |= 1;
          break;
      }
      kmer &= kmer_mask;
      ambig &= mini_kmer_mask;
    }
    return &kmer;
  }

  bool KmerScanner::ambig_kmer() {
    return !! ambig;
  }
}
