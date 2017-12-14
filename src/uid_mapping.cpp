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

#include<iostream>
#include<algorithm>
#include "uid_mapping.hpp"
#include "krakenutil.hpp"
#include "assert_helpers.h"

using namespace std;

namespace kraken {

  //static size_t INT_SIZE=sizeof(uint32_t);
  //static size_t UID_BLOCK_SIZE=2*INT_SIZE;
  static uint32_t max_uid = -1;

  uint32_t uid_mapping(
      map< TaxidSet, uint32_t>& Taxids_to_UID_map, 
      vector< const TaxidSet*  >& UID_to_taxids_vec, 
      uint32_t taxid, 
      uint32_t kmer_uid, 
      uint32_t& current_uid,
      ofstream& UID_map_file) {

    vector<uint32_t> taxid_set;
    if (kmer_uid == 0) {
      taxid_set.push_back(taxid);
    } else {
      if (kmer_uid > UID_to_taxids_vec.size()) {
        // This can happen when set_lcas is called more than once on a database (ie not all values start w/ 0)
        cerr << "kmer_uid ("<< kmer_uid <<") greater than UID vector size ("<< UID_to_taxids_vec.size()<<")!!" << endl;
        exit(1);
      }
      taxid_set = *(UID_to_taxids_vec[kmer_uid-1]);
      auto it = std::lower_bound( taxid_set.begin(), taxid_set.end(), taxid); // find proper position in descending order
      if (it == taxid_set.end() || *it != taxid) {
        // add the taxid to the set, in the right position such that it remains sorted
         taxid_set.insert( it, taxid ); // insert before iterator it
      } else {
        // the taxid is already part of the set for kmer_uid, return kmer_uid
        return kmer_uid;
      }
    }

    // This taxid is not part of kmer_uids set, but is this new taxon_set already assigned to another UID?
    // Try inserting ..
    auto insert_res = Taxids_to_UID_map.insert( { std::move(taxid_set), current_uid + 1 } );
    if (!insert_res.second) {
      // Insert unsuccessful, taxid set already has an UID
      return insert_res.first->second;
    }

    // Get a new UID
    if (max_uid <= ++current_uid) {
      cerr << "Maxxed out on UIDs!!" << endl;
      exit(1);
    }

    UID_to_taxids_vec.push_back( &(insert_res.first->first) );
    assert(UID_to_taxids_vec.size() == current_uid);

    // Write to mapping file
    // format: TAXID<uint32_t> PARENT<uint32_t>
    // read it with read_uid_mapping
    UID_map_file.write((char*)&taxid, sizeof(taxid));
    UID_map_file.write((char*)&kmer_uid, sizeof(kmer_uid));

    return current_uid;
  } // end of uid_mapping


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
      //for (auto taxid : UID_to_taxids_vec[uid-1]) {
      for (auto taxid_it = UID_to_taxids_vec[uid-1].begin(); taxid_it != UID_to_taxids_vec[uid-1].end(); ++taxid_it) { // supporting gcc 4.4
        taxid_counts[*taxid_it] += it->second;
        frac_taxid_counts[*taxid_it] += frac_count;
      }
    }
    vector<uint32_t> max_taxids;
    uint32_t max_count = 0;
    double max_frac_count = 0;
    // for (auto it : taxid_counts) { 
    for (auto itt = taxid_counts.begin(); itt != taxid_counts.end(); ++itt) {  // supporting gcc 4.4
      const auto& it = *itt;
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
      const char* fptr, const size_t fsize) {

    (void)fsize; // parameter unused

    unordered_map<uint32_t, uint32_t> taxid_counts;
    unordered_map<uint32_t, double> frac_taxid_counts;

    if (uid_hit_counts.size() == 0) {
      return(0);
    }

    for (auto it1=uid_hit_counts.begin(); it1 != uid_hit_counts.end(); ++it1) { // supporting gcc 4.4
      const auto &it = *it1;
      if (it.first == 0) {
        continue;
      }
      // TODO: Just get a uint64_t and shift the bits, probably faster
      vector<uint32_t> taxids = get_taxids_for_uid(it.first, fptr);

      double frac_count = (double)it.second / (double)taxids.size();
      for (size_t i = 0; i < taxids.size(); ++i) { // supporting gcc 4.4
        uint32_t taxid = taxids[i];
        frac_taxid_counts[taxid] += frac_count;
        taxid_counts[taxid] += it.second;
      }
    }

    if (taxid_counts.size() == 0) {
      return(0);
    }
    vector<uint32_t> max_taxids;
    uint32_t max_count = 0;
    double max_frac_count = 0;
    for (auto it1 = taxid_counts.begin(); it1 != taxid_counts.end(); ++it1) {
      const auto& it = *it1;
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
  // This version saves observed mappings in a map
  // Doesn't seem to give big runtime improvements I've hope for, so far 
  uint32_t resolve_uids3(
      const unordered_map<uint32_t, uint32_t> &uid_hit_counts,
      const unordered_map<uint32_t, uint32_t> &parent_map,
      unordered_map<uint32_t, vector<uint32_t> > &uid_dict,
      const char* fptr, const size_t fsize) {

    (void)fsize; // parameter unused

    unordered_map<uint32_t, uint32_t> taxid_counts;
    unordered_map<uint32_t, double> frac_taxid_counts;

    if (uid_hit_counts.size() == 0) {
      return(0);
    }

    for (auto it1=uid_hit_counts.begin(); it1 != uid_hit_counts.end(); ++it1) { // supporting gcc 4.4
      const auto &it = *it1;
      if (it.first == 0) {
        continue;
      }
      // TODO: Just get a uint64_t and shift the bits, probably faster
      vector<uint32_t> taxids = get_taxids_for_uid_from_map(it.first, fptr, uid_dict);

      double frac_count = (double)it.second / (double)taxids.size();
      for (size_t i = 0; i < taxids.size(); ++i) { // supporting gcc 4.4
        uint32_t taxid = taxids[i];
        frac_taxid_counts[taxid] += frac_count;
        taxid_counts[taxid] += it.second;
      }
    }

    if (taxid_counts.size() == 0) {
      return(0);
    }
    vector<uint32_t> max_taxids;
    uint32_t max_count = 0;
    double max_frac_count = 0;
    for (auto it1 = taxid_counts.begin(); it1 != taxid_counts.end(); ++it1) {
      const auto& it = *it1;
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

}

vector<uint32_t> get_taxids_for_uid(const uint32_t uid, const char* fptr) {
  size_t int_size = sizeof(int);
  size_t block_size = sizeof(int)*2;
  // TODO: Just get a uint64_t and shift the bits, probably faster
  uint32_t taxid  = *(uint32_t*)(fptr+(uid-1)*block_size);
  uint32_t parent_uid = *(uint32_t*)(fptr+(uid-1)*block_size + int_size);

  vector<uint32_t> taxids = {taxid};
  while (parent_uid != 0) {
    // TODO: Consider checking if the accessed meory is out of range. 
      // assert((parent_uid-1)*block_size <= fsize);
      // if (offset >= fsize) {
      //   cerr << "It seems you are trying to access a block after the file end: \n" <<
      //      " fptr: " << fptr << "; uid: " << next_uid << "; " << " addr: " << (offset + INT_SIZE) << endl;
      //  exit(1);
      //}
    taxid  = *(uint32_t*)(fptr+(parent_uid-1)*block_size);
    parent_uid = *(uint32_t*)(fptr+(parent_uid-1)*block_size + int_size);
    taxids.push_back(taxid);
  }
  //std::sort(taxids.begin(), taxids.end());
  return(taxids);
}

vector<uint32_t> get_taxids_for_uid_from_map(const uint32_t uid, const char* fptr, unordered_map<uint32_t, vector<uint32_t> >& uid_map ) {
  auto it = uid_map.find(uid);
  if (it != uid_map.end()) {
    return it->second;
  } 
  vector<uint32_t> taxids = get_taxids_for_uid(uid, fptr);
  uid_map[uid] = taxids;
  return(taxids);
}

