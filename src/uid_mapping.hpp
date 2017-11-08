
#ifndef UID_MAPPING_H
#define UID_MAPPING_H

#include<vector>
#include<map>
#include<unordered_map>
#include<fstream>
using namespace std;


// Takes the current UID kmer_uid, and checks whether
// - taxid is in taxon set T specified in UID_to_taxids_vec[kmer_uid]?
//   - yes: return kmer_uid
//   - no: is there a set (T,taxid) in Taxids_to_UID_map?
//     - yes: return the uid of that set
//     - no:  
//       - increment current_uid by one and set this as the set uid
//       - add the set to Taxids_to_UID_map and UID_to_taxids_vec
//       - write the mapping to UID_map_file
//

//using TaxidSet = typename std::vector<uint32_t>;
typedef std::vector<uint32_t> TaxidSet;

namespace kraken {


uint32_t uid_mapping(
      map< TaxidSet, uint32_t>& Taxids_to_UID_map, 
      vector< const TaxidSet* >& UID_to_taxids_vec, 
      uint32_t taxid, 
      uint32_t kmer_uid, 
      uint32_t& current_uid,
      ofstream& UID_map_file);


uint32_t resolve_uids(
      const unordered_map<uint32_t, uint32_t> &uid_hit_counts,
      const unordered_map<uint32_t, uint32_t> &parent_map,
      const vector< vector<uint32_t> > &UID_to_taxids_vec);

uint32_t resolve_uids2(
      const unordered_map<uint32_t, uint32_t> &uid_hit_counts,
      const unordered_map<uint32_t, uint32_t> &parent_map,
      const char* fptr, const size_t fsize);
}

vector<uint32_t> get_taxids_for_uid(const uint32_t uid, const char* fptr);

vector<uint32_t> get_taxids_for_uid_from_map(uint32_t uid, char* fptr, unordered_map<uint32_t, vector<uint32_t> >& uid_map );

#endif
