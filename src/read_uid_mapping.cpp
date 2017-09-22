
#include "kraken_headers.hpp"
#include "quickfile.hpp"
#include <unordered_map>
#include <algorithm>

using namespace std;
using namespace kraken;

inline
vector<uint32_t> get_taxids_for_uid(uint32_t uid, char* fptr) {
  size_t int_size = sizeof(int);
  size_t block_size = sizeof(int)*2;
  // TODO: Just get a uint64_t and shift the bits, probably faster
  uint32_t taxid  = *(uint32_t*)(fptr+(uid-1)*block_size);
  uint32_t parent_uid = *(uint32_t*)(fptr+(uid-1)*block_size + int_size);
  
  vector<uint32_t> taxids = {taxid};
  while (parent_uid != 0) {
    taxid  = *(uint32_t*)(fptr+(parent_uid-1)*block_size);
    parent_uid = *(uint32_t*)(fptr+(parent_uid-1)*block_size + int_size);
    taxids.push_back(taxid);
  }
  std::sort(taxids.begin(), taxids.end());
  return(taxids);
}

inline
vector<uint32_t> get_taxids_for_uid_from_map(uint32_t uid, char* fptr, unordered_map<uint32_t, vector<uint32_t> >& uid_map ) {
  auto it = uid_map.find(uid);
  if (it != uid_map.end()) {
    return it->second;
  } 
  vector<uint32_t> taxids = get_taxids_for_uid(uid, fptr);
  uid_map[uid] = taxids;
  return(taxids);
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: read_uid_mapping <uid mappingfile> [<uid>]"
        "The file is supposed to have lines terminated by '\n'."
         << std::endl;
    return 1;
  }
  char *filename = argv[1];
  kraken::QuickFile UID_to_TaxID_map_file;
  UID_to_TaxID_map_file.open_file(filename);

  char* fptr = UID_to_TaxID_map_file.ptr();
  if (argc == 2) {
    vector< vector <uint32_t> > UIDs_to_taxids;
    uint32_t UID = 1;
    size_t int_size = sizeof(UID);
    size_t i = 0;
    for (size_t pos = 0; pos < UID_to_TaxID_map_file.size(); pos += 2*int_size) {
      uint32_t* taxid_ptr  = (uint32_t*)(fptr+pos);
      uint32_t* parent_uid = (uint32_t*)(fptr+pos+int_size);
      //UIDs_to_taxids.push_back( { UIDs_to_taxids[] } );
      //pos += int_size;
      cout << ++i << '\t' << *taxid_ptr << '\t' << *parent_uid << endl;
    }
  } else {
    unordered_map<uint32_t, vector<uint32_t> > UID_to_TaxID_map;
    for (int i=2; i <argc; ++i) {
      uint32_t UID = atol(argv[i]);
      vector<uint32_t> taxids = get_taxids_for_uid(UID, UID_to_TaxID_map, fptr);
      cout << UID << '\t';
      for (auto t : taxids) {
        cout << t << ' ';
      }
      cout << endl;
    }
  }

  return 0;
}
