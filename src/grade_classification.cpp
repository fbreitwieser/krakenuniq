/*
 * Copyright 2017, Florian Breitwieser
 * licnsed under GPLv3
 */

#include "taxdb.h"
#include "quickfile.hpp"
#include <iostream>
#include <fstream>
#include <unordered_map>

using namespace std;

unordered_map<string, uint32_t> read_seqid_mapping(string filename) {
  unordered_map<string, uint32_t> ID_to_taxon_map;
  ifstream map_file(filename.c_str());
  if (map_file.rdstate() & ifstream::failbit) {
    err(EX_NOINPUT, "can't open %s", filename.c_str());
  }
  string line, seq_id;
  uint32_t taxid;

  while (map_file.good()) {
    getline(map_file, line);
    if (line.empty())
      break;
    istringstream iss(line);
    iss >> seq_id >> taxid;
    ID_to_taxon_map[seq_id] = taxid;
  }
  map_file.close();
  return ID_to_taxon_map;
}

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "Usage: grade_classification taxDB seqid2taxid.map classification_file\n";
    return 1;
  }
  TaxonomyDB<uint32_t, uint32_t> taxdb = TaxonomyDB<uint32_t, uint32_t>(argv[1], false);
  unordered_map<string, uint32_t> seqid_map = read_seqid_mapping(argv[2]);
  cerr << "Read " << seqid_map.size() << " taxa mappings" << endl;
  
  ifstream k_file(argv[3]);
  if (k_file.rdstate() & ifstream::failbit) {
    err(EX_NOINPUT, "can't open %s", argv[3]);
  }
  string line, classification_state, read_id, seq_id;
  uint32_t taxid;
  uint32_t seq_taxid;

  while (k_file.good()) {
    getline(k_file, line);
    if (line.empty())
      continue;
    istringstream iss(line);
    iss >> classification_state >> read_id >> taxid;
    seq_id = read_id.substr(read_id.find_last_of("_")+1);
    auto it = seqid_map.find(seq_id);
    if (it == seqid_map.end()) {
      cerr << "ERROR: Couldn't find taxid for " << seq_id << endl;
    } else {
      seq_taxid = it->second;
      size_t distance_between_taxids;
      string lowest_common_rank;
      seq_taxid = taxdb.getTaxIDAtRank(seq_taxid, "species");
      taxid = taxdb.getTaxIDAtRank(taxid, "species");
      pair<uint32_t, int> lca_taxid_dist = taxdb.getLowestCommonAncestor(seq_taxid, taxid);
      string lca_rank = taxdb.getRank(lca_taxid_dist.first);
      cout << seq_taxid << '\t'  << taxid << '\t' << lca_rank << '\t' << lca_taxid_dist.first << '\t' << lca_taxid_dist.second << endl;
    }
  }
  k_file.close();


}
