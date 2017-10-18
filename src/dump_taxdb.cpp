#include "taxdb.h"
#include "quickfile.hpp"
#include <iostream>
#include <fstream>
#include <unordered_map>

using namespace std;

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "Usage: dump_taxdb taxDB names.dmp nodes.dmp\n";
    return 1;
  }

	cerr << "Reading taxonomy database from " << argv[1] << ", writing nodes dump to " << argv[3] << " and names dump to " << argv[2] << "." << endl;
  TaxonomyDB<uint32_t, uint32_t> taxdb {(string)argv[1]};
  ofstream names_file(argv[2]);
  names_file.exceptions(ifstream::failbit | ifstream::badbit);
  ofstream nodes_file(argv[3]);
  nodes_file.exceptions(ifstream::failbit | ifstream::badbit);

  for (const auto &taxon : taxdb.taxIDsAndEntries) {
    std::string scientificName;
    nodes_file << taxon.second.taxonomyID 
      << "\t|\t" << taxon.second.parentTaxonomyID
      << "\t|\t" << taxon.second.rank
      << endl; // there are further columns, but Kraken does not care about them
    
    names_file << taxon.second.taxonomyID 
      << "\t|\t" << taxon.second.scientificName
      << "\t|\t" 
      << "\t|\t" << "scientific name" << endl;
  }
  names_file.close();
  nodes_file.close();
}
