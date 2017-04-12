/* Original work Copyright 2013 David Ainsworth
 * Modified work copyright 2017 Florian Breitwieser
 *
 * The original file is part of SLAM
 * The modified file is part of a modified Kraken version
 *
 * SLAM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SLAM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.

 * You should have received a copy of the GNU Affero General Public License
 * along with SLAM.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef TAXD_DB_H_
#define TAXD_DB_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include "report-cols.h"

using namespace std;

void log_msg (const std::string& s);

template<typename T> uint64_t string_to_T(std::string str);

template <typename T> 
inline uint64_t reads(const T read_count);

inline uint64_t reads(const uint64_t read_count);

std::vector<std::string> in_betweens(const std::string &s, const char start_char, const char end_char, size_t start_at = 0);

std::vector<std::string> tokenise(const std::string &s, const std::string& delimiter, size_t max_fields = 0, size_t end_chars = 0);


std::vector<std::string> get_fields(const std::string &s, const std::string& delimiter, std::vector<size_t> fields); 

template<typename TAXID, typename READCOUNTS>
class TaxonomyEntry {
 public:
  TAXID taxonomyID = 0;
  TAXID parentTaxonomyID = 0;
  std::string rank;
  std::string scientificName;

  TaxonomyEntry() {}

  TaxonomyEntry(TAXID taxonomyID_, std::string scientificName_) :
	  taxonomyID(taxonomyID_), scientificName(scientificName_) {}

  TaxonomyEntry(TAXID taxonomyID_, TAXID parentTaxonomyID_, std::string rank_) :
	  taxonomyID(taxonomyID_), parentTaxonomyID(parentTaxonomyID_), rank(rank_) {}

  TaxonomyEntry(TAXID taxonomyID_, TAXID parentTaxonomyID_, std::string rank_, std::string scientificName_) :
	  taxonomyID(taxonomyID_), parentTaxonomyID(parentTaxonomyID_), rank(rank_), scientificName(scientificName_) {}

  inline bool operator==(const TaxonomyEntry& other) const; 
  TaxonomyEntry* parent = nullptr;
  std::vector<TaxonomyEntry*> children;

  READCOUNTS read_counts = READCOUNTS();
  READCOUNTS read_counts_children = READCOUNTS();

  bool used = false;
  uint64_t genomeSize = 0;
  uint64_t genomeSizeOfChildren = 0;
  uint64_t numBelow = 0;
};

//template<>
//TaxonomyEntry<uint32_t, uint64_t>::TaxonomyEntry () {
//	read_counts = 0;
//	read_counts_children = 0;
//}

template<typename TAXID, typename READCOUNTS>
struct TaxonomyEntryPtr_comp {
	bool operator() ( const TaxonomyEntry<TAXID,READCOUNTS>* a, const TaxonomyEntry<TAXID,READCOUNTS>* b) const;
};


template<typename TAXID, typename READCOUNTS>
class TaxonomyDB {
 public:
  TaxonomyDB(const std::string inFileName);
  TaxonomyDB();
  //std::unordered_map<std::string, TAXID> seqIDsAndTaxIds;
  std::unordered_map<TAXID, TaxonomyEntry<TAXID,READCOUNTS> > taxIDsAndEntries;
  void parseNamesDump(const std::string namesDumpFileName);
  void parseNodesDump(const std::string nodesDumpFileName);
  TAXID getTaxIDAtRank(const TAXID taxID, const std::string& rank) const;
  std::string getScientificName(const TAXID taxID) const;
  std::string getRank(const TAXID taxID) const;
  TAXID getLowestCommonAncestor(const std::vector<TAXID>& taxIDs) const;
  TAXID getParentTaxID(const TAXID taxID) const;
  std::unordered_map<TAXID, TAXID> getParentMap() const;
  std::unordered_map<std::string, TAXID> getScientificNameMap() const;
  std::string getLineage(TAXID taxonomyID) const;
  std::string getMetaPhlAnLineage(TAXID taxonomyID) const;
  char* getIndexFileName(const TAXID hostTaxID) const;
  void readTaxonomyIndex(const std::string inFileName);
  void writeTaxonomyIndex(std::ostream & outs) const;
  void writeTaxonomyIndex(std::ostream & outs,
                          const std::string namesDumpFileName,
                          const std::string nodesDumpFileName);
  bool isSubSpecies(TAXID taxonomyID) const;
  int isBelowInTree(TAXID upper, TAXID lower) const;
  void addCounts(const TAXID taxid, const READCOUNTS& read_counts_);
  void fillCounts(const std::unordered_map<TAXID, READCOUNTS>& taxon_counts);
  void createPointers();
  void printReport();
};


template<typename TAXID, typename READCOUNTS>
class TaxReport {
private:
	std::ostream& _reportOfb;
	TaxonomyDB<TAXID,READCOUNTS> & _taxdb;
	std::vector<REPORTCOLS> _report_cols;
	uint64_t _total_n_reads;
	bool _show_zeros;

	void printLine(TaxonomyEntry<TAXID,READCOUNTS>& tax, unsigned depth);

public:
	TaxReport(std::ostream& _reportOfb, TaxonomyDB<TAXID,READCOUNTS> & taxdb, bool _show_zeros);

	void printReport(std::string format, std::string rank);
	void printReport(TaxonomyEntry<TAXID,READCOUNTS>& tax, unsigned depth);
};


  // Return lowest common ancestor of a and b
  // LCA(0,x) = LCA(x,0) = x
  // Default ancestor is 1 (root of tree)
uint32_t lca(std::unordered_map<uint32_t, uint32_t> &parent_map, uint32_t a, uint32_t b);

template<typename K,typename V>
inline
V find_or_use_default(const std::unordered_map<K, V>& my_map, const K& query, const V default_value);

//////////////////////////// DEFINITIONS
void log_msg (const std::string& s) {
	std::cerr << s << "\n";
}

template<typename T>
uint64_t string_to_T(string str) {
  stringstream stream(str);
  T result;
  stream >> result;
  return result;
}

template <typename T>
inline
uint64_t reads(const T read_count) {
	cerr << "No reads function for type!! " << endl;
	throw ;
	return(0);
}



inline
uint64_t reads(const uint64_t read_count) {
	return(read_count);
}

std::vector<std::string> in_betweens(const std::string &s, const char start_char, const char end_char, size_t start_at) {
    std::vector<std::string> tokens;
	size_t i = 0;
	size_t next_end = start_at-1;
    
	for (size_t next_start = s.find(start_char, next_end + 1); \
		 next_start != string::npos;
         next_start = s.find(start_char, next_end + 1), ++i) {

		next_end = s.find(end_char, next_start + 1);
		if (next_end == string::npos)
			throw std::runtime_error("unmatched start and end!");

        tokens.push_back(s.substr(next_start+1, next_end-1));
    }

    return tokens;
}



std::vector<std::string> tokenise(const std::string &s, const std::string& delimiter, size_t max_fields, size_t end_chars) {
    std::vector<std::string> tokens(max_fields);
    size_t delim_length = delimiter.length();
    size_t last = 0;
    size_t i = 0;

    for (size_t next = s.find(delimiter, last);
         (max_fields > 0 && i < max_fields) && next != string::npos;
         next = s.find(delimiter, last), ++i) {
        tokens[i] = s.substr(last, next-last);
        last = next + delim_length;
    }
    if (max_fields > 0 && i < max_fields) {
        tokens[max_fields-1] = s.substr(last, s.length()-last-end_chars);
    }

    return tokens;
}

std::vector<std::string> get_fields(const std::string &s, const std::string& delimiter, vector<size_t> fields) {
    std::vector<std::string> tokens;
	tokens.reserve(fields.size());
    size_t delim_length = delimiter.length();
    size_t last = 0;
    size_t i = 0;
	size_t current_field = 0;

    for (size_t next = s.find(delimiter, last);
         tokens.size() < fields.size() && next != string::npos;
         next = s.find(delimiter, last), ++i) {
		if (i == fields[current_field]) {
          tokens.push_back(s.substr(last, next-last));
           ++current_field;
		}
        last = next + delim_length;
    }

    return tokens;
}



//template<>
//TaxonomyEntry<uint32_t, uint64_t>::TaxonomyEntry () {
//	read_counts = 0;
//	read_counts_children = 0;
//}
template<typename TAXID, typename READCOUNTS>
bool TaxonomyEntryPtr_comp<TAXID,READCOUNTS>::operator() ( const TaxonomyEntry<TAXID,READCOUNTS>* a, const TaxonomyEntry<TAXID,READCOUNTS>* b) const {
	        return ((reads(a->read_counts)+reads(a->read_counts_children)) > (reads(b->read_counts)+reads(b->read_counts_children)));
			    }


template<typename TAXID, typename READCOUNTS>
std::unordered_map<std::string, TAXID> TaxonomyDB<TAXID,READCOUNTS>::getScientificNameMap() const {
	std::unordered_map<std::string, TAXID> scientificNameMap;
	for (const auto & tax : taxIDsAndEntries) {
		scientificNameMap[tax.second.scientificName] = tax.first;
    }
	return scientificNameMap;
}

template<typename TAXID, typename READCOUNTS>
unordered_map<TAXID, TAXID> TaxonomyDB<TAXID,READCOUNTS>::getParentMap() const {
	unordered_map<TAXID, TAXID> Parent_map;
	for (const auto & tax : taxIDsAndEntries) {
		if (tax.first != 0)
			Parent_map[tax.first] = tax.second.parentTaxonomyID;
    }
    Parent_map[1] = 1;
	return Parent_map;
}

template<typename TAXID, typename READCOUNTS>
void TaxonomyDB<TAXID,READCOUNTS>::createPointers() {
  for (auto& tax : taxIDsAndEntries) {
  if (tax.second.parentTaxonomyID != tax.first) {
    auto parentIt = taxIDsAndEntries.find(tax.second.parentTaxonomyID);
    if (parentIt != taxIDsAndEntries.end()) {
      tax.second.parent = &(parentIt->second);
      parentIt->second.children.push_back(&tax.second);
    }
  }
  }
}

template<typename TAXID, typename READCOUNTS>
TaxonomyDB<TAXID,READCOUNTS>::TaxonomyDB() { }

template<typename TAXID, typename READCOUNTS>
TaxonomyDB<TAXID,READCOUNTS>::TaxonomyDB(const std::string inFileName) {
  log_msg("Building taxonomy index");
  readTaxonomyIndex(inFileName);
  createPointers();
  log_msg("Built a taxonomy tree with " + std::to_string(taxIDsAndEntries.size()) +
      " nodes");
}

template<typename TAXID, typename READCOUNTS>
void TaxonomyDB<TAXID,READCOUNTS>::parseNodesDump(const std::string nodesDumpFileName) {
  std::ifstream nodesDumpFile(nodesDumpFileName);
  if (!nodesDumpFile.is_open())
    throw std::runtime_error("unable to open nodes file");
  std::string line;

  TAXID taxonomyID;
  TAXID parentTaxonomyID;
  std::string rank;

  while (nodesDumpFile.good()) {
    getline(nodesDumpFile, line);
    std::vector<std::string> tokens = tokenise(line, "\t|\t", 3, 2);
    if (tokens.size() < 3) {
	  continue;
	}

	taxonomyID = string_to_T<TAXID>(tokens[0]);
    parentTaxonomyID = string_to_T<TAXID>(tokens[1]);
    rank = tokens[2];

    auto entryIt = taxIDsAndEntries.find(taxonomyID);
	if (entryIt == taxIDsAndEntries.end()) {
	  taxIDsAndEntries[taxonomyID] = TaxonomyEntry<TAXID,READCOUNTS>(taxonomyID, parentTaxonomyID, rank);
	} else {
      entryIt->second.parentTaxonomyID = parentTaxonomyID;
      entryIt->second.rank = rank;
    }
  }
}

template<typename TAXID, typename READCOUNTS>
void TaxonomyDB<TAXID,READCOUNTS>::parseNamesDump(const std::string namesDumpFileName) {
  std::ifstream namesDumpFile(namesDumpFileName);
  if (!namesDumpFile.is_open())
    throw std::runtime_error("unable to open names file");
  std::string line;

  TAXID taxonomyID;
  std::string scientificName;
  while (namesDumpFile.good()) {
    getline(namesDumpFile, line);
    std::vector<std::string> tokens = tokenise(line, "\t|\t", 4, 2);
    if (tokens.size() < 4 || tokens[3] != "scientific name") {
	  continue;
	}
    taxonomyID = string_to_T<TAXID>(tokens[0]);
    scientificName = tokens[1];

    auto entryIt = taxIDsAndEntries.find(taxonomyID);
	if (entryIt == taxIDsAndEntries.end()) {
	  taxIDsAndEntries[taxonomyID] = TaxonomyEntry<TAXID,READCOUNTS>(taxonomyID, scientificName);
	} else {
      entryIt->second.scientificName = scientificName;
    }
  }
}

template<typename TAXID, typename READCOUNTS>
void TaxonomyDB<TAXID,READCOUNTS>::writeTaxonomyIndex(std::ostream & outs,
									const std::string namesDumpFileName,
                                    const std::string nodesDumpFileName) {
  parseNodesDump(nodesDumpFileName);
  parseNamesDump(namesDumpFileName);
  writeTaxonomyIndex(outs);
}

template<typename KeyType, typename ValueType>
std::vector<KeyType> getSortedKeys(const std::unordered_map<KeyType, ValueType>& unordered) {
  std::vector<KeyType> keys;
  keys.reserve (unordered.size());
  for (auto& it : unordered) {
	      keys.push_back(it.first);
  }
  std::sort (keys.begin(), keys.end());
  return keys;
}

template<typename TAXID, typename READCOUNTS>
void TaxonomyDB<TAXID,READCOUNTS>::writeTaxonomyIndex(std::ostream & outs) const {
  for (TAXID& key : getSortedKeys(taxIDsAndEntries)) {
	const auto& entry = taxIDsAndEntries.at(key);
    outs << key << "\t" << entry.parentTaxonomyID << "\t"
            << entry.scientificName << "\t" << entry.rank << "\n";
  }
}



template<typename TAXID, typename READCOUNTS>
void TaxonomyDB<TAXID,READCOUNTS>::readTaxonomyIndex(const std::string inFileName) {
  std::ifstream inFile(inFileName);
  if (!inFile.is_open())
    throw std::runtime_error("unable to open taxonomy index file " + inFileName);

  TAXID taxonomyID, parentTaxonomyID;
  std::string scientificName, rank;

  std::string line;
  while (!inFile.eof()) {
	inFile >> taxonomyID >> parentTaxonomyID;
	inFile.get(); // read tab
	std::getline(inFile, scientificName, '\t');
	std::getline(inFile, rank, '\n');
    TaxonomyEntry<TAXID,READCOUNTS> newEntry(taxonomyID, parentTaxonomyID, rank, scientificName);

	//cerr << "inserting " << taxonomyID << ";" << parentTaxonomyID << ";" << rank << ";" << scientificName << endl;
    taxIDsAndEntries.insert({
      taxonomyID, newEntry
    });
  }
  taxIDsAndEntries.insert({
	0, {0, 0, "no rank", "unclassified" }
  });
}

template<typename TAXID, typename READCOUNTS>
TAXID TaxonomyDB<TAXID,READCOUNTS>::getLowestCommonAncestor(
    const std::vector<TAXID>& taxIDs) const {
  if (taxIDs.size() == 0) {
    return 0;
  }
  std::vector<std::vector<READCOUNTS> > paths;
  for (auto& taxID : taxIDs) {
    bool good = true;
    std::vector<READCOUNTS> path;
    TAXID tempTaxID = taxID;
    while (tempTaxID != 0) {
      path.push_back(tempTaxID);
      tempTaxID = getParentTaxID(tempTaxID);
    }
    if (good) paths.push_back(path);
  }
  if (paths.size() == 0) {
    return 0;
  }
  for (auto& path : paths)
    std::reverse(path.begin(), path.end());
  std::sort(paths.begin(), paths.end(),
            [](std::vector<READCOUNTS> i, std::vector<READCOUNTS> j) {
    return i.size() < j.size();
  });
  TAXID consensus = 0;
  for (unsigned i = 0; i < paths[0].size(); i++) {
    TAXID temp = 0;
    for (auto& path : paths) {
      if (temp == 0)
        temp = path[i];
      else if (temp != path[i]) {
        return consensus;
      }
    }
    consensus = temp;
  }
  return consensus;
}

template<typename TAXID, typename READCOUNTS>
TAXID TaxonomyDB<TAXID,READCOUNTS>::getParentTaxID(const TAXID taxID) const {
  auto entry = taxIDsAndEntries.find(taxID);
  if (entry != taxIDsAndEntries.end() && entry->second.parentTaxonomyID != 1)
    return entry->second.parentTaxonomyID;
  else
    return 0;
}

template<typename TAXID, typename READCOUNTS>
std::string TaxonomyDB<TAXID,READCOUNTS>::getScientificName(const TAXID taxID) const {
  auto entry = taxIDsAndEntries.find(taxID);
  if (entry != taxIDsAndEntries.end()) {
    return entry->second.scientificName;
  } else
    return std::string();
}

template<typename TAXID, typename READCOUNTS>
std::string TaxonomyDB<TAXID,READCOUNTS>::getRank(const TAXID taxID) const {
  auto entry = taxIDsAndEntries.find(taxID);
  if (entry != taxIDsAndEntries.end()) {
    return entry->second.rank;
  } else
    return std::string();
}

template<typename TAXID, typename READCOUNTS>
std::string TaxonomyDB<TAXID,READCOUNTS>::getLineage(TAXID taxonomyID) const {
  std::string lineage;
  while (true) {
    // 131567 = Cellular organisms
    if (taxonomyID != 131567) {
      if (lineage.size()) lineage.insert(0, "; ");
      lineage.insert(0, getScientificName(taxonomyID));
      if (getRank(taxonomyID) == "species") lineage.clear();
    }
    taxonomyID = getParentTaxID(taxonomyID);
    if (taxonomyID == 0) {
      if (lineage.size()) lineage.append(".");
      break;
    }
  }
  return lineage;
}

template<typename TAXID, typename READCOUNTS>
std::string TaxonomyDB<TAXID,READCOUNTS>::getMetaPhlAnLineage(TAXID taxonomyID) const {
  std::string rank = getRank(taxonomyID);
  if (rank == "superphylum") return std::string();
  std::string lineage;
  while (true) {
    // 131567 = Cellular organisms
    if (taxonomyID != 131567) {
      std::string rank = getRank(taxonomyID);
      if (rank == "species") {
        lineage.insert(0, "|s__");
        lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "genus") {
        lineage.insert(0, "|g__");
        lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "family") {
        lineage.insert(0, "|f__");
        lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "order") {
        lineage.insert(0, "|o__");
        lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "class") {
        lineage.insert(0, "|c__");
        lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "phylum") {
        lineage.insert(0, "|p__");
        lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "superkingdom") {
        lineage.insert(0, "k__");
        lineage.insert(3, getScientificName(taxonomyID));
      }
    }
    taxonomyID = getParentTaxID(taxonomyID);
    if (taxonomyID == 0) {
      break;
    }
  }
  std::replace(lineage.begin(), lineage.end(), ' ', '_');
  return lineage;
}

template<typename TAXID, typename READCOUNTS>
TAXID TaxonomyDB<TAXID,READCOUNTS>::getTaxIDAtRank(const TAXID taxID,
                                    const std::string& rank) const {
  auto entry = taxIDsAndEntries.find(taxID);
  while (entry != taxIDsAndEntries.end() &&
         entry->second.parentTaxonomyID != 1) {
    if (entry->second.rank == rank) {
      return entry->second.taxonomyID;
    } else
      entry = taxIDsAndEntries.find(entry->second.parentTaxonomyID);
  }
  return 0;
}

template<typename TAXID, typename READCOUNTS>
int TaxonomyDB<TAXID,READCOUNTS>::isBelowInTree(TAXID upper, TAXID lower) const {
  auto entry = taxIDsAndEntries.find(lower);
  unsigned level = 0;
  while (entry != taxIDsAndEntries.end() &&
         entry->second.parentTaxonomyID != 1) {
    if (entry->first == upper) {
      return level;
    } else {
      entry = taxIDsAndEntries.find(entry->second.parentTaxonomyID);
      level++;
    }
  }
  return -1;
}

template<typename TAXID, typename READCOUNTS>
bool TaxonomyDB<TAXID,READCOUNTS>::isSubSpecies(TAXID taxonomyID) const {
  bool isSubSpecies = false;
  auto entry = taxIDsAndEntries.find(taxonomyID);
  int numLevels = 0;
  while (entry != taxIDsAndEntries.end() &&
         entry->second.parentTaxonomyID != 1) {
    if (entry->second.rank == "species") {
      if (numLevels > 0) {
        isSubSpecies = true;
      }
      break;
    } else
      entry = taxIDsAndEntries.find(entry->second.parentTaxonomyID);
    numLevels++;
  }
  return isSubSpecies;
}

template<typename TAXID, typename READCOUNTS>
void TaxonomyDB<TAXID,READCOUNTS>::addCounts(const TAXID taxid, const READCOUNTS& read_counts_) {
	auto it = taxIDsAndEntries.find(taxid);
		if (it == taxIDsAndEntries.end()) {
			cerr << "No taxonomy entry for " << taxid << "!!" << endl;
			return;
		}
		TaxonomyEntry<TAXID,READCOUNTS>* tax = &it->second;
		//cerr << taxid << " rc before: " << tax->read_counts << endl;
		tax->read_counts += read_counts_;
		//cerr << taxid << " rc after:  " << tax->read_counts << endl;

		while (tax->parent != nullptr) {
			tax = tax->parent;
			tax->read_counts_children += read_counts_;
		}
}

template<typename TAXID, typename READCOUNTS>
void TaxonomyDB<TAXID,READCOUNTS>::fillCounts(const unordered_map<TAXID, READCOUNTS>& taxon_counts) {
	for (auto& elem : taxon_counts) {
		addCounts(elem.first, elem.second);
	 }

	for (auto& tax : taxIDsAndEntries) {
		std::sort(tax.second.children.begin(), tax.second.children.end(),TaxonomyEntryPtr_comp<TAXID,READCOUNTS>());
	}
}


template<typename TAXID, typename READCOUNTS>
TaxReport<TAXID,READCOUNTS>::TaxReport(std::ostream& reportOfb, TaxonomyDB<TAXID,READCOUNTS>& taxdb, bool show_zeros) : _reportOfb(reportOfb), _taxdb(taxdb), _show_zeros(show_zeros) {
	_report_cols = {REPORTCOLS::PERCENTAGE, REPORTCOLS::NUM_READS_CLADE, REPORTCOLS::NUM_READS, REPORTCOLS::NUM_UNIQUE_KMERS, REPORTCOLS::NUM_KMERS, REPORTCOLS::TAX_RANK, REPORTCOLS::TAX_ID, REPORTCOLS::SPACED_NAME};
}

template<typename TAXID, typename READCOUNTS>
void TaxReport<TAXID,READCOUNTS>::printReport(std::string format, std::string rank) {
	_total_n_reads =
			reads(_taxdb.taxIDsAndEntries.at(0).read_counts) +
			reads(_taxdb.taxIDsAndEntries.at(0).read_counts_children) +
			reads(_taxdb.taxIDsAndEntries.at(1).read_counts) +
			reads(_taxdb.taxIDsAndEntries.at(1).read_counts_children);// +
	if (_total_n_reads == 0) {
		std::cerr << "total number of reads is zero - not creating a report!" << endl;
		return;
	}

	if (format == "kraken") {
		// A: print number of unidentified reads
		printReport(_taxdb.taxIDsAndEntries.at(0),0u);
		// B: print normal results
		printReport(_taxdb.taxIDsAndEntries.at(1),0u);
		// C: Print Unclassified stuff
		//printReport(_taxdb.taxIDsAndEntries.at(-1),0u);
	} else {
		// print stuff at a certain level ..
		//_uid_abundance;
		//_taxinfo

	}
}

template<typename TAXID, typename READCOUNTS>
void TaxReport<TAXID,READCOUNTS>::printReport(TaxonomyEntry<TAXID,READCOUNTS>& tax, unsigned depth) {
	if (_show_zeros || (reads(tax.read_counts)+reads(tax.read_counts_children)) > 0) {
		printLine(tax, depth);
		for (auto child : tax.children)
			printReport(*child, depth+1);
	}
}

template<typename TAXID, typename READCOUNTS>
void TaxReport<TAXID,READCOUNTS>::printLine(TaxonomyEntry<TAXID,READCOUNTS>& tax, unsigned depth) {
	for (auto& col : _report_cols) {
		switch (col) {
		case REPORTCOLS::NAME:        _reportOfb << tax.scientificName ; break;
		case REPORTCOLS::SPACED_NAME:       _reportOfb << string(2*depth, ' ') + tax.scientificName; break;
		case REPORTCOLS::TAX_ID:     _reportOfb << (tax.taxonomyID == (uint32_t)-1? -1 : (int32_t) tax.taxonomyID); break;
		case REPORTCOLS::DEPTH:     _reportOfb << depth; break;
		case REPORTCOLS::PERCENTAGE:  _reportOfb << 100.0*(reads(tax.read_counts) + reads(tax.read_counts_children))/_total_n_reads; break;
		//case REPORTCOLS::ABUNDANCE:  _reportOfb << 100*counts.abundance[0]; break;
		//case REPORTCOLS::ABUNDANCE_LEN:  _reportOfb << 100*counts.abundance[1]; break;
		case REPORTCOLS::NUM_READS_CLADE:  _reportOfb << (reads(tax.read_counts) + reads(tax.read_counts_children)); break;
		case REPORTCOLS::NUM_READS:  _reportOfb << reads(tax.read_counts); break;
		//case REPORTCOLS::NUM_UNIQUE_KMERS: _reportOfb << tax.kmers.cardinality(); break;
		//case REPORTCOLS::NUM_KMERS: _reportOfb << tax.numKmers; break;
		//case REPORTCOLS::GENOME_SIZE: ; break;
		//case REPORTCOLS::NUM_WEIGHTED_READS: ; break;
		//case REPORTCOLS::SUM_SCORE: ; break;
		case REPORTCOLS::TAX_RANK: _reportOfb << tax.rank; break;
		default: _reportOfb << "NA";
		}
		if (&col == &_report_cols.back()) {
			_reportOfb << '\n';
		} else {
			_reportOfb << '\t';
		}
	}
}


  // Return lowest common ancestor of a and b
  // LCA(0,x) = LCA(x,0) = x
  // Default ancestor is 1 (root of tree)
uint32_t lca(unordered_map<uint32_t, uint32_t> &parent_map, uint32_t a, uint32_t b)
  {
    if (a == 0 || b == 0)
      return a ? a : b;

    // create a path from a to the root
	std::unordered_set<uint32_t> a_path;
    while (a > 0 && a != parent_map[a]) {
	  if (a == b)
		  return a;
      a_path.insert(a);
      a = parent_map[a];
    }

    // search for b in the path from a to the root
    while (b > 0 && b != parent_map[b]) {
      if (a_path.count(b) > 0)
        return b;
      b = parent_map[b];
    }
    return 1;
  }

template<typename K,typename V>
inline
V find_or_use_default(const std::unordered_map<K, V>& my_map, const K& query, const V default_value) {
	auto itr = my_map.find(query);

	if (itr == my_map.end()) {
		return default_value;
	}

	return itr->second;
}



#endif /* TAXD_DB_H_ */
