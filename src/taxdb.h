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
#include "hyperloglogplus.h"
#include "report-cols.h"

struct ReadCounts {
	uint64_t n_reads = 0;
	uint64_t n_kmers = 0;
    HyperLogLogPlusMinus<uint64_t> kmers; // unique k-mer count per taxon
	void add_kmer(uint64_t kmer) {
		++ n_kmers;
		kmers.add(kmer);
	}
    ReadCounts& operator+=(const ReadCounts& b) {
        n_reads += b.n_reads;
        n_kmers += b.n_kmers;
		kmers += kmers;
        return *this;
    }
};


void log (const std::string& s) {
	std::cerr << s << "\n";
}

template<typename T>
uint64_t string_to_T(string str) {
  stringstream stream(str);
  T result;
  stream >> result;
  return result;
}

std::vector<std::string> tokenise(const std::string &s, const std::string& delimiter, size_t max_fields, size_t end_chars) {
    std::vector<std::string> tokens(max_fields);
    size_t delim_length = delimiter.length();
    size_t last = 0;
    size_t i = 0;

    for (size_t next = s.find(delimiter, last);
         i < max_fields && next != string::npos;
         next = s.find(delimiter, last), ++i) {
        tokens[i] = s.substr(last, next-last);
        last = next + delim_length;
    }
    if (i < max_fields) {
        tokens[max_fields-1] = s.substr(last, s.length()-last-end_chars);
    }

    return tokens;
}

template<typename TAXID>
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

  inline bool operator==(const TaxonomyEntry& other) const {
    return this->taxonomyID == other.taxonomyID &&
           this->parentTaxonomyID == other.parentTaxonomyID &&
           this->scientificName == other.scientificName;
  }
  TaxonomyEntry* parent = nullptr;
  std::vector<TaxonomyEntry*> children;

  unsigned numReadsAligned = 0;
  unsigned numReadsAlignedToChildren = 0;
  bool used = false;
  uint64_t genomeSize = 0;
  uint64_t genomeSizeOfChildren = 0;
  uint64_t numBelow = 0;
  uint64_t numKmers = 0;
  HyperLogLogPlusMinus<uint64_t> kmers;
};


template<typename TAXID>
struct TaxonomyEntryPtr_comp {
	bool operator() ( const TaxonomyEntry<TAXID>* a, const TaxonomyEntry<TAXID>* b) const { 
		return ((a->numReadsAligned+a->numReadsAlignedToChildren) > (b->numReadsAligned+b->numReadsAlignedToChildren)); 
	}
};

template<typename TAXID>
class TaxonomyDB {
 public:
  TaxonomyDB(const std::string inFileName);
  TaxonomyDB() {};
  //std::unordered_map<std::string, TAXID> seqIDsAndTaxIds;
  std::unordered_map<TAXID, TaxonomyEntry<TAXID> > taxIDsAndEntries;
  void parseNamesDump(const std::string namesDumpFileName);
  void parseNodesDump(const std::string nodesDumpFileName);
  TAXID getTaxIDAtRank(const TAXID taxID, const std::string& rank) const;
  std::string getScientificName(const TAXID taxID) const;
  std::string getRank(const TAXID taxID) const;
  TAXID getLowestCommonAncestor(const std::vector<TAXID>& taxIDs) const;
  TAXID getParentTaxID(const TAXID taxID) const;
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
  void fillCounts(const unordered_map<TAXID, ReadCounts>& taxon_counts);
  void createPointers();
  void printReport();
};

template<typename TAXID>
void TaxonomyDB<TAXID>::createPointers() {
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

template<typename TAXID>
TaxonomyDB<TAXID>::TaxonomyDB(const std::string inFileName) {
  log("Building taxonomy index");
  readTaxonomyIndex(inFileName);
  createPointers();
  log("Built a taxonomy tree with " + std::to_string(taxIDsAndEntries.size()) +
      " nodes");
}

template<typename TAXID>
void TaxonomyDB<TAXID>::parseNodesDump(const std::string nodesDumpFileName) {
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
	  taxIDsAndEntries[taxonomyID] = TaxonomyEntry<TAXID>(taxonomyID, parentTaxonomyID, rank);
	} else {
      entryIt->second.parentTaxonomyID = parentTaxonomyID;
      entryIt->second.rank = rank;
    }
  }
}

template<typename TAXID>
void TaxonomyDB<TAXID>::parseNamesDump(const std::string namesDumpFileName) {
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
	  taxIDsAndEntries[taxonomyID] = TaxonomyEntry<TAXID>(taxonomyID, scientificName);
	} else {
      entryIt->second.scientificName = scientificName;
    }
  }
}

template<typename TAXID>
void TaxonomyDB<TAXID>::writeTaxonomyIndex(std::ostream & outs,
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

template<typename TAXID>
void TaxonomyDB<TAXID>::writeTaxonomyIndex(std::ostream & outs) const {
  for (TAXID& key : getSortedKeys(taxIDsAndEntries)) {
	const auto& entry = taxIDsAndEntries.at(key);
    outs << key << "\t" << entry.parentTaxonomyID << "\t"
            << entry.scientificName << "\t" << entry.rank << "\n";
  }
}



template<typename TAXID>
void TaxonomyDB<TAXID>::readTaxonomyIndex(const std::string inFileName) {
  std::ifstream inFile(inFileName);
  if (!inFile.is_open())
    throw std::runtime_error("unable to open taxonomy index file");

  TAXID taxonomyID, parentTaxonomyID;
  std::string scientificName, rank;

  std::string line;
  while (!inFile.eof()) {
	inFile >> taxonomyID >> parentTaxonomyID;
	inFile.get(); // read tab
	std::getline(inFile, scientificName, '\t');
	std::getline(inFile, rank, '\n');
    TaxonomyEntry<TAXID> newEntry(taxonomyID, parentTaxonomyID, rank, scientificName);

	//cerr << "inserting " << taxonomyID << ";" << parentTaxonomyID << ";" << rank << ";" << scientificName << endl;
    taxIDsAndEntries.insert({
      taxonomyID, newEntry
    });
  }
  taxIDsAndEntries.insert({
	0, {0, 0, "no rank", "unclassified" }
  });
}

template<typename TAXID>
TAXID TaxonomyDB<TAXID>::getLowestCommonAncestor(
    const std::vector<TAXID>& taxIDs) const {
  if (taxIDs.size() == 0) {
    return 0;
  }
  std::vector<std::vector<TAXID> > paths;
  for (auto& taxID : taxIDs) {
    bool good = true;
    std::vector<TAXID> path;
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
            [](std::vector<TAXID> i, std::vector<TAXID> j) {
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

template<typename TAXID>
TAXID TaxonomyDB<TAXID>::getParentTaxID(const TAXID taxID) const {
  auto entry = taxIDsAndEntries.find(taxID);
  if (entry != taxIDsAndEntries.end() && entry->second.parentTaxonomyID != 1)
    return entry->second.parentTaxonomyID;
  else
    return 0;
}

template<typename TAXID>
std::string TaxonomyDB<TAXID>::getScientificName(const TAXID taxID) const {
  auto entry = taxIDsAndEntries.find(taxID);
  if (entry != taxIDsAndEntries.end()) {
    return entry->second.scientificName;
  } else
    return std::string();
}

template<typename TAXID>
std::string TaxonomyDB<TAXID>::getRank(const TAXID taxID) const {
  auto entry = taxIDsAndEntries.find(taxID);
  if (entry != taxIDsAndEntries.end()) {
    return entry->second.rank;
  } else
    return std::string();
}

template<typename TAXID>
std::string TaxonomyDB<TAXID>::getLineage(TAXID taxonomyID) const {
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

template<typename TAXID>
std::string TaxonomyDB<TAXID>::getMetaPhlAnLineage(TAXID taxonomyID) const {
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

template<typename TAXID>
TAXID TaxonomyDB<TAXID>::getTaxIDAtRank(const TAXID taxID,
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

template<typename TAXID>
int TaxonomyDB<TAXID>::isBelowInTree(TAXID upper, TAXID lower) const {
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

template<typename TAXID>
bool TaxonomyDB<TAXID>::isSubSpecies(TAXID taxonomyID) const {
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

template<typename TAXID>
void TaxonomyDB<TAXID>::fillCounts(const unordered_map<TAXID, ReadCounts>& taxon_counts) {
	for (auto& elem : taxon_counts) {
		auto it = taxIDsAndEntries.find(elem.first);
		if (it == taxIDsAndEntries.end()) {
			cerr << "No taxonomy entry for " << elem.first << "!!" << endl;
			continue;
		}
		TaxonomyEntry<TAXID>* tax = &it->second;
		//cerr << "fill done: "<< elem.first << endl;
		tax->numReadsAligned += elem.second.n_reads;
		tax->numKmers += elem.second.n_kmers;
		tax->kmers += elem.second.kmers;

		//std::cerr << "adding " << elem.second.n_reads << " to " << tax->scientificName << ": ";

		while (tax->parent != nullptr) {
			tax = tax->parent;
			//std::cerr << " >> " << tax->scientificName;
			tax->numReadsAlignedToChildren += elem.second.n_reads;
			tax->numKmers += elem.second.n_kmers;
			tax->kmers += elem.second.kmers;
		}
		//std::cerr << endl;
	 }

	for (auto& tax : taxIDsAndEntries) {
		std::sort(tax.second.children.begin(), tax.second.children.end(),TaxonomyEntryPtr_comp<TAXID>());
	}
}


template<typename TAXID>
class TaxReport {
private:
	std::ostream& _reportOfb;
	TaxonomyDB<TAXID> & _taxdb;
	std::vector<REPORTCOLS> _report_cols;
	uint64_t _total_n_reads;
	bool _show_zeros;

	void printLine(TaxonomyEntry<TAXID>& tax, unsigned depth);

public:
	TaxReport(std::ostream& _reportOfb, TaxonomyDB<TAXID> & taxdb, bool _show_zeros);

	void printReport(std::string format, std::string rank);
	void printReport(TaxonomyEntry<TAXID>& tax, unsigned depth);
};

template<typename TAXID>
TaxReport<TAXID>::TaxReport(std::ostream& reportOfb, TaxonomyDB<TAXID>& taxdb, bool show_zeros) : _reportOfb(reportOfb), _taxdb(taxdb), _show_zeros(show_zeros) {
	_report_cols = {REPORTCOLS::PERCENTAGE, REPORTCOLS::NUM_READS_CLADE, REPORTCOLS::NUM_READS, REPORTCOLS::NUM_UNIQUE_KMERS, REPORTCOLS::NUM_KMERS, REPORTCOLS::TAX_RANK, REPORTCOLS::TAX_ID, REPORTCOLS::SPACED_NAME};
}

template<typename TAXID>
void TaxReport<TAXID>::printReport(std::string format, std::string rank) {
	_total_n_reads =
			_taxdb.taxIDsAndEntries.at(0).numReadsAligned +
			_taxdb.taxIDsAndEntries.at(0).numReadsAlignedToChildren +
			_taxdb.taxIDsAndEntries.at(1).numReadsAligned +
			_taxdb.taxIDsAndEntries.at(1).numReadsAlignedToChildren;// +
			//_taxdb.taxIDsAndEntries.at(-1).numReadsAligned +
			//_taxdb.taxIDsAndEntries.at(-1).numReadsAlignedToChildren; // -1 is a magic number in centrifuge for reads not matched to the taxonomy tree
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

template<typename TAXID>
void TaxReport<TAXID>::printReport(TaxonomyEntry<TAXID>& tax, unsigned depth) {

	if (_show_zeros || (tax.numReadsAligned+tax.numReadsAlignedToChildren) > 0) {
		printLine(tax, depth);

		for (auto child : tax.children) {
			printReport(*child, depth+1);
		}
	}

}

template<typename TAXID>
void TaxReport<TAXID>::printLine(TaxonomyEntry<TAXID>& tax, unsigned depth) {
	for (auto& col : _report_cols) {
		switch (col) {
		case REPORTCOLS::NAME:        _reportOfb << tax.scientificName ; break;
		case REPORTCOLS::SPACED_NAME:       _reportOfb << string(2*depth, ' ') + tax.scientificName; break;
		case REPORTCOLS::TAX_ID:     _reportOfb << (tax.taxonomyID == (uint32_t)-1? -1 : (int32_t) tax.taxonomyID); break;
		case REPORTCOLS::DEPTH:     _reportOfb << depth; break;
		case REPORTCOLS::PERCENTAGE:  _reportOfb << 100.0*(tax.numReadsAligned + tax.numReadsAlignedToChildren)/_total_n_reads; break;
		//case REPORTCOLS::ABUNDANCE:  _reportOfb << 100*counts.abundance[0]; break;
		//case REPORTCOLS::ABUNDANCE_LEN:  _reportOfb << 100*counts.abundance[1]; break;
		case REPORTCOLS::NUM_READS_CLADE:  _reportOfb << (tax.numReadsAligned + tax.numReadsAlignedToChildren); break;
		case REPORTCOLS::NUM_READS:  _reportOfb << tax.numReadsAligned; break;
		case REPORTCOLS::NUM_UNIQUE_KMERS: _reportOfb << tax.kmers.cardinality(); break;
		case REPORTCOLS::NUM_KMERS: _reportOfb << tax.numKmers; break;
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




#endif /* TAXD_DB_H_ */
