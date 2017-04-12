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

typedef uint32_t TaxId;

struct ReadCounts {
	uint64_t n_reads = 0;
	uint64_t n_kmers = 0;
    HyperLogLogPlusMinus<uint64_t> kmers; // unique k-mer count per taxon
};


void log (const std::string& s) {
	std::cerr << s << "\n";
}

std::vector<std::string> tokenise(const std::string &line, const std::string& delimiters) {
	std::vector<std::string> tokens;
	// Skip delimiters at beginning.
	std::string::size_type lastPos = line.find_first_not_of(delimiters, 0);
	std::string::size_type pos = line.find_first_of(delimiters, lastPos);
	while (std::string::npos != pos || std::string::npos != lastPos) {
	  tokens.push_back(line.substr(lastPos, pos - lastPos));
	  // Skip delimiters.  Note the "not_of"
	  lastPos = line.find_first_not_of(delimiters, pos);
	  pos = line.find_first_of(delimiters, lastPos);
	}
	return tokens;
}

class TaxonomyEntry {
 public:
  uint32_t taxonomyID = 0;
  uint32_t parentTaxonomyID = 0;
  std::string rank;
  std::string scientificName;

  TaxonomyEntry() {}
  TaxonomyEntry(uint32_t taxonomyID_, uint32_t parentTaxonomyID_, std::string rank_, std::string scientificName_) :
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

struct TaxonomyEntryPtr_comp {
	bool operator() ( const TaxonomyEntry* a, const TaxonomyEntry* b) const { 
		return ((a->numReadsAligned+a->numReadsAlignedToChildren) > (b->numReadsAligned+b->numReadsAlignedToChildren)); 
	}
};

class TaxonomyDB {
 public:
  TaxonomyDB(const std::string inFileName);
  TaxonomyDB() {};
  std::unordered_map<uint32_t, TaxonomyEntry> taxIDsAndEntries;
  void parseNamesDump(const std::string namesDumpFileName);
  void parseNodesDump(const std::string nodesDumpFileName);
  uint32_t getTaxIDAtRank(const uint32_t taxID, const std::string& rank) const;
  std::string getScientificName(const uint32_t taxID) const;
  std::string getRank(const uint32_t taxID) const;
  uint32_t getLowestCommonAncestor(const std::vector<uint32_t>& taxIDs) const;
  uint32_t getParentTaxID(const uint32_t taxID) const;
  std::string getLineage(uint32_t taxonomyID) const;
  std::string getMetaPhlAnLineage(uint32_t taxonomyID) const;
  char* getIndexFileName(const uint32_t hostTaxID) const;
  void readTaxonomyIndex(const std::string inFileName);
  void writeTaxonomyIndex(std::ostream & outs,
                          const std::string namesDumpFileName,
                          const std::string nodesDumpFileName);
  bool isSubSpecies(uint32_t taxonomyID) const;
  int isBelowInTree(uint32_t upper, uint32_t lower) const;
  void fillCounts(const unordered_map<uint32_t, ReadCounts>& taxon_counts);
  void createPointers();
  void printReport();
};


void TaxonomyDB::createPointers() {
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
TaxonomyDB::TaxonomyDB(const std::string inFileName) {
  log("Building taxonomy index");
  readTaxonomyIndex(inFileName);
  createPointers();
  log("Built a taxonomy tree with " + std::to_string(taxIDsAndEntries.size()) +
      " nodes");
}

void TaxonomyDB::parseNodesDump(const std::string nodesDumpFileName) {
  std::ifstream nodesDumpFile(nodesDumpFileName);
  if (!nodesDumpFile.is_open())
    throw std::runtime_error("unable to open nodes file");
  std::string line;
  while (nodesDumpFile.good()) {
    getline(nodesDumpFile, line);
    std::vector<std::string> tokens = tokenise(line, "\t|");
    if (tokens.size() > 2) {
      TaxonomyEntry newEntry;
      newEntry.taxonomyID = stoi(tokens[0]);
      newEntry.parentTaxonomyID = stoi(tokens[1]);
      newEntry.rank = tokens[2];
      auto entryIt = taxIDsAndEntries.insert({
        newEntry.taxonomyID, newEntry
      });
      if (!entryIt.second) {
        entryIt.first->second.taxonomyID = newEntry.taxonomyID;
        newEntry.parentTaxonomyID = stoi(tokens[1]);
      }
    }
  }
}

void TaxonomyDB::parseNamesDump(const std::string namesDumpFileName) {
  std::ifstream namesDumpFile(namesDumpFileName);
  if (!namesDumpFile.is_open())
    throw std::runtime_error("unable to open names file");
  std::string line;
  while (namesDumpFile.good()) {
    getline(namesDumpFile, line);
    std::vector<std::string> tokens = tokenise(line, "|");
    for (auto& token : tokens) {
      if (token.size() > 1) {
        if (token[0] == '\t') token.erase(0, 1);
        if (token[token.size() - 1] == '\t') token.erase(token.size() - 1, 1);
      }
    }
    if (tokens.size() > 3) {
      TaxonomyEntry newEntry;
      newEntry.taxonomyID = stoi(tokens[0]);
      //            for(auto & token : tokens)
      //                std::cout<<token<<"\n";
      if (tokens[3] == "scientific name") {
        //      std::cout<<"Found\n";
        newEntry.scientificName = tokens[1];
        //      std::cout<<newEntry.scientificName<<"\n";
      } else
        continue;
      auto entryIt = taxIDsAndEntries.insert({
        newEntry.taxonomyID, newEntry
      });
      if (!entryIt.second) {
        entryIt.first->second.scientificName = newEntry.scientificName;
      }
    }
  }
}

void TaxonomyDB::writeTaxonomyIndex(std::ostream & outs,
									const std::string namesDumpFileName,
                                    const std::string nodesDumpFileName) {
  parseNodesDump(nodesDumpFileName);
  parseNamesDump(namesDumpFileName);
  for (auto& entry : taxIDsAndEntries) {
    outs << entry.first << "\t" << entry.second.parentTaxonomyID << "\t"
            << entry.second.scientificName << "\t" << entry.second.rank << "\n";
  }
}

void TaxonomyDB::readTaxonomyIndex(const std::string inFileName) {
  std::ifstream inFile(inFileName);
  if (!inFile.is_open())
    throw std::runtime_error("unable to open taxonomy index file");

  uint32_t taxonomyID, parentTaxonomyID;
  std::string scientificName, rank;

  std::string line;
  while (!inFile.eof()) {
	inFile >> taxonomyID >> parentTaxonomyID;
	inFile.get(); // read tab
	std::getline(inFile, scientificName, '\t');
	std::getline(inFile, rank, '\n');
    TaxonomyEntry newEntry(taxonomyID, parentTaxonomyID, rank, scientificName);

	//cerr << "inserting " << taxonomyID << ";" << parentTaxonomyID << ";" << rank << ";" << scientificName << endl;
    taxIDsAndEntries.insert({
      taxonomyID, newEntry
    });
  }
  taxIDsAndEntries.insert({
	0, {0, 0, "no rank", "unclassified" }
  });
}

uint32_t TaxonomyDB::getLowestCommonAncestor(
    const std::vector<uint32_t>& taxIDs) const {
  if (taxIDs.size() == 0) {
    return 0;
  }
  std::vector<std::vector<uint32_t> > paths;
  for (auto& taxID : taxIDs) {
    bool good = true;
    std::vector<uint32_t> path;
    uint32_t tempTaxID = taxID;
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
            [](std::vector<uint32_t> i, std::vector<uint32_t> j) {
    return i.size() < j.size();
  });
  uint32_t consensus = 0;
  for (unsigned i = 0; i < paths[0].size(); i++) {
    uint32_t temp = 0;
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

uint32_t TaxonomyDB::getParentTaxID(const uint32_t taxID) const {
  auto entry = taxIDsAndEntries.find(taxID);
  if (entry != taxIDsAndEntries.end() && entry->second.parentTaxonomyID != 1)
    return entry->second.parentTaxonomyID;
  else
    return 0;
}

std::string TaxonomyDB::getScientificName(const uint32_t taxID) const {
  auto entry = taxIDsAndEntries.find(taxID);
  if (entry != taxIDsAndEntries.end()) {
    return entry->second.scientificName;
  } else
    return std::string();
}

std::string TaxonomyDB::getRank(const uint32_t taxID) const {
  auto entry = taxIDsAndEntries.find(taxID);
  if (entry != taxIDsAndEntries.end()) {
    return entry->second.rank;
  } else
    return std::string();
}

std::string TaxonomyDB::getLineage(uint32_t taxonomyID) const {
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
std::string TaxonomyDB::getMetaPhlAnLineage(uint32_t taxonomyID) const {
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

uint32_t TaxonomyDB::getTaxIDAtRank(const uint32_t taxID,
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
int TaxonomyDB::isBelowInTree(uint32_t upper, uint32_t lower) const {
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
bool TaxonomyDB::isSubSpecies(uint32_t taxonomyID) const {
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

void TaxonomyDB::fillCounts(const unordered_map<uint32_t, ReadCounts>& taxon_counts) {
	for (auto& elem : taxon_counts) {
		auto it = taxIDsAndEntries.find(elem.first);
		if (it == taxIDsAndEntries.end()) {
			cerr << "No taxonomy entry for " << elem.first << "!!" << endl;
			continue;
		}
		TaxonomyEntry* tax = &it->second;
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
		std::sort(tax.second.children.begin(), tax.second.children.end(),TaxonomyEntryPtr_comp());
	}
}


class TaxReport {
private:
	std::ostream& _reportOfb;
	TaxonomyDB & _taxdb;
	std::vector<REPORTCOLS> _report_cols;
	uint64_t _total_n_reads;
	bool _show_zeros;

	void printLine(TaxonomyEntry& tax, unsigned depth);

public:
	TaxReport(std::ostream& _reportOfb, TaxonomyDB & taxdb, bool _show_zeros);

	void printReport(std::string format, std::string rank);
	void printReport(TaxonomyEntry& tax, unsigned depth);
};

TaxReport::TaxReport(std::ostream& reportOfb, TaxonomyDB& taxdb, bool show_zeros) : _reportOfb(reportOfb), _taxdb(taxdb), _show_zeros(show_zeros) {
	_report_cols = {REPORTCOLS::PERCENTAGE, REPORTCOLS::NUM_READS_CLADE, REPORTCOLS::NUM_READS, REPORTCOLS::NUM_UNIQUE_KMERS, REPORTCOLS::NUM_KMERS, REPORTCOLS::TAX_RANK, REPORTCOLS::TAX_ID, REPORTCOLS::SPACED_NAME};
}

void TaxReport::printReport(std::string format, std::string rank) {
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

void TaxReport::printReport(TaxonomyEntry& tax, unsigned depth) {

	if (_show_zeros || (tax.numReadsAligned+tax.numReadsAlignedToChildren) > 0) {
		printLine(tax, depth);

		for (auto child : tax.children) {
			printReport(*child, depth+1);
		}
	}

}

void TaxReport::printLine(TaxonomyEntry& tax, unsigned depth) {
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
