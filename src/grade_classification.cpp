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

#include "taxdb.h"
#include "quickfile.hpp"
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>
#include <limits>

using namespace std;

//using TAXID = uint32_t;
typedef uint32_t TAXID;

unordered_map<string, uint32_t> read_seqid_mapping(string filename) {
  unordered_map<string, uint32_t> ID_to_taxon_map;
  ifstream map_file(filename.c_str());
  if (map_file.rdstate() & ifstream::failbit) {
    err(EX_NOINPUT, "can't open %s", filename.c_str());
  }
  string line, seq_id;
  uint32_t taxid;

  while (map_file >> seq_id >> taxid) {
    ID_to_taxon_map[seq_id] = taxid;
    map_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  map_file.close();
  cerr << "Read " << ID_to_taxon_map.size() << " taxa mappings" << endl;
  return ID_to_taxon_map;
}

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "Usage: grade_classification taxDB seqid2taxid.map classification_file result_file\n";
    return 1;
  }
  TaxonomyDB<uint32_t> taxdb (argv[1], false);
  unordered_map<string, uint32_t> seqid_map = read_seqid_mapping(argv[2]);
  
  ofstream out_file(argv[4]);
  unordered_set<string> all_ranks;
  unordered_map< string, size_t > rank_counts;
  unordered_map< int, set<TAXID> > simulated_taxids_at_rank;
  unordered_map< int, set<TAXID> > identified_taxids_at_rank;
  unordered_map< int, size_t > correct_reads_at_rank;
  unordered_map< int, size_t > incorrect_reads_at_rank;
  unordered_map< int, size_t > reads_at_higher_rank;
  unordered_set<uint32_t> ignored_taxa;
  size_t total_reads = 0;
  size_t unidentified_reads = 0;
  

  vector<TaxRank::RANK> ranks_of_interest = {TaxRank::RANK::assembly, TaxRank::RANK::species, TaxRank::RANK::genus, TaxRank::RANK::family, TaxRank::RANK::order}; 

  ifstream k_file(argv[3]);
  if (k_file.rdstate() & ifstream::failbit) {
    err(EX_NOINPUT, "can't open %s", argv[3]);
  }

  string line, classification_state, read_id, seq_id;
  uint32_t identified_taxid;
  uint32_t seq_taxid;

  while (k_file.good()) {
    getline(k_file, line);
    if (line.empty())
      continue;
    istringstream iss(line);
    string l;
    string classi;
    iss >> classification_state >> read_id >> identified_taxid >> l;
    iss.get();
    getline(iss,classi);

    ++total_reads;
    if (identified_taxid == 0) {
      ++unidentified_reads;
    }

    // sequence id is after the 5th underscore with random_reads.sh - find it
    size_t pos = 0;
    size_t count = 0;
    do {
      pos = read_id.find("_", pos) + 1;
      ++count;
    } while (count <= 5 && pos != std::string::npos);

    seq_id = read_id.substr(pos);
    auto it = seqid_map.find(seq_id);
    if (it == seqid_map.end()) {
      cerr << "ERROR: Couldn't find taxid for " << seq_id << endl;
      exit(1);
    } else {
      seq_taxid = it->second;
      if (!taxdb.hasTaxon(seq_taxid)) {
        if (ignored_taxa.count(seq_taxid) == 0) {
          cerr << "Ignoring taxon " << seq_taxid << " - not in database" << endl;
          ignored_taxa.insert(seq_taxid);
        }
        continue;
      }
      //cerr <<"seqid" << seq_taxid;

      // go up to species level or next proper (i.e. not 'no rank') rank for
      //  both real and assigned taxon
      if (0) {
        seq_taxid = taxdb.getTaxIDAtRank(seq_taxid, "species");
        uint32_t identified_species_taxid = taxdb.getTaxIDAtRank(identified_taxid, "species");
        if (identified_species_taxid != 0) {
          identified_taxid = identified_species_taxid;
        } else {
          identified_taxid = taxdb.getTaxIDAtNextProperRank(identified_taxid);
        }  
      }
      
      string seq_name = taxdb.getScientificName(seq_taxid);
      // getLowestCommonAncestor returns lca taxon as well as distance between the taxa
      pair<uint32_t, int> lca_taxid_dist = taxdb.getLowestCommonAncestor(seq_taxid, identified_taxid);
      string lca_rank_string = taxdb.getNextProperRank(lca_taxid_dist.first);
      // TaxRank::RANK lca_rank = TaxRank::toRank(lca_rank_string);

      TaxRank::RANK identified_rank = TaxRank::toRank(taxdb.getRank(identified_taxid));
      for (size_t i=0; i < ranks_of_interest.size(); ++i) {
        TaxRank::RANK rank = ranks_of_interest[i];
        TAXID simulated_taxid_at_rank = taxdb.getTaxIDAtRank(seq_taxid, TaxRank::toString(rank));
        TAXID identified_taxid_at_rank = taxdb.getTaxIDAtRank(identified_taxid, TaxRank::toString(rank));
        simulated_taxids_at_rank[rank].insert(simulated_taxid_at_rank);
        // only consider identifications at the rank or more specific
        //  alternative: count identifications that are further up, too
        if (identified_rank <= rank) { 
          identified_taxids_at_rank[rank].insert(identified_taxid_at_rank);
          if (simulated_taxid_at_rank == identified_taxid_at_rank) {
            ++correct_reads_at_rank[rank];
          } else {
            ++incorrect_reads_at_rank[rank];
          }
        } else {
          ++reads_at_higher_rank[rank];
        }
      }

      if (identified_taxid == 0) 
        lca_rank_string = "unidentified";

      ++rank_counts[lca_rank_string];
      out_file 
        << read_id << '\t' << seq_name << '\t' << seq_taxid << '\t'  
        << identified_taxid << '\t' << taxdb.getRank(taxdb.getTaxIDAtNextProperRank(identified_taxid)) << '\t'
        << lca_rank_string << '\t' << lca_taxid_dist.first << '\t' << lca_taxid_dist.second << '\t' << classi << '\n';
    }
  }
  k_file.close();

  char delim = '\t';

  if (0) {
    cout << "#LCA_RANK_READ_COUNTS" << endl;
    for (auto it = rank_counts.begin(); it != rank_counts.end(); ++it) {
      cout << it->first << delim << it->second << endl;
    }
    cout << endl;
  }

  cout << "#rank" << delim << "total_reads" << delim << "correct"<< delim << "incorrect"<< delim << "sensitivity" << delim << "precision"  << delim << "higher_rank" << delim << "unidentified" << endl;
  for (size_t i=0; i < ranks_of_interest.size(); ++i) {
    TaxRank::RANK rank = ranks_of_interest[i];
    size_t true_positives = correct_reads_at_rank.at(rank);
    size_t false_positives = incorrect_reads_at_rank.at(rank);
    double sensitivity = 100.0*(double)true_positives/(double)total_reads;
    double specificity = 100.0*(double)true_positives/(double)(true_positives+false_positives);
    cout << TaxRank::toString(rank) << delim << total_reads 
      << delim << true_positives
      << delim << false_positives
      << delim << sensitivity << '%'
      << delim << specificity << '%'
      << delim << reads_at_higher_rank.at(rank)
      << delim << unidentified_reads 
      << setprecision(2) << std::fixed
      << '\n';
  }

  cout << "#rank" << delim << "true_count" << delim << "correct" << delim << "incorrect" << delim << "recall" << delim << "precision" << endl;
  for (size_t i=0; i < ranks_of_interest.size(); ++i) {
    TaxRank::RANK rank = ranks_of_interest[i];
    size_t true_positives = 0;
    size_t false_positives = 0;
    
    if (identified_taxids_at_rank.find(rank) != identified_taxids_at_rank.end()) {
    for (auto it = identified_taxids_at_rank[rank].begin(); it != identified_taxids_at_rank[rank].end(); ++it) {
      const auto & tid = *it;
      if (simulated_taxids_at_rank[rank].count(tid) == 1) {
        ++true_positives;
      } else {
        ++false_positives;
      }
    }
    }

    double sensitivity = 100.0*(double)true_positives/(double)simulated_taxids_at_rank[rank].size();
    double specificity = 100.0*(double)true_positives/(double)(true_positives+false_positives);

    cout << TaxRank::toString(rank)
      << delim << simulated_taxids_at_rank[rank].size()
      << delim << true_positives
      << delim << false_positives 
      << setprecision(2) << std::fixed
      << delim << sensitivity << '%'
      << delim << specificity << '%'
      << '\n';
  }
}
