/*
 * Copyright 2017-2018, Florian Breitwieser
 *
 * This file is part of the KrakenUniq taxonomic sequence classification system.
 *
 * KrakenUniq is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KrakenUniq is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KrakenUniq.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "taxdb.hpp"
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <sysexits.h>
#include <cstring>
#include <getopt.h>

using namespace std;

bool do_print_LCA = false;
bool do_print_Rank = false;
bool do_print_Lineage = true;
string return_rank;

void process_taxID(uint32_t taxID);
void process_taxIDs(const vector<uint32_t>& taxIDs);


void print_LCA(uint32_t taxID);
void print_Rank(uint32_t taxID);
void print_Lineage(uint32_t taxID);

size_t parse_command_line(int argc, char **argv);
void usage(int exit_code=EX_USAGE);

TaxonomyDB<uint32_t> taxdb;

int main(int argc, char **argv) {
  size_t optind = parse_command_line(argc, argv);
  uint32_t taxID;

  for (;optind < static_cast<unsigned>(argc); ++optind) {
    taxID = atoll(argv[optind]);
    process_taxID(taxID);
  }

  std::string line;
  while (!getline(std::cin, line).eof()) {
    stringstream ss(line);
    vector<uint32_t> taxIDs;
	while (ss >> taxID) {
      process_taxID(taxID);
     //taxIDs.push_back(taxID);
	}
    //process_taxIDs(taxIDs);
  }
  exit(0);
}

void process_taxID(uint32_t taxID) {
  if (do_print_Lineage) {
      cout << taxID << '\t' << taxdb.getMetaPhlAnLineage(taxID) << endl;
  } else {
      //cout << taxdb.getTaxIDAtRank(taxID, return_rank) << '\n';
    //cout << taxdb.getEntry(taxdb.getLowestCommonAncestor(taxIDs)).rank << endl;
  }
}

size_t parse_command_line(int argc, char **argv) {
  int opt;
  long long sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);

  while ((opt = getopt(argc, argv, "r:lL")) != -1) {
    switch (opt) {
      case 'r':
		do_print_Rank = true;
        return_rank = optarg;
        break;
	  case 'l':
		do_print_LCA = true;
		break;
	  case 'L':
		do_print_Lineage = true;
	    break; 
      default:
        usage();
        break;
    }
  }

  if (argv[optind] == NULL) {
    printf("Mandatory argument taxDB missing\n");
    exit(1);
  }

  taxdb.readTaxonomyIndex(argv[optind++], false);
  return optind;
}

void usage(int exit_code) {
  cerr << "Usage: query_taxdb [options] taxDB [taxIDs]" << endl
       << endl
       << "Options: (*mandatory)" << endl
       << "  -l           Show LCA of taxIDs" << endl
       << "  -r RANK      Map taxIDs to RANK" << endl
	   << "  -L           Print lineage of taxIDs" << endl
       << "  -h           Print this message" << endl
       << endl;
  exit(exit_code);
}

