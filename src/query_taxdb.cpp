/*
 * Copyright 2017, Florian Breitwieser
 *
 * This file is part of the Kraken taxonomic sequence classification system.
 *
 * Kraken is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Kraken is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kraken.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "taxdb.h"
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <sysexits.h>
#include <cstring>
#include <getopt.h>

using namespace std;

string return_rank;

void process_taxID(char mode, uint32_t taxID);
void process_taxIDs(char mode, vector<uint32_t> taxIDs);
size_t parse_command_line(int argc, char **argv);
void usage(int exit_code=EX_USAGE);

TaxonomyDB<uint32_t, uint32_t> taxdb;

int main(int argc, char **argv) {
  size_t optind = parse_command_line(argc, argv);

  string line;
  uint32_t taxID;
  char mode = *argv[optind++];
	for (;optind < argc; ++optind) {
	  if (strcmp(argv[optind],"-") == 0) {
	    // read STDIN
	    if (mode == 'l') {
	    while (getline(std::cin, line)) {
	      stringstream ss(line);
	      vector<uint32_t> taxIDs;
	      while (ss >> taxID) {
	        taxIDs.push_back(taxID);
	      }
	      process_taxIDs(mode,taxIDs);
	    }
	    }
	    while (std::cin >> taxID) {
	      process_taxID(mode,taxID);
	    }
	  } else {
	    taxID = atol(argv[optind]);
	    process_taxID(mode,taxID);
	  }
	}

	exit(1);
}
void process_taxIDs(char mode, vector<uint32_t> taxIDs) {
  switch (mode) {

  case 'r':
      if (!return_rank.empty()) {
        cout << taxdb.getTaxIDAtRank(taxIDs[0], return_rank) << '\n';
      }
      break;
  case 'l':
    cout << taxdb.getEntry(taxdb.getLowestCommonAncestor(taxIDs)).rank << endl;
    break;
  default:
    usage();
    break;
  }
}


void process_taxID(char mode, uint32_t taxID) {
  switch (mode) {
  case 'r':
    if (!return_rank.empty()) {
      cout << taxdb.getTaxIDAtRank(taxID, return_rank) << '\n';
    }
    break;
  case 'l':
  default:
    usage();
    break;
  }
}

size_t parse_command_line(int argc, char **argv) {
  int opt;
  long long sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);

  while ((opt = getopt(argc, argv, "r:m:")) != -1) {
    switch (opt) {
      case 'r':
        return_rank = optarg;
        break;
      default:
        usage();
        break;
    }
  }

  if (argv[optind] == NULL || argv[optind + 1] == NULL) {
    printf("Mandatory argument(s) missing\n");
    exit(1);
  }

  taxdb.readTaxonomyIndex(argv[optind++], false);
  return optind;
}

void usage(int exit_code) {
  cerr << "Usage: query_taxdb [options] taxDB mode [taxIDs]" << endl
       << endl
       << "Options: (*mandatory)" << endl
       << "  -m mode      Mode: l for LCA, r for rank" << endl
       << "  -r rank      Output parent rank of taxIDs" << endl
       << "  -h           Print this message" << endl
       << endl;
  exit(exit_code);
}

