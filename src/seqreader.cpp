/*
 * Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
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

#include "kraken_headers.hpp"
#include "seqreader.hpp"

using namespace std;

namespace kraken {
  FastaReader::FastaReader(string filename) {
    file.open(filename.c_str());
    if (file.rdstate() & ifstream::failbit) {
      err(EX_NOINPUT, "can't open %s", filename.c_str());
    }
    valid = true;
  }

  DNASequence FastaReader::next_sequence() {
    DNASequence dna;

    if (! valid || ! file.good()) {
      valid = false;
      return dna;
    }
    string line;

    if (linebuffer.empty()) {
      getline(file, line);
    }
    else {
      line = linebuffer;
      linebuffer.clear();
    }

    if (line[0] != '>') {
      warnx("malformed fasta file - expected header char > not found");
      valid = false;
      return dna;
    }
    dna.header_line = line.substr(1);
    istringstream seq_id(dna.header_line);
    seq_id >> dna.id;
    
    ostringstream seq_ss;

    while (file.good()) {
      getline(file, line);
      if (line[0] == '>') {
        linebuffer = line;
        break;
      }
      else {
        seq_ss << line;
      }
    }
    dna.seq = seq_ss.str();

    if (dna.seq.empty()) {
      warnx("malformed fasta file - zero-length record (%s)", dna.id.c_str());
      valid = false;
      return dna;
    }

    return dna;
  }

  bool FastaReader::is_valid() {
    return valid;
  }

  FastqReader::FastqReader(string filename) {
    file.open(filename.c_str());
    if (file.rdstate() & ifstream::failbit) {
      err(EX_NOINPUT, "can't open %s", filename.c_str());
    }
    valid = true;
  }

  DNASequence FastqReader::next_sequence() {
    DNASequence dna;

    if (! valid || ! file.good()) {
      valid = false;
      return dna;
    }

    string line;
    getline(file, line);
    if (line.empty()) {
      valid = false;  // Sometimes FASTQ files have empty last lines
      return dna;
    }
    if (line[0] != '@') {
      if (line[0] != '\r')
        warnx("malformed fastq file - sequence header (%s)", line.c_str());
      valid = false;
      return dna;
    }
    dna.header_line = line.substr(1);
    istringstream line_ss(dna.header_line);
    
    line_ss >> dna.id;
    getline(file, dna.seq);

    getline(file, line);
    if (line.empty() || line[0] != '+') {
      if (line[0] != '\r')
        warnx("malformed fastq file - quality header (%s)", line.c_str());
      valid = false;
      return dna;
    }
    getline(file, dna.quals);

    return dna;
  }

  bool FastqReader::is_valid() {
    return valid;
  }
} // namespace
