/*
 * Copyright 2013, Derrick Wood <dwood@cs.umd.edu>
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
#include "quickfile.hpp"
#include "seqreader.hpp"

using namespace std;

namespace kraken {
  FastaReader::FastaReader(string filename) {
    file.open_file(filename);
    ptr = file.ptr();
    valid = true;
  }

  DNASequence FastaReader::next_sequence() {
    DNASequence dna;

    if (! valid || ! ptr || (size_t)(ptr - file.ptr()) >= file.size()) {
      valid = false;
      return dna;
    }

    size_t remaining = file.size() - (ptr - file.ptr());
    char *end_ptr = ptr;
    size_t end_rem = remaining;
    while (end_ptr) {
      end_ptr = (char *) memchr(end_ptr+1, '>', end_rem-1);
      if (! end_ptr || *(end_ptr - 1) == '\n')
        break;
      end_rem = remaining - (end_ptr - ptr);
    }

    size_t record_len = end_ptr ? end_ptr - ptr : remaining;
    string fasta_record(ptr, record_len);
    ptr += record_len;

    istringstream iss(fasta_record);
    string line;

    getline(iss, line);
    if (line[0] != '>')
      errx(EX_DATAERR, "malformed fasta file");
    istringstream line_ss(line.substr(1));
    
    line_ss >> dna.id;
    dna.seq = "";
    while (iss.good()) {
      getline(iss, line);
      if (line.empty())
        break;
      dna.seq.append(line);
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
    istringstream line_ss(line.substr(1));
    
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
