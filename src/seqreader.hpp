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

#ifndef SEQREADER_HPP
#define SEQREADER_HPP

#include "kraken_headers.hpp"

namespace kraken {
  typedef struct {
    std::string id;
    std::string header_line;  // id + optional description
    std::string seq;
    std::string quals;
  } DNASequence;

  class DNASequenceReader {
    public:
    virtual DNASequence next_sequence() = 0; 
    virtual bool is_valid() = 0;
    virtual ~DNASequenceReader() {}
  };

  class FastaReader : public DNASequenceReader {
    public:
    FastaReader(std::string filename);
    DNASequence next_sequence();
    bool is_valid();

    private:
    std::ifstream file;
    std::string linebuffer;
    bool valid;
  };

  class FastqReader : public DNASequenceReader {
    public:
    FastqReader(std::string filename);
    DNASequence next_sequence();
    bool is_valid();

    private:
    std::ifstream file;
    bool valid;
  };
}

#endif
