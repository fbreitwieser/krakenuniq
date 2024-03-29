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

#ifndef KRAKENDB_HPP
#define KRAKENDB_HPP

#include "kraken_headers.hpp"
#include <unordered_map>
#include <map>

namespace kraken {
  class KrakenDB;

  class KrakenDBIndex {
    public:
    KrakenDBIndex();
    // ptr points to mmap'ed existing file opened in read or read/write mode
    KrakenDBIndex(char *ptr);

    uint8_t index_type();
    uint8_t indexed_nt();
    uint64_t mmap_at(uint64_t idx);
    uint64_t *get_array_with_db_chunks();
    uint64_t at_with_db_chunks(uint64_t idx);
    uint64_t *get_array();
    uint64_t at(uint64_t idx);

    friend KrakenDB;

    private:
    uint8_t idx_type;
    char *fptr;
    uint8_t nt;

    char *data;
    uint64_t data_offset;
  };

  class KrakenDB {
    public:

    char *get_ptr();            // Return the file pointer
    char *get_pair_ptr();       // Return pointer to start of pairs
    KrakenDBIndex *get_index(); // Return ptr to assoc'd index obj
    uint8_t get_k();            // how many nt are in each key?
    uint64_t get_key_bits();    // how many bits are in each key?
    uint64_t get_key_len();     // how many bytes does each key occupy?
    uint64_t get_val_len();     // how many bytes does each value occupy?
    uint64_t get_key_ct();      // how many key/value pairs are there?
    uint64_t pair_size();       // how many bytes does each pair occupy?

    size_t header_size();  // Jellyfish uses variable header sizes
    uint32_t *kmer_query(uint64_t kmer);  // return ptr to pair w/ kmer

    // perform search over last range to speed up queries
    uint32_t *kmer_query(uint64_t kmer, uint64_t *last_bin_key,
                         int64_t *min_pos, int64_t *max_pos,
                         bool retry_on_failure=true);

    uint32_t *kmer_query_with_db_chunks(uint64_t kmer);  // return ptr to pair w/ kmer

    // perform search over last range to speed up queries
    uint32_t *kmer_query_with_db_chunks(uint64_t kmer, uint64_t *last_bin_key,
                                        int64_t *min_pos, int64_t *max_pos,
                                        bool retry_on_failure=true);

    // return a count of k-mers for all taxons
    std::map<uint32_t,uint64_t> count_taxons();
    
    // return "bin key" for kmer, based on index
    // If idx_nt not specified, use index's value
    uint64_t bin_key(uint64_t kmer, uint64_t idx_nt);
    uint64_t bin_key(uint64_t kmer);

    // Code from Jellyfish, rev. comp. of a k-mer with n nt.
    // If n is not specified, use k in DB, otherwise use first n nt in kmer
    uint64_t reverse_complement(uint64_t kmer, uint8_t n);
    uint64_t reverse_complement(uint64_t kmer);

    // Return lexicographically smallest of kmer/revcom(kmer)
    // If n is not specified, use k in DB, otherwise use first n nt in kmer
    uint64_t canonical_representation(uint64_t kmer, uint8_t n);
    uint64_t canonical_representation(uint64_t kmer);

    void make_index(std::string index_filename, uint8_t nt);

    void set_index(KrakenDBIndex *i_ptr);

    size_t filesize() const;

    // Null constructor
    KrakenDB();

    // ptr points to start of mmap'ed DB in read or read/write mode
    KrakenDB(char *ptr, size_t filesize = 0);

    ~KrakenDB();

    uint32_t chunks() const;
    void load_chunk(const uint32_t db_chunk_id);
    void prepare_chunking(const uint64_t max_bytes_for_db);
    bool is_minimizer_in_chunk(const uint64_t minimizer, const uint32_t db_chunk_id) const;

    private:

    size_t _filesize;
    char *fptr;
    KrakenDBIndex *index_ptr;
    uint8_t k;
    uint64_t key_bits;
    uint64_t key_len;
    uint64_t val_len;
    uint64_t key_ct;

    uint64_t upper_bound(const uint64_t first, const uint64_t last);

    uint32_t _chunks;
    std::vector<uint64_t> idx_chunk_bounds;
    std::vector<uint64_t> dbx_chunk_bounds;

    char *data;
    uint64_t data_size;
    uint64_t data_offset;
  };
}

#endif
