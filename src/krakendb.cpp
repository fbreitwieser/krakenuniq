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
#include "krakendb.hpp"
#include "quickfile.hpp"

using std::string;
using std::vector;

namespace kraken {

// File type code for Jellyfish/Kraken DBs
static const char * DATABASE_FILE_TYPE = "JFLISTDN";

// File type code on Kraken DB index
// Next byte determines # of indexed nt
static const char * KRAKEN_INDEX_STRING = "KRAKIDX";

// File type code for Kraken DB index (v2)
// v2 index corresponds to DB sorted by scrambled order
// Next byte determines # of indexed nt
static const char * KRAKEN_INDEX2_STRING = "KRAKIX2";

// XOR mask for minimizer bin keys (allows for better distribution)
// scrambles minimizer sort order
static const uint64_t INDEX2_XOR_MASK = 0xe37e28c4271b5a2dULL;

// Basic constructor
KrakenDB::KrakenDB() {
  fptr = NULL;
  index_ptr = NULL;
  key_ct = 0;
  val_len = 0;
  key_len = 0;
  key_bits = 0;
  k = 0;
}

// Assumes ptr points to start of a readable mmap'ed file
KrakenDB::KrakenDB(char *ptr) {
  index_ptr = NULL;
  fptr = ptr;
  if (strncmp(ptr, DATABASE_FILE_TYPE, strlen(DATABASE_FILE_TYPE)))
    errx(EX_DATAERR, "database in improper format");
  memcpy(&key_bits, ptr + 8, 8);
  memcpy(&val_len, ptr + 16, 8);
  memcpy(&key_ct, ptr + 48, 8);
  if (val_len != 4)
    errx(EX_DATAERR, "can only handle 4 byte DB values");
  k = key_bits / 2;
  key_len = key_bits / 8 + !! (key_bits % 8);
}

// Creates an index, indicating starting positions of each bin
// Bins contain k-mer/taxon pairs with k-mers that share a bin key
void KrakenDB::make_index(string index_filename, uint8_t nt) {
  uint64_t entries = 1ull << (nt * 2);
  vector<uint64_t> bin_counts(entries);
  char *ptr = get_pair_ptr();
  #pragma omp parallel for schedule(dynamic,400)
  for (uint64_t i = 0; i < key_ct; i++) {
    uint64_t kmer = 0;
    memcpy(&kmer, ptr + i * pair_size(), key_len);
    uint64_t b_key = bin_key(kmer, nt);
    #pragma omp atomic
    bin_counts[b_key]++;
  }

  uint64_t *bin_offsets = new uint64_t[ entries + 1 ];
  bin_offsets[0] = 0;
  for (uint64_t i = 1; i <= entries; i++)
    bin_offsets[i] = bin_offsets[i-1] + bin_counts[i-1];

  QuickFile idx_file(index_filename, "w",
    strlen(KRAKEN_INDEX2_STRING) + 1 + sizeof(*bin_offsets) * (entries + 1));
  char *idx_ptr = idx_file.ptr();
  memcpy(idx_ptr, KRAKEN_INDEX2_STRING, strlen(KRAKEN_INDEX2_STRING));
  idx_ptr += strlen(KRAKEN_INDEX2_STRING);
  memcpy(idx_ptr++, &nt, 1);
  memcpy(idx_ptr, bin_offsets, sizeof(*bin_offsets) * (entries + 1));
}

// Simple accessor
char *KrakenDB::get_ptr() {
  return fptr;
}

// Returns start of k-mer/taxon pair array (skips header)
char *KrakenDB::get_pair_ptr() {
  return fptr == NULL ? NULL : fptr + header_size();
}

// Simple accessor
KrakenDBIndex *KrakenDB::get_index() {
  return index_ptr;
}

// Associates the index with this database
void KrakenDB::set_index(KrakenDBIndex *i_ptr) {
  index_ptr = i_ptr;
}

// Simple accessors/convenience methods
uint8_t KrakenDB::get_k() { return k; }
uint64_t KrakenDB::get_key_bits() { return key_bits; }
uint64_t KrakenDB::get_key_len() { return key_len; }
uint64_t KrakenDB::get_val_len() { return val_len; }
uint64_t KrakenDB::get_key_ct() { return key_ct; }
uint64_t KrakenDB::pair_size() { return key_len + val_len; }
size_t KrakenDB::header_size() { return 72 + 2 * (4 + 8 * key_bits); }

// Bin key: each k-mer is made of several overlapping m-mers, m < k
// The bin key is the m-mer whose canonical representation is "smallest"
// ("smallest" can refer to lexico. ordering or some other ordering)
uint64_t KrakenDB::bin_key(uint64_t kmer, uint64_t idx_nt) {
  uint8_t nt = idx_nt;
  uint64_t xor_mask = INDEX2_XOR_MASK;
  uint64_t mask = 1 << (nt * 2);
  mask--;
  xor_mask &= mask;
  uint64_t min_bin_key = ~0;
  for (uint64_t i = 0; i < key_bits / 2 - nt + 1; i++) {
    uint64_t temp_bin_key = xor_mask ^ canonical_representation(kmer & mask, nt);
    if (temp_bin_key < min_bin_key)
      min_bin_key = temp_bin_key;
    kmer >>= 2;
  }
  return min_bin_key;
}

// Separate functions to avoid a conditional in the function
// This probably isn't necessary...
uint64_t KrakenDB::bin_key(uint64_t kmer) {
  uint8_t nt = index_ptr->indexed_nt();
  uint8_t idx_type = index_ptr->index_type();
  uint64_t xor_mask = idx_type == 1 ? 0 : INDEX2_XOR_MASK;
  uint64_t mask = 1 << (nt * 2);
  mask--;
  xor_mask &= mask;
  uint64_t min_bin_key = ~0;
  for (uint64_t i = 0; i < key_bits / 2 - nt + 1; i++) {
    uint64_t temp_bin_key = xor_mask ^ canonical_representation(kmer & mask, nt);
    if (temp_bin_key < min_bin_key)
      min_bin_key = temp_bin_key;
    kmer >>= 2;
  }
  return min_bin_key;
}

// Code mostly from Jellyfish 1.6 source
uint64_t KrakenDB::reverse_complement(uint64_t kmer, uint8_t n) {
  kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
  kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
  kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
  kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
  kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
  return (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - (n << 1));
}

// Code mostly from Jellyfish 1.6 source
uint64_t KrakenDB::reverse_complement(uint64_t kmer) {
  kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
  kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
  kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
  kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
  kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
  return (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - (k << 1));
}

// Lexicographically smallest of k-mer and reverse comp. of k-mer
uint64_t KrakenDB::canonical_representation(uint64_t kmer, uint8_t n) {
  uint64_t revcom = reverse_complement(kmer, n);
  return kmer < revcom ? kmer : revcom;
}

uint64_t KrakenDB::canonical_representation(uint64_t kmer) {
  uint64_t revcom = reverse_complement(kmer, k);
  return kmer < revcom ? kmer : revcom;
}

// perform search over last range to speed up queries
// NOTE: retry_on_failure implies all pointer params are non-NULL
uint32_t *KrakenDB::kmer_query(uint64_t kmer, uint64_t *last_bin_key,
                               int64_t *min_pos, int64_t *max_pos,
                               bool retry_on_failure)
{
  int64_t min, max, mid;
  uint64_t comp_kmer;
  uint64_t b_key;
  char *ptr = get_pair_ptr();
  size_t pair_sz = pair_size();

  // Use provided values if they exist and are valid
  if (retry_on_failure && *min_pos <= *max_pos) {
    b_key = *last_bin_key;
    min = *min_pos;
    max = *max_pos;
  }
  else {
    b_key = bin_key(kmer);
    min = index_ptr->at(b_key);
    max = index_ptr->at(b_key + 1) - 1;
    // Invalid min/max values + retry_on_failure means min/max need to be
    // initialized and set in caller
    if (retry_on_failure) {
      *last_bin_key = b_key;
      *min_pos = min;
      *max_pos = max;
    }
  }

  // Binary search with large window
  while (min + 15 <= max) {
    mid = min + (max - min) / 2;
    comp_kmer = 0;
    memcpy(&comp_kmer, ptr + pair_sz * mid, key_len);
    comp_kmer &= (1ull << key_bits) - 1;  // trim any excess
    if (kmer > comp_kmer)
      min = mid + 1;
    else if (kmer < comp_kmer)
      max = mid - 1;
    else
      return (uint32_t *) (ptr + pair_sz * mid + key_len);
  }
  // Linear search once window shrinks
  for (mid = min; mid <= max; mid++) {
    comp_kmer = 0;
    memcpy(&comp_kmer, ptr + pair_sz * mid, key_len);
    comp_kmer &= (1ull << key_bits) - 1;  // trim any excess
    if (kmer == comp_kmer)
      return (uint32_t *) (ptr + pair_sz * mid + key_len);
  }

  uint32_t *answer = NULL;
  // ROF implies the provided values might be out of date
  // If they are, we'll update them and search again
  if (retry_on_failure) {
    b_key = bin_key(kmer);
    // If bin key hasn't changed, search fails
    if (b_key == *last_bin_key)
      return NULL;
    min = index_ptr->at(b_key);
    max = index_ptr->at(b_key + 1) - 1;
    // Recursive call w/ adjusted search params and w/o retry
    answer = kmer_query(kmer, &b_key, &min, &max, false);
    // Update caller's search params due to bin key change
    if (last_bin_key != NULL) {
      *last_bin_key = b_key;
      *min_pos = min;
      *max_pos = max;
    }
  }
  return answer;
}

// Binary search w/in the k-mer's bin
uint32_t *KrakenDB::kmer_query(uint64_t kmer) {
  return kmer_query(kmer, NULL, NULL, NULL, false);
}

KrakenDBIndex::KrakenDBIndex() {
  fptr = NULL;
  idx_type = 1;
  nt = 0;
}

KrakenDBIndex::KrakenDBIndex(char *ptr) {
  fptr = ptr;
  idx_type = 1;
  if (strncmp(ptr, KRAKEN_INDEX_STRING, strlen(KRAKEN_INDEX_STRING))) {
    idx_type = 2;
    if (strncmp(ptr, KRAKEN_INDEX2_STRING, strlen(KRAKEN_INDEX2_STRING)))
      errx(EX_DATAERR, "illegal Kraken DB index format");
  }
  ptr += strlen(KRAKEN_INDEX_STRING);
  memcpy(&nt, ptr, 1);
}

// Index version (v2 uses different minimizer sort order)
uint8_t KrakenDBIndex::index_type() {
  return idx_type;
}

// How long are bin keys (i.e., what is minimizer length in bp?)
uint8_t KrakenDBIndex::indexed_nt() {
  return nt;
}

// Return start of index array (skips header)
uint64_t *KrakenDBIndex::get_array() {
  return (uint64_t *) (fptr + strlen(KRAKEN_INDEX_STRING) + 1);
}

// Convenience method, allows for testing guard
uint64_t KrakenDBIndex::at(uint64_t idx) {
  uint64_t *array = get_array();
  #ifdef TESTING
  if (idx > 1 + (1ull << (nt * 2)))
    errx(EX_SOFTWARE, "KrakenDBIndex::at() called with illegal index");
  #endif
  return array[idx];
}

} // namespace
