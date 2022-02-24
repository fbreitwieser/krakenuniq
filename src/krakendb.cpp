/*
 * Portions (c) 2017-2018, Florian Breitwieser <fbreitwieser@jhu.edu> as part of KrakenUniq
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
#include <unordered_map>

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
  _filesize = 0;
}

// Assumes ptr points to start of a readable mmap'ed file
KrakenDB::KrakenDB(char *ptr, size_t filesize) {
  _filesize = filesize;
  index_ptr = NULL;
  fptr = ptr;
  if (ptr == NULL) {
    errx(EX_DATAERR, "pointer is NULL");
  }
  if (strncmp(ptr, DATABASE_FILE_TYPE, strlen(DATABASE_FILE_TYPE))) {
    errx(EX_DATAERR,"database in improper format - found %s", string(ptr, strlen(DATABASE_FILE_TYPE)).c_str());
  }
  memcpy(&key_bits, ptr + 8, 8);
  memcpy(&val_len, ptr + 16, 8);
  memcpy(&key_ct, ptr + 48, 8);
  if (val_len != 4)
    errx(EX_DATAERR, "can only handle 4 byte DB values");
  k = key_bits / 2;
  key_len = key_bits / 8 + !! (key_bits % 8);
  std::cerr << "Loaded database with " << key_ct << " keys with k of " << (size_t)k << " [val_len " << val_len << ", key_len " << key_len << "]." << std::endl;
}

// destructor
KrakenDB::~KrakenDB() {
  munmap(data, data_size);
}

size_t KrakenDB::filesize() const {
  return _filesize;
}

//using std::map to have the keys sorted
std::map<uint32_t,uint64_t> KrakenDB::count_taxons() {
  char *ptr = get_pair_ptr();
  size_t pair_sz = pair_size();

  std::map<uint32_t, uint64_t> taxon_counts;
  if (ptr == NULL) { 
    std::cerr << "Kraken database pointer is NULL [pair_sz: " << pair_sz << ", key_ct: "<<key_ct<<", key_len: "<< key_len<<"]!" << std::endl;
    exit(1);
  }
  for (uint64_t i = 0; i < key_ct; i++) {
    if (i % 10000000 == 1) {
      fprintf(stderr, "\r %.2f %%", 100*(double)i/(double)key_ct );
    }
    uint32_t* taxon = (uint32_t *) (ptr + pair_sz * i + key_len);
    if (taxon == NULL) {
        std::cerr << "taxon is NULL (i is " << i << " and key_ct is " << key_ct << ")" << std::endl;
    } else {
        uint32_t taxon_i = *taxon;
        ++taxon_counts[taxon_i];
    }
  }
  fprintf(stderr, "\r");
  return taxon_counts;
}


// Creates an index, indicating starting positions of each bin
// Bins contain k-mer/taxon pairs with k-mers that share a bin key
void KrakenDB::make_index(string index_filename, uint8_t nt) {
  uint64_t entries = 1ull << (nt * 2);
  vector<uint64_t> bin_counts(entries);
  char *ptr = get_pair_ptr();

#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic,400)
#endif
  for (uint64_t i = 0; i < key_ct; i++) {
    uint64_t kmer = 0;
    memcpy(&kmer, ptr + i * pair_size(), key_len);
    uint64_t b_key = bin_key(kmer, nt);
#ifdef _OPENMP
    #pragma omp atomic
#endif
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

// perform search over last range to speed up queries
// NOTE: retry_on_failure implies all pointer params are non-NULL
uint32_t *KrakenDB::kmer_query_with_db_chunks(uint64_t kmer, uint64_t *last_bin_key,
                               int64_t *min_pos, int64_t *max_pos,
                               bool retry_on_failure)
{
  int64_t min, max, mid;
  uint64_t comp_kmer;
  uint64_t b_key;
  size_t pair_sz = pair_size();

  // Use provided values if they exist and are valid
  if (retry_on_failure && *min_pos <= *max_pos) {
    b_key = *last_bin_key;
    min = *min_pos;
    max = *max_pos;
  }
  else {
    b_key = bin_key(kmer);
    min = index_ptr->at_with_db_chunks(b_key);
    max = index_ptr->at_with_db_chunks(b_key + 1) - 1;
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
    memcpy(&comp_kmer, data + pair_sz * mid - data_offset, key_len);
    comp_kmer &= (1ull << key_bits) - 1;  // trim any excess
    if (kmer > comp_kmer)
      min = mid + 1;
    else if (kmer < comp_kmer)
      max = mid - 1;
    else
      return (uint32_t *) (data + pair_sz * mid - data_offset + key_len);
  }
  // Linear search once window shrinks
  for (mid = min; mid <= max; mid++) {
    comp_kmer = 0;
    memcpy(&comp_kmer, data + pair_sz * mid - data_offset, key_len);
    comp_kmer &= (1ull << key_bits) - 1;  // trim any excess
    if (kmer == comp_kmer)
      return (uint32_t *) (data + pair_sz * mid - data_offset + key_len);
  }

  uint32_t *answer = NULL;
  // ROF implies the provided values might be out of date
  // If they are, we'll update them and search again
  if (retry_on_failure) {
    b_key = bin_key(kmer);
    // If bin key hasn't changed, search fails
    if (b_key == *last_bin_key)
      return NULL;
    min = index_ptr->at_with_db_chunks(b_key);
    max = index_ptr->at_with_db_chunks(b_key + 1) - 1;
    // Recursive call w/ adjusted search params and w/o retry
    answer = kmer_query_with_db_chunks(kmer, &b_key, &min, &max, false);
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
uint32_t *KrakenDB::kmer_query_with_db_chunks(uint64_t kmer) {
  return kmer_query_with_db_chunks(kmer, NULL, NULL, NULL, false);
}

uint32_t KrakenDB::chunks() const {
  return this->_chunks;
}

void KrakenDB::load_chunk(const uint32_t db_chunk_id) {
  char* db_chunk_start = get_pair_ptr() + dbx_chunk_bounds[db_chunk_id] * pair_size();
  data_offset = dbx_chunk_bounds[db_chunk_id] * pair_size();
  const size_t db_chunk_len = (dbx_chunk_bounds[db_chunk_id + 1] - dbx_chunk_bounds[db_chunk_id]) * pair_size();
  memcpy(data, db_chunk_start, db_chunk_len);

  char* id_chunk_start = index_ptr->fptr + strlen(KRAKEN_INDEX_STRING) + 1 + (idx_chunk_bounds[db_chunk_id] * 8);
  index_ptr->data_offset = idx_chunk_bounds[db_chunk_id] * 8;
  const size_t id_chunk_len = (idx_chunk_bounds[db_chunk_id + 1] - idx_chunk_bounds[db_chunk_id]) * 8;

  index_ptr->data = data + db_chunk_len;
  memcpy(index_ptr->data, id_chunk_start, id_chunk_len);
}

// returns first element LARGER
// half-open interval. first is inclusive, last is exclusive
// the return value has to be treated as inclusive (closed interval)
uint64_t KrakenDB::upper_bound(uint64_t first, const uint64_t last)
{
  uint64_t it;
  uint64_t count, step;
  count = last - first;

  const uint64_t orig_first = first;
  const uint64_t data_orig_first = index_ptr->mmap_at(orig_first);

  // printf("upper_bound(%llu, %llu)\n", orig_first, last);
  // float max_size = 1.0f * data_size / (1024 * 1024 * 1024);
  // printf("- max_size = %f\n", max_size);

  while (count > 0) {
    it = first;
    step = count / 2;
    it += step;

    uint64_t size_index = (it + 1 - orig_first) * 8;
    uint64_t size_data = (index_ptr->mmap_at(it + 1) - data_orig_first) * pair_size();
    // float total_size_gb = 1.0f * (size_index + size_data) / (1024 * 1024 * 1024);

    // printf("- it = %llu (size: %f)\n", it, total_size_gb);
    if (size_index + size_data <= data_size) { // if (!(data_size < size_index + size_data)) // if (!comp(value, *it))
      first = ++it;
      count -= step + 1;
    }
    else
      count = step;
  }
  return first;
}

void KrakenDB::prepare_chunking(const uint64_t max_bytes_for_db) {
  data = (char*) mmap(0, max_bytes_for_db, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  if (data == (void *) -1)
  {
    printf("ERROR: Could not allocated memory for loading database chunks (error code %i: %s)\n", errno, strerror(errno));
    exit(-1);
  }
  data_size = max_bytes_for_db;
  mlock(data, data_size);

  this->_chunks = 0;

  const uint64_t max_idx_pos = 1ull << (2 * index_ptr->indexed_nt()); // excl bound (half-open interval)
  // const uint64_t max_dbx_pos = index_ptr->mmap_at(max_idx_pos); // excl bound (half-open interval)
  uint64_t idx_pos = 0;
  uint64_t dbx_pos = 0; // == index_ptr->mmap_at(0);
  idx_chunk_bounds.push_back(idx_pos);
  dbx_chunk_bounds.push_back(dbx_pos);

  while (idx_pos < max_idx_pos)
  {
    // upper_bound returns element idx_pos that is the first that "is greater" (i.e., that does not fit into "data")
    // in other terms, idx_pos closes our half-open interval (excl. upper bound), hence we don't need to decrement it
    idx_pos = upper_bound(idx_pos, max_idx_pos);
    // get starting position for next larger minimizer (that does not fit into "data")
    // in other terms, dbx_pos closes our half-open interval (excl. upper bound), hence we don't need to decrement it
    dbx_pos = index_ptr->mmap_at(idx_pos);

    // uint64_t idx_size = (idx_pos - idx_chunk_bounds.back()) * 8;
    // uint64_t dbx_size = (dbx_pos - dbx_chunk_bounds.back()) * pair_size();
    // printf("- %" PRIu64 ": %s\n", idx_size + dbx_size, (idx_size + dbx_size <= data_size ? "ok" : "NOT OK, does not fit into data"));

    // especially the last chunk can have 0 k-mers but a huge part of minimizers in it (that did not occur in any k-mer in the database)
    // in that case we can skip that chunk!
    if (dbx_pos == dbx_chunk_bounds.back())
      continue;

    idx_chunk_bounds.push_back(idx_pos);
    dbx_chunk_bounds.push_back(dbx_pos);

    this->_chunks++;
  }

  // printf("Chunk bounds\n");
  // for (uint32_t i = 0; i < this->_chunks + 1; ++i)
  //   printf("%llu\t%llu\n", idx_chunk_bounds[i], dbx_chunk_bounds[i]);

  // printf("Chunk sizes\n");
  // for (uint32_t i = 1; i < this->_chunks + 1; ++i) {
    // uint64_t idx_size = (idx_chunk_bounds[i] - idx_chunk_bounds[i - 1]) * 8;
    // uint64_t dbx_size = (dbx_chunk_bounds[i] - dbx_chunk_bounds[i - 1]) * pair_size();
    // float idx_size_gb = 1.0f * idx_size / (1024 * 1024 * 1024);
    // float dbx_size_gb = 1.0f * dbx_size / (1024 * 1024 * 1024);
    // printf("Chunk %u: %.2f GB\t%.2f GB (upper bounds: %llu\t%llu)\n", i - 1, idx_size_gb, dbx_size_gb, idx_chunk_bounds[i], dbx_chunk_bounds[i]);
    // if (idx_size + dbx_size > data_size) {
    //   printf("ERROR: %llu > %llu\n", idx_size + dbx_size, data_size);
    //   exit(1);
    // }
  // }
}

bool KrakenDB::is_minimizer_in_chunk(const uint64_t minimizer, const uint32_t db_chunk_id) const {
    return idx_chunk_bounds[db_chunk_id] <= minimizer && minimizer < idx_chunk_bounds[db_chunk_id + 1];
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

// Convenience method, allows for testing guard
uint64_t KrakenDBIndex::mmap_at(uint64_t idx) {
  uint64_t *array = (uint64_t*) (fptr + strlen(KRAKEN_INDEX_STRING) + 1);
  #ifdef TESTING
  if (idx > 1 + (1ull << (nt * 2)))
    errx(EX_SOFTWARE, "KrakenDBIndex::at() called with illegal index");
  #endif
  return array[idx];
}

// Return start of index array (skips header)
uint64_t *KrakenDBIndex::get_array_with_db_chunks() { // TODO: remove wrapper
  return (uint64_t *) data;
}

uint64_t KrakenDBIndex::at_with_db_chunks(uint64_t idx) {
  uint64_t *array = get_array_with_db_chunks();
#ifdef TESTING
  if ((idx - (data_offset / 8)) > 1 + (1ull << (nt * 2)))
errx(EX_SOFTWARE, "KrakenDBIndex::at() called with illegal index");
#endif
  return array[idx - (data_offset / 8)];
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
