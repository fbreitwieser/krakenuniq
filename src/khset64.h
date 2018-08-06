/* 
 * Copyright 2018 Daniel N Baker 
 *
 * The file is part of KrakenHLL
 *
 * KrakenHLL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KrakenHLL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.

 * You should have received a copy of the GNU Affero General Public License
 * along with KrakenHLL.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef KHSET64_H_
#define KHSET64_H_

#include "third_party/khash.h"
#include <cstring>
#include <algorithm>
#include <iostream>
#include <cstdint>

using namespace std;
KHASH_SET_INIT_INT64(set64);
using u64 = uint64_t;
struct khset64_t: khash_t(set64) {
    khset64_t() {std::memset(this, 0, sizeof(*this));}

    khset64_t(const khset64_t& k) {
	  this->n_buckets = k.n_buckets;
	  khash_t(set64)::size = k.size();
	  this->n_occupied = k.n_occupied;
	  this->upper_bound = k.upper_bound;
	  this->flags = (khint32_t*)kmalloc(n_buckets * sizeof(khint32_t));	
	  std::copy_n(k.flags, n_buckets, this->flags);
	  cerr << "khset copy const flags" << endl;
	  this->keys = (u64*)kmalloc(n_buckets * sizeof(u64));
	  std::copy_n(k.keys, n_buckets, this->keys);
	  //std::copy_n(k.vals, n_buckets, this->vals);
	  cerr << "khset copy const end" << endl;
	}

    khset64_t(khset64_t&& k) {
	  cerr << "khset rval const start" << endl;
	  this->n_buckets = k.n_buckets;
	  khash_t(set64)::size = k.size();
	  this->n_occupied = k.n_occupied;
	  this->upper_bound = k.upper_bound;
	  std::copy(k.flags, k.flags+size(), this->flags);
	  std::copy(k.keys, k.keys+size(), this->keys);
	  std::copy(k.vals, k.vals+size(), this->vals);
	  cerr << "khset rval const end" << endl;
	}

    khset64_t& operator=(khset64_t&& k) {
	  cerr << "khset rval = start" << endl;
	  this->n_buckets = k.n_buckets;
	  khash_t(set64)::size = k.size();
	  this->n_occupied = k.n_occupied;
	  this->upper_bound = k.upper_bound;
	  std::copy(k.flags, k.flags+size(), this->flags);
	  std::copy(k.keys, k.keys+size(), this->keys);
	  std::copy(k.vals, k.vals+size(), this->vals);
	  cerr << "khset rval = start" << endl;
	  return *this;
	}
    /*~khset64_t() {
			std::free(this->flags); 
			std::free(this->vals); 
			std::free(this->keys);
	}*/
    operator khash_t(set64) &() {return *reinterpret_cast<khash_t(set64) *>(this);}
    operator khash_t(set64) *() {return reinterpret_cast<khash_t(set64) *>(this);}
    operator const khash_t(set64) &() const {return *reinterpret_cast<const khash_t(set64) *>(this);}
    operator const khash_t(set64) *() const {return reinterpret_cast<const khash_t(set64) *>(this);}
    khash_t(set64) *operator->() {return static_cast<khash_t(set64) *>(*this);}
    const khash_t(set64) *operator->() const {return static_cast<const khash_t(set64) *>(*this);}
    void insert(u64 val) {
        int khr;
        kh_put(set64, *this, val, &khr);
    }
    void clear() {kh_clear(set64, this);}
    bool contains(u64 x) const {return kh_get(set64, this, x) != kh_end(this);}
    size_t size() const {return kh_size(static_cast<const khash_t(set64) *>(this));}
};

khset64_t& operator+=(khset64_t& left, const khset64_t& right) {
  for (auto it = kh_begin(right); it != kh_end(right); ++it) {
	if (kh_exist(right, it))
      left.insert(kh_key(right,it));
  }
  return left;
}

#endif