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

KHASH_SET_INIT_INT64(set64);
using u64 = uint64_t;
struct khset64_t: khash_t(set64) {
    khset64_t() {std::memset(this, 0, sizeof(*this));}
    ~khset64_t() {std::free(this->flags); std::free(this->vals); std::free(this->keys);}
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
    left.insert(kh_val(right, it));
  }
  return left;
}

#endif