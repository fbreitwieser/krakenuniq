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

#include <cstdint>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include "third_party/khash.h"

KHASH_SET_INIT_INT(set)
KHASH_SET_INIT_INT64(set64);

using u32 = ::std::uint32_t;
using u64 = std::uint64_t;

#define __FE__ \
    template<typename Functor>\
    void for_each(const Functor &func) {\
        for(khiter_t ki = 0; ki < this->n_buckets; ++ki)\
            if(kh_exist(this, ki))\
                func(this->keys[ki]);\
    }\
    template<typename Functor>\
    void for_each(const Functor &func) const {\
        for(khiter_t ki = 0; ki < this->n_buckets; ++ki)\
            if(kh_exist(this, ki))\
                func(this->keys[ki]);\
    }


// Steal everything, take no prisoners.
#define MOVE_DEC(t) \
   t(t &&other) {std::memcpy(this, &other, sizeof(*this)); std::memset(&other, 0, sizeof(other));}\
   t& operator=(t &&other) {std::memcpy(this, &other, sizeof(*this)); std::memset(&other, 0, sizeof(other)); return (*this);}

#define COPY_DEC(t) \
    t(const t &other) {\
        if(other.size()) {\
            std::memcpy(this, &other, sizeof(*this));\
            auto memsz = other.size() * sizeof(*keys);\
            keys = static_cast<decltype(keys)>(std::malloc(memsz));\
            if(!keys) throw std::bad_alloc();\
            std::memcpy(keys, other.keys, memsz);\
            memsz = __ac_fsize(other.size()) * sizeof(u32);\
            flags = static_cast<u32 *>(std::malloc(memsz));\
            if(!flags) throw std::bad_alloc();\
            std::memcpy(flags, other.flags, memsz);\
        } else std::memset(this, 0, sizeof(*this));\
    }

struct khset32_t: khash_t(set) {
    khset32_t() {std::memset(this, 0, sizeof(*this));}
    ~khset32_t() {std::free(this->flags); std::free(this->keys);}
    COPY_DEC(khset32_t)
    MOVE_DEC(khset32_t)
    // For each
    __FE__
    operator khash_t(set) &() {return *reinterpret_cast<khash_t(set) *>(this);}
    operator khash_t(set) *() {return reinterpret_cast<khash_t(set) *>(this);}
    operator const khash_t(set) &() const {return *reinterpret_cast<const khash_t(set) *>(this);}
    operator const khash_t(set) *() const {return reinterpret_cast<const khash_t(set) *>(this);}
    khash_t(set) *operator->() {return static_cast<khash_t(set) *>(*this);}
    const khash_t(set) *operator->() const {return static_cast<const khash_t(set) *>(*this);}
    auto insert(u32 val) {
        int khr;
        return kh_put(set, this, val, &khr);
    }
    template<typename ItType>
    void insert(ItType i1, ItType i2) {
        while(i1 < i2) insert(*i1++);
    }
    void clear() {kh_clear(set, this);}
    bool contains(u32 x) const {return kh_get(set, this, x) != kh_end(this);}
    size_t size() const {return kh_size(static_cast<const khash_t(set) *>(this));}
};


struct khset64_t: khash_t(set64) {
    khset64_t() {std::memset(this, 0, sizeof(*this));}
    ~khset64_t() {std::free(this->flags); std::free(this->keys);}
    COPY_DEC(khset64_t)
    MOVE_DEC(khset64_t)
    // For each
    __FE__
    operator khash_t(set64) &() {return *reinterpret_cast<khash_t(set64) *>(this);}
    operator khash_t(set64) *() {return reinterpret_cast<khash_t(set64) *>(this);}
    operator const khash_t(set64) &() const {return *reinterpret_cast<const khash_t(set64) *>(this);}
    operator const khash_t(set64) *() const {return reinterpret_cast<const khash_t(set64) *>(this);}
    khash_t(set64) *operator->() {return static_cast<khash_t(set64) *>(*this);}
    const khash_t(set64) *operator->() const {return static_cast<const khash_t(set64) *>(*this);}
    auto insert(u64 val) {
        int khr;
        return kh_put(set64, *this, val, &khr);
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

#undef __FE__
#undef COPY_DEC
#undef MOVE_DEC
#endif