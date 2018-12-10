/* 
 * Copyright (c) 2018 Daniel N Baker, <dnb@jhu.edu>
 *
 * The file is part of KrakenUniq
 *
 * KrakenUniq is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KrakenUniq is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.

 * You should have received a copy of the GNU Affero General Public License
 * along with KrakenUniq.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef KHSET_H
#define KHSET_H

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include "klib/khash.h"

KHASH_SET_INIT_INT(set)
KHASH_SET_INIT_INT64(set64)

using u32 = ::std::uint32_t;
using u64 = ::std::uint64_t;

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
// re memset reinterpret_cast:
//   khset is not trivial, so we cannot clear it with memset(3), 
//   but zeroing objects like that is safe and should be find for our purposes 
#define MOVE_DEC(t) \
   t(t &&other) { \
       std::memcpy(this, &other, sizeof(*this)); \
       std::memset(reinterpret_cast<void*>(&other), 0, sizeof(other));\
   }

#define COPY_DEC(t) \
    t(const t &other) {\
        if(other.size()) {\
            std::memcpy(this, &other, sizeof(*this));\
            auto memsz = other.capacity() * sizeof(*keys);\
            keys = static_cast<decltype(keys)>(std::malloc(memsz));\
            if(!keys) throw std::bad_alloc();\
            std::memcpy(keys, other.keys, memsz);\
            memsz = __ac_fsize(other.capacity()) * sizeof(u32);\
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
    void insert(u32 val) {
        int khr;
        kh_put(set, this, val, &khr);
    }
    template<typename ItType>
    void insert(ItType i1, ItType i2) {
        while(i1 < i2) insert(*i1++);
    }
    void clear() {kh_clear(set, this);}
    bool contains(u32 x) const {return kh_get(set, this, x) != kh_end(this);}
    size_t size() const {return kh_size(static_cast<const khash_t(set) *>(this));}
    size_t capacity() const {return this->n_buckets;}
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
    void insert(u64 val) {
        int khr;
        kh_put(set64, *this, val, &khr);
    }
    void clear() {kh_clear(set64, this);}
    bool contains(u64 x) const {return kh_get(set64, this, x) != kh_end(this);}
    size_t size() const {return kh_size(static_cast<const khash_t(set64) *>(this));}
    size_t capacity() const {return this->n_buckets;}
};
#undef __FE__
#undef COPY_DEC
#undef MOVE_DEC

template<typename T>
size_t capacity(const T &a) {
    return a.capacity();
}

khset64_t&operator+=(khset64_t &a, const khset64_t &b) {
   b.for_each([&](u64 k){a.insert(k);});
   return a;
}

#endif // ndef KHSET_H
