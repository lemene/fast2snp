#pragma once


#include "parallel_hashmap/phmap.h"
#include "utils.hpp"


static inline uint64_t hash_64(uint64_t key)
{
  key = (~key + (key << 21));
  key = key ^ key >> 24;
  key = ((key + (key << 3)) + (key << 8));
  key = key ^ key >> 14;
  key = ((key + (key << 2)) + (key << 4));
  key = key ^ key >> 28;
  key = (key + (key << 31));
  return key;
}


class spinlock_mutex {
  public:
     spinlock_mutex& operator=(const spinlock_mutex& rhs) = delete;

     void lock() noexcept {
        for (;;) {
           if (!lock_.exchange(true, std::memory_order_acquire))
              break;
           while (lock_.load(std::memory_order_relaxed))
              __builtin_ia32_pause();
        }
     }

     void unlock() noexcept {
        lock_.store(false, std::memory_order_release);
     }

  private:
     alignas(4 * sizeof(std::max_align_t)) std::atomic_bool lock_ = false;
};

// 重载 unsigned __int128 的输出
inline std::ostream& operator<<(std::ostream& os, unsigned __int128 x) {
    if (x == 0) return os << '0';
    
    char buf[34];  // 128位最多32个十六进制字符 + 结束符
    char* p = buf + 34;
    *--p = '\0';
    
    const char hex_chars[] = "0123456789abcdef";
    
    while (x > 0) {
        *--p = hex_chars[x & 0xF];  // 取低4位
        x >>= 4;                     // 右移4位
    }
    
    return os << p;
}

struct UllHash {
    std::size_t operator()(const KmerInt& key) const {
            // combine low and high 64-bit parts
            uint64_t lo = static_cast<uint64_t>(key);
            uint64_t hi = static_cast<uint64_t>(key >> 64);
            return hash_64(lo) ^ (hash_64(hi) << 1);
    }
};


// GraphIndex definition: Key is Canonical (k-1)-mer, Value packs 8 directional
// edge counts (bytes 0-3: Right A/C/G/T, bytes 4-7: Left A/C/G/T). A zero byte
// means the edge is absent.
using PackedEdgeCounts = unsigned __int128;

static inline int baseToCountSlot(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

static inline uint16_t getPackedCount(PackedEdgeCounts packed, int slot) {
    return static_cast<uint16_t>((packed >> (slot * 16)) & 0xFFFFULL);
}

static inline void setPackedCount(PackedEdgeCounts& packed, int slot, uint16_t count) {
    const PackedEdgeCounts shift = static_cast<PackedEdgeCounts>(slot) * 16ULL;
    const PackedEdgeCounts mask = (PackedEdgeCounts)0xFFFF << shift;
    packed = (packed & ~mask) | (static_cast<PackedEdgeCounts>(count) << shift);
}

static inline int countNonZeroSlots(const PackedEdgeCounts& packed, int startSlot, int endSlot) {
    int count = 0;
    for (int slot = startSlot; slot < endSlot; ++slot) {
        count += getPackedCount(packed, slot) > 0 ? 1 : 0;
    }
    return count;
}

static inline int countValidSlots(PackedEdgeCounts packed, int startSlot, int endSlot, int minCount, int maxCount) {
    int count = 0;
    for (int slot = startSlot; slot < endSlot; ++slot) {
        auto c = getPackedCount(packed, slot);
        if (c >= minCount &&  c <= maxCount) {
            count ++;
        }
    }
    return count;
}

using GraphIndex = phmap::parallel_flat_hash_map<
    KmerInt, PackedEdgeCounts, UllHash,
    std::equal_to<KmerInt>,
    std::allocator<std::pair<const KmerInt, PackedEdgeCounts>>, 10, spinlock_mutex>;