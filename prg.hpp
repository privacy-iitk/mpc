#pragma once
#include "dpf/prg.h"
#include <vector>
#include <emmintrin.h>   // for __m128i, _mm_cvtsi128_si64, etc.
#include <wmmintrin.h>   // for AES intrinsics
#include <type_traits>

namespace crypto {

// Re-export the low-level AES-PRG
using dpf::PRG;

// ----------------------------------------------------------------------
// Higher-level helper: expand a seed into a vector<T>
// Each AES block (128 bits) produces two 64-bit outputs
// ----------------------------------------------------------------------
template<typename T>
inline void fill_vector_with_prg(std::vector<T>& vec, const AES_KEY& key, const __m128i& seed) {
    static_assert(std::is_integral_v<T>, "fill_vector_with_prg requires an integral type");

    const size_t blocks = (vec.size() + 1) / 2;
    std::vector<__m128i> buf(blocks);

    PRG(key, seed, buf.data(), static_cast<uint32_t>(blocks), 0);

    size_t idx = 0;
    for (size_t i = 0; i < blocks && idx < vec.size(); ++i) {
        uint64_t lo = _mm_cvtsi128_si64(buf[i]);
        uint64_t hi = _mm_extract_epi64(buf[i], 1);
        vec[idx++] = static_cast<T>(lo);
        if (idx < vec.size()) vec[idx++] = static_cast<T>(hi);
    }
}

} // namespace crypto
