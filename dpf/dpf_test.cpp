#include <assert.h>
#include <bsd/stdlib.h>
#include <iostream>
#include <emmintrin.h>   // SSE2
#include <tmmintrin.h>   // SSSE3
#include <cstring>
#include <chrono>        // For timing
#include "dpf.hpp"
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <algorithm>     // for std::min
#include <future>
#include <array>
#include <utility>

using namespace dpf;

// helper traits to detect std::array<__m128i, N>
template<typename T>
struct is_array_of_m128i : std::false_type { static constexpr size_t N = 0; };

template<size_t N>
struct is_array_of_m128i<std::array<__m128i, N>> : std::true_type {
    static constexpr size_t N_v = N;
};

// Build a __m128i with two uint64_t words (high, low)
static inline __m128i make_m128i_from_u64pair(uint64_t hi, uint64_t lo) {
    return _mm_set_epi64x((long long)hi, (long long)lo);
}

// Generic maker: construct a single leaf_t value from two 64-bit words (hi, lo).
template<typename LeafT>
LeafT make_leaf_from_u64_pair(uint64_t hi, uint64_t lo) {
    if constexpr (std::is_same_v<LeafT, uint8_t>) {
        // keep lowest byte
        return static_cast<uint8_t>(lo & 0xFFu);
    } else if constexpr (std::is_same_v<LeafT, uint32_t>) {
        return static_cast<uint32_t>(lo & 0xFFFFFFFFu);
    } else if constexpr (std::is_same_v<LeafT, uint64_t>) {
        return static_cast<uint64_t>(lo);
    } else if constexpr (std::is_same_v<LeafT, __m128i>) {
        return make_m128i_from_u64pair(hi, lo);
    } else if constexpr (is_array_of_m128i<LeafT>::value) {
        // For an array<__m128i, N> build an array where element 0 holds (hi, lo)
        // and the remainder are zero.
        LeafT out;
        const size_t N = is_array_of_m128i<LeafT>::N_v;
        // zero-initialize
        for (size_t i = 0; i < N; ++i) out[i] = _mm_setzero_si128();
        out[0] = make_m128i_from_u64pair(hi, lo);
        return out;
    } else {
        static_assert(!std::is_same_v<LeafT, LeafT>, "Unsupported LeafT in make_leaf_from_u64_pair");
    }
}

// Convenience: set two target slots (index 0 and 1) like your original code.
// If Leaf is not an array, it sets a single value (index 0).
template<typename LeafT>
void set_target_values(LeafT &target_value,
                       uint64_t hi0, uint64_t lo0,
                       uint64_t hi1 = 0, uint64_t lo1 = 0)
{
    if constexpr (is_array_of_m128i<LeafT>::value) {
        // std::array<__m128i, N> case: fill indices 0 and 1 (if available)
        constexpr size_t N = is_array_of_m128i<LeafT>::N_v;
        target_value[0] = make_m128i_from_u64pair(hi0, lo0);
        if (N > 1) target_value[1] = make_m128i_from_u64pair(hi1, lo1);
        // zero rest (optional)
        for (size_t i = 2; i < N; ++i) target_value[i] = _mm_setzero_si128();
    } else {
        // scalar or single __m128i: set one value (use hi0/lo0)
        target_value = make_leaf_from_u64_pair<LeafT>(hi0, lo0);
    }
}



// -----------------------------------------------------------------------------
// Helper: generate three DPF key-pairs in parallel (or throw on error).
// Returns std::array of 3 std::pair<key_t,key_t> where each pair is (k0,k1).
// Copies prgkey into each task to avoid potential non-reentrant key-schedule use.
// -----------------------------------------------------------------------------
template<typename leaf_t, typename node_t, typename prgkey_t>
auto gen_three_parallel(const prgkey_t &prgkey,
                        std::size_t nitems,
                        uint64_t target_ind,
                        const leaf_t &target_value)
-> std::array<std::pair< dpf_key<leaf_t, node_t, prgkey_t>,
                        dpf_key<leaf_t, node_t, prgkey_t> >, 3>
{
    using key_t  = dpf_key<leaf_t, node_t, prgkey_t>;
    using pair_t = std::pair<key_t, key_t>;

    std::array<std::future<pair_t>, 3> futs;

    // Launch three async tasks; capture prgkey and target_value by value.
    futs[0] = std::async(std::launch::async,
                         [prgkey, nitems, target_ind, target_value]() -> pair_t {
                             return key_t::gen(prgkey, nitems, target_ind, target_value);
                         });

    futs[1] = std::async(std::launch::async,
                         [prgkey, nitems, target_ind, target_value]() -> pair_t {
                             return key_t::gen(prgkey, nitems, target_ind, target_value);
                         });

    futs[2] = std::async(std::launch::async,
                         [prgkey, nitems, target_ind, target_value]() -> pair_t {
                             return key_t::gen(prgkey, nitems, target_ind, target_value);
                         });

    // Retrieve futures into local temporaries (moveable)
    pair_t r0 = futs[0].get();   //moves or constructs into r0
    pair_t r1 = futs[1].get();
    pair_t r2 = futs[2].get();

    // Construct the array from the three temporaries (use move to be explicit)
    return std::array<pair_t, 3>{ std::move(r0), std::move(r1), std::move(r2) };
}


// Build an __m128i with exactly one lane set and all other bytes zero.
// lane_bytes: 1, 4, or 8
// k: lane index
// value: placed into the lane
static inline __m128i build_m128i_single_lane(int lane_bytes, int k, uint64_t value)
{
    uint8_t tmp[16];
    std::memset(tmp, 0, sizeof tmp);

    if (lane_bytes == 1) {
        if ((unsigned)k >= 16u) {
            assert(!"k out of range for lane_bytes==1");
            return _mm_setzero_si128();
        }
        tmp[k] = static_cast<uint8_t>(value & 0xFFu);
    }
    else if (lane_bytes == 4) {
        if ((unsigned)k >= 4u) {
            assert(!"k out of range for lane_bytes==4");
            return _mm_setzero_si128();
        }
        uint32_t v32 = static_cast<uint32_t>(value & 0xFFFFFFFFu);
        std::memcpy(tmp + (k * 4), &v32, sizeof(v32));
    }
    else if (lane_bytes == 8) {
        if ((unsigned)k >= 2u) {
            assert(!"k out of range for lane_bytes==8");
            return _mm_setzero_si128();
        }
        uint64_t v64 = value;
        std::memcpy(tmp + (k * 8), &v64, sizeof(v64));
    }
    else {
        assert(!"lane_bytes must be 1, 4, or 8");
        return _mm_setzero_si128();
    }

    return _mm_loadu_si128(reinterpret_cast<const __m128i*>(tmp));
}

void print_m128i(__m128i v, const char *label)
{
    alignas(16) uint8_t bytes[16];
    _mm_storeu_si128((__m128i*)bytes, v);

    const uint64_t *u64 = (const uint64_t*)bytes;
    const uint32_t *u32 = (const uint32_t*)bytes;

    printf("=== %s ===\n", label);
    printf("u64: x[0]=0x%016" PRIx64 ", x[1]=0x%016" PRIx64 "\n",
           u64[0], u64[1]);

    printf("u32: x[0]=0x%08" PRIx32 ", x[1]=0x%08" PRIx32
           ", x[2]=0x%08" PRIx32 ", x[3]=0x%08" PRIx32 "\n",
           u32[0], u32[1], u32[2], u32[3]);

    printf("u8 :");
    for (int i = 0; i < 16; i++) {
        printf(" %02" PRIx8, bytes[i]);
    }
    printf("\n");
}

int main(int argc, char *argv[])
{
    
    uint64_t target_ind = 11;

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <log2(nitems)>\n";
        return 1;
    }

    uint64_t log_nitems = std::stoull(argv[1]);
    uint64_t nitems = 1ULL << log_nitems;

    AES_KEY prgkey;

    typedef std::array<__m128i, 6> leaf_t;
    leaf_t target_value;
  
    target_value[0] = _mm_set_epi64x(100, 300);
    target_value[1] = _mm_set_epi64x(200, 400);

    leaf_t *output0 = new leaf_t[nitems];
    leaf_t *output1 = new leaf_t[nitems];

    uint8_t *_t0 = new uint8_t[nitems];
    uint8_t *_t1 = new uint8_t[nitems];

    auto start = std::chrono::high_resolution_clock::now();

    auto [dpfkey0, dpfkey1] =
        dpf_key<leaf_t, __m128i, AES_KEY>::gen(prgkey, nitems, target_ind, target_value);

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed_us = std::chrono::duration<double, std::milli>(end - start).count();
    std::cout << "Gen time = " << elapsed_us << " ms\n\n\n\n";

    // -----------------------------------------------------------------
    // Use the helper to generate three DPF key-pairs in parallel
    // -----------------------------------------------------------------
    auto start_parallel = std::chrono::high_resolution_clock::now();

    auto keys = gen_three_parallel<leaf_t, __m128i, AES_KEY>(prgkey, nitems, target_ind, target_value);

    // extract k0_1, k1_1, k0_2, k1_2, k0_3, k1_3 for later use
    auto k0_1 = std::move(keys[0].first);
    auto k1_1 = std::move(keys[0].second);
    auto k0_2 = std::move(keys[1].first);
    auto k1_2 = std::move(keys[1].second);
    auto k0_3 = std::move(keys[2].first);
    auto k1_3 = std::move(keys[2].second);

    auto end_parallel = std::chrono::high_resolution_clock::now();
    double elapsed_us_parallel =
        std::chrono::duration<double, std::milli>(end_parallel - start_parallel).count();

    std::cout << "Gen time_parallel = " << elapsed_us_parallel << " ms\n\n\n\n";

    // prepare three distinct output buffers and temps
    leaf_t *output00 = new leaf_t[nitems];
    leaf_t *output01 = new leaf_t[nitems];
    leaf_t *output02 = new leaf_t[nitems];

    uint8_t *tmp0 = new uint8_t[nitems];
    uint8_t *tmp1 = new uint8_t[nitems];
    uint8_t *tmp2 = new uint8_t[nitems];

    auto t0 = std::chrono::high_resolution_clock::now();

    auto fa = std::async(std::launch::async, [&] {
        __evalinterval(k0_1, 0, nitems - 1, output00, tmp0);
    });
    auto fb = std::async(std::launch::async, [&] {
        __evalinterval(k1_1, 0, nitems - 1, output01, tmp1);
    });
    auto fc = std::async(std::launch::async, [&] {
        __evalinterval(k0_3, 0, nitems - 1, output02, tmp2);
    });

    fa.get();
    fb.get();
    fc.get();

    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed_ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::cout << "Three evalintervals in parallel took " << elapsed_ms << " ms\n";

    auto start_eval = std::chrono::high_resolution_clock::now();
    __evalinterval(dpfkey0, 0, nitems - 1, output0, _t0);
    auto end_eval = std::chrono::high_resolution_clock::now();

    double elapsed_us_eval =
        std::chrono::duration<double, std::milli>(end_eval - start_eval).count();
    std::cout << "Evalfull time = " << elapsed_us_eval << " ms\n\n\n\n\n";

    std::cout << "\n\n ----------------------------------------------------------- \n\n";

    __evalinterval(dpfkey1, 0, nitems - 1, output1, _t1);

    for (size_t j = 0; j < nitems; ++j) {
        if (output00[j][0][0] != output01[j][0][0]) {
            std::cout << j << " :--> "
                      << (output0[j][0][0] ^ output1[j][0][0]) << " <---> "
                      << (output0[j][1][0] ^ output1[j][1][0]) << std::endl;

            std::cout << j << " :--> "
                      << (output0[j][0][1] ^ output1[j][0][1]) << " <---> "
                      << (output0[j][1][1] ^ output1[j][1][1]) << std::endl;

            std::cout << j << "(flags ): -> "
                      << (int)_t0[j] << " <> " << (int)_t1[j] << std::endl;
        }
    }

#ifdef VERBOSE
    for (size_t j = 0; j < nitems; ++j) {
        std::cout << j << "(flags ): -> "
                  << (int)_t0[j] << " <> " << (int)_t1[j] << std::endl;
    }

    __m128i a = build_m128i_single_lane(1, 15, 0xFF);
    __m128i b = build_m128i_single_lane(4, 2, 0xDEADBEEF);
    __m128i c = build_m128i_single_lane(8, 1, 0x0123456789ABCDEFULL);

    std::cout << "\n\n ----------------------------------------------------------- \n\n";
    print_m128i(a, "uint8_t");
    std::cout << "\n\n ----------------------------------------------------------- \n\n";
    print_m128i(b, "uint32_t");
    std::cout << "\n\n ----------------------------------------------------------- \n\n";
    print_m128i(c, "uint64_t");
    std::cout << "\n\n ----------------------------------------------------------- \n\n";
#endif

    delete[] output0;
    delete[] output1;
    delete[] _t0;
    delete[] _t1;

    // cleanup for parallel buffers
    delete[] output00;
    delete[] output01;
    delete[] output02;
    delete[] tmp0;
    delete[] tmp1;
    delete[] tmp2;

    return 0;
}
