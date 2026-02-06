#include <assert.h>
#include <bsd/stdlib.h>
#include <iostream>

#include <emmintrin.h>   // SSE2
#include <tmmintrin.h>   // SSSE3

#include <cstring>
#include <chrono>
#include <future>
#include <tuple>
#include <array>
#include <utility>
#include <algorithm>

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>

#include "types.hpp"
#include "dpf.hpp"

using namespace dpf;

// -----------------------------------------------------------------------------
// Runs three __evalinterval calls in parallel and waits for completion.
// -----------------------------------------------------------------------------
template<typename KeyT, typename LeafT>
void eval_three_parallel(
    KeyT k0, KeyT k1, KeyT k2,
    size_t lo, size_t hi,
    LeafT *output0, uint8_t *tmp0,
    LeafT *output1, uint8_t *tmp1,
    LeafT *output2, uint8_t *tmp2)
{
    auto fa = std::async(std::launch::async,
        [&k0, lo, hi, output0, tmp0]() {
            __evalinterval(k0, lo, hi, output0, tmp0);
        });

    auto fb = std::async(std::launch::async,
        [&k1, lo, hi, output1, tmp1]() {
            __evalinterval(k1, lo, hi, output1, tmp1);
        });

    auto fc = std::async(std::launch::async,
        [&k2, lo, hi, output2, tmp2]() {
            __evalinterval(k2, lo, hi, output2, tmp2);
        });

    fa.get();
    fb.get();
    fc.get();
}

// -----------------------------------------------------------------------------
// Generate three DPF key pairs in parallel
// -----------------------------------------------------------------------------
template<typename leaf_t, typename node_t, typename prgkey_t>
auto gen_three_parallel(
    const prgkey_t &prgkey,
    std::size_t nitems,
    uint64_t target_ind,
    const leaf_t &target_value)
-> std::array<std::pair<
       dpf_key<leaf_t, node_t, prgkey_t>,
       dpf_key<leaf_t, node_t, prgkey_t>
   >, 3>
{
    using key_t  = dpf_key<leaf_t, node_t, prgkey_t>;
    using pair_t = std::pair<key_t, key_t>;

    std::array<std::future<pair_t>, 3> futs;

    futs[0] = std::async(std::launch::async,
        [prgkey, nitems, target_ind, target_value]() {
            return key_t::gen(prgkey, nitems, target_ind, target_value);
        });

    futs[1] = std::async(std::launch::async,
        [prgkey, nitems, target_ind, target_value]() {
            return key_t::gen(prgkey, nitems, target_ind, target_value);
        });

    futs[2] = std::async(std::launch::async,
        [prgkey, nitems, target_ind, target_value]() {
            return key_t::gen(prgkey, nitems, target_ind, target_value);
        });

    pair_t r0 = futs[0].get();
    pair_t r1 = futs[1].get();
    pair_t r2 = futs[2].get();

    return { std::move(r0), std::move(r1), std::move(r2) };
}

// -----------------------------------------------------------------------------
// Build an __m128i with exactly one lane set
// -----------------------------------------------------------------------------
static inline __m128i build_m128i_single_lane(int lane_bytes, int k, uint64_t value)
{
    uint8_t tmp[16];
    std::memset(tmp, 0, sizeof tmp);

    if (lane_bytes == 1) {
        assert(k < 16);
        tmp[k] = static_cast<uint8_t>(value);
    }
    else if (lane_bytes == 4) {
        assert(k < 4);
        uint32_t v32 = static_cast<uint32_t>(value);
        std::memcpy(tmp + (k * 4), &v32, sizeof(v32));
    }
    else if (lane_bytes == 8) {
        assert(k < 2);
        uint64_t v64 = value;
        std::memcpy(tmp + (k * 8), &v64, sizeof(v64));
    }
    else {
        assert(!"lane_bytes must be 1, 4, or 8");
    }

    return _mm_loadu_si128(reinterpret_cast<const __m128i*>(tmp));
}

void print_m128i_all(__m128i v, const char *label)
{
    alignas(16) uint8_t bytes[16];
    _mm_storeu_si128((__m128i*)bytes, v);

    const uint64_t *u64 = (const uint64_t*)bytes;
    const uint32_t *u32 = (const uint32_t*)bytes;

    printf("=== %s ===\n", label);
    printf("u64: [%016" PRIx64 ", %016" PRIx64 "]\n", u64[0], u64[1]);
    printf("u32: [%08" PRIx32 ", %08" PRIx32 ", %08" PRIx32 ", %08" PRIx32 "]\n",
           u32[0], u32[1], u32[2], u32[3]);

    printf("u8 :");
    for (int i = 0; i < 16; i++)
        printf(" %02" PRIx8, bytes[i]);
    printf("\n");
}



#pragma once
#include <type_traits>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <emmintrin.h>

// ============================================================================
// Low-level printer for __m128i
// ============================================================================
inline void print_m128i(const __m128i &v)
{
    alignas(16) uint8_t bytes[16];
    _mm_storeu_si128((__m128i*)bytes, v);

    printf("__m128i: [");
    for (int i = 0; i < 16; i++) {
        printf("%02x", bytes[i]);
        if (i < 15) printf(" ");
    }
    printf("]");
}

// ============================================================================
// Generic print for scalar integral types
// ============================================================================
template<typename T>
inline void print_scalar(const T &x)
{
    if constexpr (std::is_same_v<T, uint8_t>) {
        printf("%u", (unsigned)x);
    }
    else if constexpr (std::is_same_v<T, uint32_t>) {
        printf("%u", x);
    }
    else if constexpr (std::is_same_v<T, uint64_t>) {
        printf("%llu", (unsigned long long)x);
    }
    else {
        static_assert(sizeof(T) == 0, "Unsupported scalar type in print_scalar");
    }
}

// ============================================================================
// Master leaf printer: works for all required types
// ============================================================================
template<typename T>
void print_leaf(const T &x)
{
    // ------------------------------------------------------------------------
    // Case 1: __m128i
    // ------------------------------------------------------------------------
    if constexpr (std::is_same_v<T, __m128i>) {
        print_m128i(x);
    }

    // ------------------------------------------------------------------------
    // Case 2: scalar integrals uint8_t / uint32_t / uint64_t
    // ------------------------------------------------------------------------
    else if constexpr (std::is_integral_v<T>) {
        print_scalar(x);
    }

    // ------------------------------------------------------------------------
    // Case 3: std::array<__m128i, K>
    // ------------------------------------------------------------------------
    else if constexpr (
        std::is_same_v<typename T::value_type, __m128i> &&
        std::is_array_v<T> == false &&
        std::tuple_size<T>::value > 0
    ) {
        printf("[ ");
        for (size_t i = 0; i < x.size(); i++) {
            print_m128i(x[i]);
            if (i + 1 < x.size()) printf(", ");
        }
        printf(" ]");
    }

    // ------------------------------------------------------------------------
    // Case 4: std::array<U, K> for scalar U
    // ------------------------------------------------------------------------
    else if constexpr (
        std::tuple_size<T>::value > 0 &&
        std::is_integral_v<typename T::value_type>
    ) {
        printf("[ ");
        for (size_t i = 0; i < x.size(); i++) {
            print_scalar(x[i]);
            if (i + 1 < x.size()) printf(", ");
        }
        printf(" ]");
    }

    // ------------------------------------------------------------------------
    // Unsupported
    // ------------------------------------------------------------------------
    else {
        static_assert(sizeof(T) == 0,
                      "print_leaf(): unsupported leaf type");
    }
}




// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    AES_KEY prgkey;

    uint64_t target_ind = 155;

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <log2(nitems)>\n";
        return 1;
    }

    uint64_t log_nitems = std::stoull(argv[1]);
    uint64_t nitems = 1ULL << log_nitems;

    

    using leaf_t =  __m128i; //std::array<__m128i, 2>;
    leaf_t target_value;
    set_target_values(target_value, 100, 300);
 
    leaf_t *output0 = new leaf_t[nitems];
    leaf_t *output1 = new leaf_t[nitems];
    uint8_t *_t0 = new uint8_t[nitems];
    uint8_t *_t1 = new uint8_t[nitems];

    auto start = std::chrono::high_resolution_clock::now();
    auto [dpfkey0, dpfkey1] =
        dpf_key<leaf_t, __m128i, AES_KEY>::gen(prgkey, nitems, target_ind, target_value);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Gen time = "
              << std::chrono::duration<double, std::milli>(end - start).count()
              << " ms\n\n";

    auto start_parallel = std::chrono::high_resolution_clock::now();
    auto keys = gen_three_parallel<leaf_t, __m128i, AES_KEY>(
                    prgkey, nitems, target_ind, target_value);
    auto end_parallel = std::chrono::high_resolution_clock::now();

    std::cout << "Gen time_parallel = "
              << std::chrono::duration<double, std::milli>(end_parallel - start_parallel).count()
              << " ms\n\n";

    auto k0_1 = std::move(keys[0].first);
    auto k1_1 = std::move(keys[0].second);
    auto k0_2 = std::move(keys[1].first);
    auto k1_2 = std::move(keys[1].second);
    auto k0_3 = std::move(keys[2].first);
    auto k1_3 = std::move(keys[2].second);

    // Three parallel eval buffers
    leaf_t *output00 = new leaf_t[nitems];
    leaf_t *output01 = new leaf_t[nitems];
    leaf_t *output02 = new leaf_t[nitems];

    uint8_t *tmp0 = new uint8_t[nitems];
    uint8_t *tmp1 = new uint8_t[nitems];
    uint8_t *tmp2 = new uint8_t[nitems];

    auto t0 = std::chrono::high_resolution_clock::now();
    eval_three_parallel(
        std::move(k0_1), std::move(k1_1), std::move(k0_3),
        0, nitems - 1,
        output00, tmp0,
        output01, tmp1,
        output02, tmp2);
    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "Three evalintervals in parallel took "
              << std::chrono::duration<double, std::milli>(t1 - t0).count()
              << " ms\n";

    auto start_eval = std::chrono::high_resolution_clock::now();
    __evalinterval_mpc(dpfkey0, 0, nitems - 1, output0, _t0);
    auto end_eval = std::chrono::high_resolution_clock::now();

    std::cout << "Evalfull time = "
              << std::chrono::duration<double, std::milli>(end_eval - start_eval).count()
              << " ms\n\n";

    std::cout << "\n\n ----------------------------------------------------------- \n\n";

    __evalinterval_mpc(dpfkey1, 0, nitems - 1, output1, _t1);

    for (size_t j = 0; j < nitems; ++j) {
        if (!leaf_equal(output0[j], output1[j])) {
            std::cout << j << "(flags): → "
                      << (int)_t0[j] << " <> " << (int)_t1[j] << "\n";

            std::cout << j << "(flags): → "  << std::endl;
             
             print_leaf(leaf_xor(output0[j], output1[j]));
             std::cout << "\n\n\n";
             print_leaf(output0[j]);
             std::cout << " <> "; 
             print_leaf(output1[j]);
             std::cout << std::endl;
        }
    }

    delete[] output0;
    delete[] output1;
    delete[] _t0;
    delete[] _t1;

    delete[] output00;
    delete[] output01;
    delete[] output02;
    delete[] tmp0;
    delete[] tmp1;
    delete[] tmp2;

    return 0;
}
