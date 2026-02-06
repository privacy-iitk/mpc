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



 
#include <type_traits>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <emmintrin.h>




awaitable<void> run_party(
    boost::asio::io_context& io,
    Role role,
    const size_t log_nitems
) {
    using leaf_t = __m128i;
    
    const size_t nitems = 1ULL << log_nitems;
    const size_t target_ind = 44;

    __m128i* DB = (__m128i*)aligned_alloc(16, nitems * sizeof(__m128i));
    arc4random_buf(DB, nitems * sizeof(__m128i));

    __m128i* ProofDB = (__m128i*)aligned_alloc(16, nitems * sizeof(__m128i));
    arc4random_buf(ProofDB, nitems * sizeof(__m128i));

    // std::array<uint8_t, log_nitems> a{};
    // std::array<uint8_t, log_nitems> b{};
    uint8_t* a = static_cast<uint8_t*>(malloc(log_nitems));
uint8_t* b = static_cast<uint8_t*>(malloc(log_nitems));

if (!a || !b) {
    throw std::bad_alloc();
}

    // fresh randomness for masking
    std::mt19937_64 gen2(999);   // fixed seed for reproducibility
    std::uniform_int_distribution<int> bitdist(0, 1);

    for (size_t i = 0; i < log_nitems; ++i) {
        size_t bitpos = log_nitems - 1 - i;   // MSB first
        uint8_t bit = (target_ind >> bitpos) & 1;
        a[i] = bitdist(gen2);
        b[i] = a[i] ^ bit;
    }

    size_t reconstructed = 0;
    for (size_t i = 0; i < log_nitems; ++i) {
        reconstructed |=
            (size_t(a[i] ^ b[i]) << (log_nitems - 1 - i));
    }


    assert(reconstructed == target_ind);



        AES_KEY prgkey;
        leaf_t target_value;
        set_target_values(target_value, 100, 300);

        auto [k0, k1] =
            dpf_key<leaf_t, __m128i, AES_KEY>::gen(
                prgkey, nitems, target_ind, target_value);
    // =================================================
    // P0
    // =================================================
    if (role == Role::P0) {
        auto s_peer = co_await make_server(io, 9200);
        NetPeer peer(Role::P1, std::move(s_peer));

        MPCContext ctx(Role::P0, peer, &peer);

        // dummy sync
        co_await peer.send(uint8_t{1});

        leaf_t*  _output = new leaf_t[nitems];
        uint8_t* _t      = new uint8_t[nitems];

        /* ================= Eval MPC timing ================= */
        auto start_eval = std::chrono::high_resolution_clock::now();

        co_await __evalinterval_mpc(
            ctx, peer,
            k0,
            0, nitems - 1,
            _output, _t, a, 8
        );

        auto end_eval = std::chrono::high_resolution_clock::now();
        auto eval_ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(end_eval - start_eval).count();

        std::cout << "Time for evalinterval_mpc: " << eval_ms << " ms\n";
        /* =================================================== */


        /* ================= Audit timing ================= */
        auto start_audit = std::chrono::high_resolution_clock::now();

        __m128i acc0 = _mm_setzero_si128();
        __m128i acc1 = _mm_setzero_si128();
        __m128i acc2 = _mm_setzero_si128();
        __m128i acc3 = _mm_setzero_si128();

        size_t j = 0;
        for (; j + 3 < nitems; j += 4) {
            if (_t[j+0]) acc0 = _mm_xor_si128(acc0, ProofDB[j+0]);
            if (_t[j+1]) acc1 = _mm_xor_si128(acc1, ProofDB[j+1]);
            if (_t[j+2]) acc2 = _mm_xor_si128(acc2, ProofDB[j+2]);
            if (_t[j+3]) acc3 = _mm_xor_si128(acc3, ProofDB[j+3]);
        }

        for (; j < nitems; ++j) {
            if (_t[j]) acc0 = _mm_xor_si128(acc0, ProofDB[j]);
        }

        __m128i result =
            _mm_xor_si128(
                _mm_xor_si128(acc0, acc1),
                _mm_xor_si128(acc2, acc3)
            );


        auto end_audit = std::chrono::high_resolution_clock::now();
        auto audit_ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(end_audit - start_audit).count();

        // prevent DCE
        asm volatile("" :: "x"(result));

        std::cout << "Time to audit: " << audit_ms << " ms\n";
        /* ================================================= */


        /* ================= Write timing ================= */
        auto start_write = std::chrono::high_resolution_clock::now();

        for (size_t j = 0; j < nitems; ++j) {
            DB[j] ^= _output[j];
        }

        auto end_write = std::chrono::high_resolution_clock::now();
        auto write_ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(end_write - start_write).count();

        std::cout << "Time to write: " << write_ms << " ms\n";
        /* ================================================= */

        delete[] _output;
        delete[] _t;

        co_return;
    }


    // =================================================
    // P1
    // =================================================
    if (role == Role::P1) {
        auto s_peer = co_await connect_with_retry(io, "127.0.0.1", 9200);
        NetPeer peer(Role::P0, std::move(s_peer));

        MPCContext ctx(Role::P1, peer, &peer);
 
        // dummy sync
        auto dummy = co_await peer.recv<uint8_t>();
        std::cout << "[P1] dummy = " << int(dummy) << std::endl;
    
        leaf_t *_output = new leaf_t[nitems];
        uint8_t *_t = new uint8_t[nitems];
                auto start_eval = std::chrono::high_resolution_clock::now();
        co_await __evalinterval_mpc(
            ctx, peer,
            k1,
            0, nitems - 1,
            _output, _t, b, 8
        );

         auto end_eval = std::chrono::high_resolution_clock::now();
        auto eval_ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(end_eval - start_eval).count();

        std::cout << "Time for evalinterval_mpc: " << eval_ms << " ms\n";

        auto start = std::chrono::high_resolution_clock::now();

        __m128i acc0 = _mm_setzero_si128();
        __m128i acc1 = _mm_setzero_si128();
        __m128i acc2 = _mm_setzero_si128();
        __m128i acc3 = _mm_setzero_si128();

        size_t j = 0;
        for (; j + 3 < nitems; j += 4) {
            if (_t[j+0]) acc0 = _mm_xor_si128(acc0, ProofDB[j+0]);
            if (_t[j+1]) acc1 = _mm_xor_si128(acc1, ProofDB[j+1]);
            if (_t[j+2]) acc2 = _mm_xor_si128(acc2, ProofDB[j+2]);
            if (_t[j+3]) acc3 = _mm_xor_si128(acc3, ProofDB[j+3]);
        }

        for (; j < nitems; ++j) {
            if (_t[j]) acc0 = _mm_xor_si128(acc0, ProofDB[j]);
        }

        __m128i result =
            _mm_xor_si128(
                _mm_xor_si128(acc0, acc1),
                _mm_xor_si128(acc2, acc3)
            );

        auto end = std::chrono::high_resolution_clock::now();
        auto ns = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        std::cout << "Time: " << ns << " ms\n";

                auto start2 = std::chrono::high_resolution_clock::now();

     
        for (size_t j = 0; j < nitems; ++j) {
            DB[j] ^= _output[j];
        }

        auto end2 = std::chrono::high_resolution_clock::now();
        auto ns2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();

        std::cout << "Time to write: " << ns2 << " ms\n";

        // for(size_t jj = 0; jj < nitems; ++jj)
        // {
        //     leaf_t rec_out;
        //     co_await (peer << _output[jj]) ;
        //     co_await (peer >> rec_out);
        //     rec_out ^= _output[jj];;
        //     if(rec_out[0] != 0 || rec_out[1] !=0)
        //        std::cout << jj << " -> " << _output[jj] << std::endl;

        // }

        co_return;
    }
}



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

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "usage: ./dpf_mpc_test [p0|p1] <num_items>\n";
        return 1;
    }

    Role role;
    std::string r = argv[1];

    size_t log_nitems = 0;
    try {
        log_nitems = std::stoull(argv[2]);
    } catch (const std::exception& e) {
        std::cerr << "Invalid log_nitems: " << argv[2] << "\n";
        return 1;
    }

    // optional sanity check
    if (log_nitems == 0) {
        std::cerr << "log_nitems must be > 0\n";
        return 1;
    }

 if (r == "p0")
        role = Role::P0;
    else if (r == "p1")
        role = Role::P1;
    else {
        std::cerr << "invalid role (use p0 or p1)\n";
        return 1;
    }

    boost::asio::io_context io;

    // -------------------------------------------------
    // SAME function, SAME code path
    // -------------------------------------------------
    boost::asio::co_spawn(
        io,
        run_party(io, role, log_nitems),
        boost::asio::detached
    );

    io.run();
    return 0;
}



// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------
// int main(int argc, char *argv[])
// {
//     AES_KEY prgkey;

//     uint64_t target_ind = 155;

//     if (argc != 2) {
//         std::cerr << "Usage: " << argv[0] << " <log2(nitems)>\n";
//         return 1;
//     }

//     uint64_t log_nitems = std::stoull(argv[1]);
//     uint64_t nitems = 1ULL << log_nitems;

    

//     using leaf_t =  __m128i; //std::array<__m128i, 2>;
//     leaf_t target_value;
//     set_target_values(target_value, 100, 300);
 
//     leaf_t *output0 = new leaf_t[nitems];
//     leaf_t *output1 = new leaf_t[nitems];
//     uint8_t *_t0 = new uint8_t[nitems];
//     uint8_t *_t1 = new uint8_t[nitems];

//     auto start = std::chrono::high_resolution_clock::now();
    
//     auto [dpfkey0, dpfkey1] =
//         dpf_key<leaf_t, __m128i, AES_KEY>::gen(prgkey, nitems, target_ind, target_value);
    
//     auto end = std::chrono::high_resolution_clock::now();
//     std::cout << "Gen time = "
//               << std::chrono::duration<double, std::milli>(end - start).count()
//               << " ms\n\n";



//     __evalinterval_mpc(dpfkey1, 0, nitems - 1, output1, _t1);

//     for (size_t j = 0; j < nitems; ++j) {
//         if (!leaf_equal(output0[j], output1[j])) {
//             std::cout << j << "(flags): → "
//                       << (int)_t0[j] << " <> " << (int)_t1[j] << "\n";

//             std::cout << j << "(flags): → "  << std::endl;
             
//              print_leaf(leaf_xor(output0[j], output1[j]));
//              std::cout << "\n\n\n";
//              print_leaf(output0[j]);
//              std::cout << " <> "; 
//              print_leaf(output1[j]);
//              std::cout << std::endl;
//         }
//     }

// //    delete[] output0;
//     delete[] output1;
//   //  delete[] _t0;
//     delete[] _t1;

//     // delete[] output00;
//     // delete[] output01;
//     // delete[] output02;
//     // delete[] tmp0;
//     // delete[] tmp1;
//     // delete[] tmp2;

//     return 0;
// }
