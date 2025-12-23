#include <assert.h>
#include <bsd/stdlib.h>
#include <iostream>
#include <emmintrin.h>   // SSE2
#include <tmmintrin.h>   // SSSE3
#include <cstring>
#include <chrono>        // For timing
#include "dpf.h"

using namespace dpf;

__m128i *proofdb;
__m128i *BB0, *BB1;

__m128i dotProductVector(const __m128i* proofdb, const uint8_t* _t0, size_t nitems, size_t offset = 0) {
    __m128i accumulator = _mm_setzero_si128();

    for (size_t i = 0; i < nitems; ++i) {
        uint8_t mask_byte = -(_t0[i] & 1);  // 0x00 if 0, 0xFF if 1
        __m128i mask = _mm_set1_epi8(static_cast<char>(mask_byte));
        __m128i masked = _mm_and_si128(proofdb[(i + offset) % nitems], mask);
        accumulator = _mm_xor_si128(accumulator, masked);
    }

    return accumulator;
}

int main(int argc, char *argv[])
{
    typedef __m128i leaf_t;

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <log2(nitems)>\n";
        return 1;
    }

    uint64_t log_nitems = std::stoull(argv[1]);
    uint64_t nitems = 1ULL << log_nitems;

    AES_KEY prgkey;

    leaf_t target_value = _mm_set_epi64x(100, 500); // high, low

    proofdb = new leaf_t[nitems];
    std::cout << "nitems = " << nitems << std::endl;
    
    for (size_t j = 0; j < nitems; ++j) {
        uint64_t lo = 100 * j;
        uint64_t hi = 200 * j;
        proofdb[j] = _mm_set_epi64x(hi, lo); // high, low
    }
    
    auto start_P2 = std::chrono::high_resolution_clock::now();
    uint64_t alpha = 534;

    auto [dpfkey0, dpfkey1] = dpf_key<leaf_t, __m128i, AES_KEY>::gen(prgkey, nitems, alpha, target_value);

    uint64_t alpha0, alpha1;
    arc4random_buf(&alpha0, sizeof(alpha0));
    alpha1 = alpha - alpha0;

    auto end_P2 = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double, std::milli> duration_ms_P2 = end_P2 - start_P2;
    std::cout << "P2 time: " << duration_ms_P2.count() << " ms\n";


    BB0 = new leaf_t[nitems]();
    BB1 = new leaf_t[nitems]();

    leaf_t *output0 = new leaf_t[nitems];
    leaf_t *output1 = new leaf_t[nitems];

    uint8_t *_t0 = new uint8_t[nitems];
    uint8_t *_t1 = new uint8_t[nitems];

    uint64_t target_ind = 567;
    uint64_t target_ind0, target_ind1;
    uint64_t proof_submitted0, proof_submitted1;

    arc4random_buf(&proof_submitted0, sizeof(proof_submitted0));
    arc4random_buf(&proof_submitted1, sizeof(proof_submitted1));
    arc4random_buf(&target_ind0, sizeof(target_ind0));
    target_ind1 = target_ind - target_ind0;

    uint64_t offset = (alpha0 - target_ind0) - (alpha1 - target_ind1);

 
    
    // Timing dotProductVector
    auto start_P0_preproc = std::chrono::high_resolution_clock::now();
    
    __evalinterval(0, 0, nitems - 1, output0, _t0);

    auto end_P0_preproc = std::chrono::high_resolution_clock::now();
    
    __m128i x0 = dotProductVector(proofdb, _t0, nitems, offset); // Performed by P0
    

    std::chrono::duration<double, std::milli> duration_ms_P0_preproc = end_P0_preproc - start_P0_preproc;
    std::cout << "P0 time (preprocessing): " << duration_ms_P0_preproc.count() << " ms\n";

    for (size_t j = 0; j + 3 < nitems; j += 4) {
        BB0[j + 0] = _mm_add_epi8(BB0[j + 0], output0[(j + 0 + offset) % nitems]);
        BB0[j + 1] = _mm_add_epi8(BB0[j + 1], output0[(j + 1 + offset) % nitems]);
        BB0[j + 2] = _mm_add_epi8(BB0[j + 2], output0[(j + 2 + offset) % nitems]);
        BB0[j + 3] = _mm_add_epi8(BB0[j + 3], output0[(j + 3 + offset) % nitems]);
    }

 


 
    auto end_P0 = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double, std::milli> duration_ms_total = end_P0 - start_P0_preproc;
    std::cout << "P0 time: " << duration_ms_total.count() << " ms\n";

    std::chrono::duration<double, std::milli> duration_ms_P0_online = end_P0 - end_P0_preproc;
    std::cout << "P0 time (online): " << duration_ms_P0_online.count() << " ms\n";
    
    std::chrono::duration<double, std::milli> duration_ms_P0_total = end_P0 - start_P0_preproc;
    std::cout << "P0 time (total): " << duration_ms_P0_total.count() << " ms\n";


    auto start_P1 = std::chrono::high_resolution_clock::now();
    __evalinterval(dpfkey1, 0, nitems - 1, output1, _t1);

    __m128i x1 = dotProductVector(proofdb, _t1, nitems, offset); // Performed by P1

    __m128i x = _mm_xor_si128(x0, x1);

    std::cout << "x = " << _mm_extract_epi64(x, 0) << " " << _mm_extract_epi64(x, 1) << std::endl;

    __m128i check_proof0 = _mm_xor_si128(x0, _mm_set_epi64x(0, proof_submitted0));
    __m128i check_proof1 = _mm_xor_si128(x1, _mm_set_epi64x(0, proof_submitted1));

    std::cout << "check_proof0 = " << _mm_extract_epi64(check_proof0, 0) << " "
              << _mm_extract_epi64(check_proof0, 1) << std::endl;
    std::cout << "check_proof1 = " << _mm_extract_epi64(check_proof1, 0) << " "
              << _mm_extract_epi64(check_proof1, 1) << std::endl;


    
    for (size_t j = 0; j + 3 < nitems; j += 4) {
        BB1[j + 0] = _mm_add_epi8(BB1[j + 0], output1[(j + 0 + offset) % nitems]);
        BB1[j + 1] = _mm_add_epi8(BB1[j + 1], output1[(j + 1 + offset) % nitems]);
        BB1[j + 2] = _mm_add_epi8(BB1[j + 2], output1[(j + 2 + offset) % nitems]);
        BB1[j + 3] = _mm_add_epi8(BB1[j + 3], output1[(j + 3 + offset) % nitems]);
    }
    auto end_P1 = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double, std::milli> duration_ms_total_P1 = end_P1 - start_P1;
    std::cout << "P1 time: " << duration_ms_total_P1.count() << " ms\n";
    // Cleanup
    delete[] proofdb;
    delete[] BB0;
    delete[] BB1;
    delete[] output0;
    delete[] output1;
    delete[] _t0;
    delete[] _t1;

    return 0;
}
