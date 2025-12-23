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

uint64_t dotProductVector(const uint64_t* X, const uint64_t* Y, size_t nitems)
{
  uint64_t dot = 0;

  for(size_t j = 0; j < nitems; ++j)
  {
    dot += X[j] * Y[j]; 
  }
  return dot;
}


uint64_t *X0, *Y0, *X1, *Y1, alpha; 
uint64_t *DB0, *DB1, *DB0_tilde, *DB1_tilde;
uint64_t *e0, *e1, *e0_tilde, *e1_tilde, *e0_plus, *e1_plus, *f0, *f1, * neg_e0, * neg_e1;

uint64_t * X0_new;



void refresh_blinds(uint64_t* DB0,  uint64_t* DB0_tilde, uint64_t * X0,  uint64_t * update_vector0, uint64_t * update_vector1, size_t nitems)
{
    for(size_t j = 0; j < nitems; ++j)
    {
        DB0[j] = DB0[j] + update_vector0[j]; // Dpne by P0
        X0[j]  = X0[j] - update_vector0[j]; // Done by P0

        DB0_tilde[j] = DB0_tilde[j] - update_vector1[j] + update_vector1[j]; // Done by P1

    }
}


int main(int argc, char *argv[])
{
    typedef __m128i leaf_t;
    uint64_t target_ind = 5;
    uint64_t target_ind2 = 5;
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <log2(nitems)>\n";
        return 1;
    }



    uint64_t log_nitems = std::stoull(argv[1]);
    uint64_t nitems = 1ULL << log_nitems;

    X0 = new uint64_t[nitems]; // Blinds of the database
    X1 = new uint64_t[nitems];  // Blinds of the database
    Y0 = new uint64_t[nitems];
    Y1 = new uint64_t[nitems]; 
    DB0 = new uint64_t[nitems];
    DB1 = new uint64_t[nitems];
    DB0_tilde = new uint64_t[nitems];
    DB1_tilde = new uint64_t[nitems];
    e0_plus = new uint64_t[nitems];
    e1_plus = new uint64_t[nitems];
    e0 = new uint64_t[nitems];
    e1 = new uint64_t[nitems];
    neg_e0 = new uint64_t[nitems];
    neg_e1 = new uint64_t[nitems];
    f0 = new uint64_t[nitems];
    f1 = new uint64_t[nitems];
    e0_tilde = new uint64_t[nitems]; 
    e1_tilde = new uint64_t[nitems];

    for(size_t j = 0; j < nitems; ++j)
    {
        arc4random_buf(&X0[j], sizeof(uint64_t));
        arc4random_buf(&X1[j], sizeof(uint64_t));
        arc4random_buf(&Y0[j], sizeof(uint64_t));
        arc4random_buf(&Y1[j], sizeof(uint64_t));
        
    }

    for(size_t j = 0; j < nitems; ++j)
    {
        DB0[j] = 0;
        DB1[j] = j;
        DB0_tilde[j] = DB0[j] + X0[j];
        DB1_tilde[j] = DB1[j] + X1[j];

    }

    AES_KEY prgkey;

    leaf_t target_value = _mm_set_epi64x(100, 500); // high, low

     uint64_t alpha = 534;
    leaf_t *output0 = new leaf_t[nitems];
    leaf_t *output1 = new leaf_t[nitems];

    uint8_t *_t0 = new uint8_t[nitems];
    uint8_t *_t1 = new uint8_t[nitems];

     
    auto [dpfkey0, dpfkey1]       = dpf_key<leaf_t, __m128i, AES_KEY>::gen(prgkey, nitems, target_ind, target_value);
    auto [a_dpfkey0, a_dpfkey1]   = dpf_key<leaf_t, __m128i, AES_KEY>::gen(prgkey, nitems, target_ind, target_value);
    auto [b_dpfkey0, b_dpfkey1]   = dpf_key<leaf_t, __m128i, AES_KEY>::gen(prgkey, nitems, target_ind, target_value);
    
    auto [c_dpfkey0, c_dpfkey1]   = dpf_key<leaf_t, __m128i, AES_KEY>::gen(prgkey, nitems, alpha, target_value);
    auto [d_dpfkey0, d_dpfkey1]   = dpf_key<leaf_t, __m128i, AES_KEY>::gen(prgkey, nitems, alpha, target_value);
    auto [e_dpfkey0, e_dpfkey1]   = dpf_key<leaf_t, __m128i, AES_KEY>::gen(prgkey, nitems, alpha, target_value);
    auto [f_dpfkey0, f_dpfkey1]   = dpf_key<leaf_t, __m128i, AES_KEY>::gen(prgkey, nitems, alpha, target_value);

    __evalinterval(dpfkey0, 0, nitems - 1, output0, _t0);
    __evalinterval(dpfkey1, 0, nitems - 1, output1, _t1);
    

    for(size_t j = 0; j < nitems; ++j)
    {   
        e0[j] = _t0[j];
        e1[j] = -_t1[j];
        
        if(j == target_ind) 
        {
            if (e0[j] == 0) e1[j] = 1;
            if (e1[j] == 0) e0[j] = 1;
        }
        if(_t0[j] != _t1[j])
            std::cout << "j = " << j << " -> " << e0[j] + e1[j] << " ----? " << (int) _t0[j] << " " << (int) _t1[j] << std::endl << std::endl;
    } 




    for(size_t j = 0; j < nitems; ++j)
    { 
        neg_e0[j] = -e0[j];
        neg_e1[j] = -e1[j];
        if(j == target_ind2) 
        {
            f1[j] = 1 - f0[j];
        }
        else
        {
          f1[j] = 0 - f0[j];  
        }

        e0_tilde[j] = e0[j] + Y0[j];
        e1_tilde[j] = e1[j] + Y1[j];

        e0_plus[j] = e0[j] + e1_tilde[j];
        e1_plus[j] = e1[j] + e0_tilde[j];
    }   



    uint64_t Y_i_star0, Y_i_star1;

    arc4random_buf(&Y_i_star0, sizeof(uint64_t));

    Y_i_star1 = Y0[target_ind2] + Y1[target_ind2] - Y_i_star0;
    uint64_t cancellation_correction0 = dotProductVector(f0, X0, nitems);
    uint64_t cancellation_correction1 = dotProductVector(f1, X1, nitems);

    arc4random_buf(&alpha, sizeof(alpha));
    auto gamma0 = dotProductVector(X0, Y1, nitems); // Held by P0
    auto gamma1 = dotProductVector(X1, Y0, nitems); // Held by P1
    
    std::cout << "gamma0 = " << gamma0 << std::endl;
    std::cout << "gamma1 = " << gamma1 << std::endl;

    auto cancellation_term0 = gamma0 + alpha;
    auto cancellation_term1 = gamma1 - alpha;


    uint64_t dot0 =  dotProductVector(DB0, e0_plus, nitems) - dotProductVector(Y0, DB1_tilde, nitems) + cancellation_term0;
    uint64_t dot1 =  dotProductVector(DB1, e1_plus, nitems) - dotProductVector(Y1, DB0_tilde, nitems) + cancellation_term1;

    uint64_t dot = dot0 + dot1;

    std::cout << "dot = " << dot << std::endl;
    std::cout << "cancellation_term0 = " << cancellation_term0 << std::endl;
    std::cout << "cancellation_term1 = " << cancellation_term1 << std::endl;

    refresh_blinds(e0,  e0_tilde, Y0,  f0, f1,  nitems);
    refresh_blinds(e1,  e1_tilde, Y1,  f1, f0,  nitems);



    refresh_blinds(e0,  e0_tilde, Y0,  neg_e0, neg_e1,  nitems);
    refresh_blinds(e1,  e1_tilde, Y1,  neg_e1, neg_e0,  nitems);


   uint64_t cancellation_correction0_new = cancellation_term0 - X0[target_ind2]  + dotProductVector(X0, f0, nitems) + X0[target_ind]  + dotProductVector(X0, neg_e0, nitems);
     uint64_t cancellation_correction1_new = cancellation_term1 - X1[target_ind2]  + dotProductVector(X1, f1, nitems)+ X1[target_ind]  + dotProductVector(X1, neg_e1, nitems);


    uint64_t assert1 = dotProductVector(X0, Y1,nitems) + dotProductVector(X1, Y0, nitems);
    uint64_t assert2 = cancellation_correction0_new + cancellation_correction1_new;
    
    std::cout << assert1 << " <<<>>>>> " << assert2 << std::endl;
    
    for(size_t j =  0; j < nitems; ++j)
    {
        uint64_t tmp = e0[j] + Y0[j];
        //std::cout << tmp << " <> " << e0_tilde[j] << std::endl;
        uint64_t reconstruction = e0[j] + e1[j];
        //std::cout << "reconstruction = " << reconstruction << std::endl;
        if(reconstruction != 0) std::cout << "j (reconstruction) = " << j << std::endl;
    }
    
    for(size_t j = 0; j < nitems; ++j)
    {
      e0_plus[j] = e0[j] + e1_tilde[j];
      e1_plus[j] = e1[j] + e0_tilde[j];
    }




  
     dot0 =  dotProductVector(DB0, e0_plus, nitems) - dotProductVector(Y0, DB1_tilde, nitems) + cancellation_correction0_new;
     dot1 =  dotProductVector(DB1, e1_plus, nitems) - dotProductVector(Y1, DB0_tilde, nitems) + cancellation_correction1_new;

     dot = dot0 + dot1;



    std::cout << "dot new = " << dot << std::endl;

    auto gamma0_new = dotProductVector(X0, Y1, nitems); // Held by P0
    auto gamma1_new = dotProductVector(X1, Y0, nitems); // Held by P1
    
    std::cout << "gamma0_mew = " << gamma0_new << std::endl;
    std::cout << "gamma1_new = " << gamma1_new << std::endl;



    proofdb = new leaf_t[nitems];
    std::cout << "nitems = " << nitems << std::endl;
    
    for (size_t j = 0; j < nitems; ++j) {
        uint64_t lo = 100 * j;
        uint64_t hi = 200 * j;
        proofdb[j] = _mm_set_epi64x(hi, lo); // high, low
    }
    
    auto start_P2 = std::chrono::high_resolution_clock::now();
  
    uint64_t alpha0, alpha1;
    arc4random_buf(&alpha0, sizeof(alpha0));
    alpha1 = alpha - alpha0;

    auto end_P2 = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double, std::milli> duration_ms_P2 = end_P2 - start_P2;
    std::cout << "P2 time: " << duration_ms_P2.count() << " ms\n";


    BB0 = new leaf_t[nitems]();
    BB1 = new leaf_t[nitems]();

  
    
    uint64_t target_ind0, target_ind1;
    uint64_t proof_submitted0, proof_submitted1;

    arc4random_buf(&proof_submitted0, sizeof(proof_submitted0));
    arc4random_buf(&proof_submitted1, sizeof(proof_submitted1));
    arc4random_buf(&target_ind0, sizeof(target_ind0));
    target_ind1 = target_ind - target_ind0;

    uint64_t offset = (alpha0 - target_ind0) - (alpha1 - target_ind1);

 
    
    // Timing dotProductVector
    auto start_P0_preproc = std::chrono::high_resolution_clock::now();
    
    __evalinterval(dpfkey0, 0, nitems - 1, output0, _t0);

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
    
    for(size_t j = 0; j < nitems; ++j)
    {
        if(output0[j][0] != output1[j][0] ||  output0[j][0] != output1[j][0])
        {
            std::cout << j << " :--> " << output0[j][0] << " || " << output0[j][1] << " <---> " << output1[j][0] << " || " << output1[j][1] << std::endl;
            std::cout << j << "(flags ): -> " << (int) _t0[j] << " <> " << (int) _t1[j] << std::endl;
        }
        else
        {
            assert(_t0[j] == _t1[j]);
        }
        
    }

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
