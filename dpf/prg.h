/* Copyright (C) 2019  Anonymous
 *
 * This is a pre-release version of the DPF++ library distributed anonymously
 * for peer review. A public release of the software will be published under the
 * LPGL v2.1 license in the near future. Please do not redistribute this version
 * of the software.
 */

#ifndef DPFPP_PRG_H__
#define DPFPP_PRG_H__

#include <cstring>
#include <cassert>
#include <cstdint>
#include <immintrin.h> // for __m128i intrinsics
//#include <openssl/aes.h> // if using OpenSSL AES functions
#include "aes.h"

namespace dpf
{

template<typename node_t, typename prgkey_t>
inline void PRG(const prgkey_t & prgkey, const node_t seed, void * outbuf, const uint32_t len, const uint32_t from = 0);

template<>
inline void PRG(const AES_KEY & prgkey, const __m128i seed, void * outbuf, const uint32_t len, const uint32_t from)
{
	__m128i * outbuf128 = reinterpret_cast<__m128i *>(outbuf);
	for (size_t i = 0; i < len; ++i)
	{
		outbuf128[i] = _mm_xor_si128(seed, _mm_set_epi64x(0, from+i));
	}
	AES_ecb_encrypt_blks(outbuf128, static_cast<unsigned int>(len), &prgkey);
	for (size_t i = 0; i < len; ++i)
	{
		outbuf128[i] = _mm_xor_si128(outbuf128[i], _mm_set_epi64x(0, from+i));
		outbuf128[i] = _mm_xor_si128(outbuf128[i], seed);
	}
} // PRG



inline void PRG_safe(const AES_KEY &prgkey,
                     const __m128i seed,
                     void *outbuf,
                     const uint32_t len_blocks,
                     const uint32_t from)
{
    assert(outbuf != nullptr);

    const size_t blocks = static_cast<size_t>(len_blocks);
    if (blocks == 0) return;

    // Work in bytes for the intrinsics, but cast to __m128i* for AES_ecb_encrypt_blks
    uint8_t *out = reinterpret_cast<uint8_t *>(outbuf);

    // Pre-encrypt XOR with seed and counter
    for (size_t i = 0; i < blocks; ++i) {
        __m128i tweak = _mm_set_epi64x(0, static_cast<long long>(from + i));
        __m128i pt = _mm_xor_si128(seed, tweak);
        // Unaligned store to avoid crashes if outbuf isn't 16-byte aligned
        _mm_storeu_si128(reinterpret_cast<__m128i *>(out + i * 16), pt);
    }

    // AES over the blocks; AES_ecb_encrypt_blks wants __m128i*
    AES_ecb_encrypt_blks(
        reinterpret_cast<__m128i *>(out),
        static_cast<unsigned int>(blocks),
        &prgkey
    );

    // Post-encrypt XORs
    for (size_t i = 0; i < blocks; ++i) {
        __m128i ct = _mm_loadu_si128(
            reinterpret_cast<const __m128i *>(out + i * 16)
        );
        __m128i tweak = _mm_set_epi64x(0, static_cast<long long>(from + i));
        ct = _mm_xor_si128(ct, tweak);
        ct = _mm_xor_si128(ct, seed);
        _mm_storeu_si128(reinterpret_cast<__m128i *>(out + i * 16), ct);
    }
}



} // namespace dpf

#endif




// #include <cassert>
// #include <cstdint>
// #include <immintrin.h> // for __m128i intrinsics
// #include <openssl/aes.h> // if using OpenSSL AES functions

