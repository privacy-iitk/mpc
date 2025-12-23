 
#ifndef DPFPP_DPF_H__
#define DPFPP_DPF_H__
#pragma once
#include <type_traits>  // std::is_same<>
#include <limits>       // std::numeric_limits<>
#include <climits>      // CHAR_BIT
#include <cmath>        // std::log2, std::ceil, std::floor
#include <stdexcept>    // std::runtime_error
#include <array>        // std::array<>
#include <iostream>     // std::istream and std::ostream
#include <vector>       // std::vector<>
#include <memory>       // std::shared_ptr<>
#include <utility>      // std::move
#include <algorithm>    // std::copy

#include <bsd/stdlib.h> // arc4random_buf
#include <x86intrin.h>  // SSE and AVX intrinsics

#include "aes.h"
#include "block.h"
#include "prg.h"

#define L 0
#define R 1
 
namespace dpf
{

template<typename leaf_t = __m128i, typename node_t = __m128i, typename prgkey_t = AES_KEY>
struct dpf_key;

template<typename leaf_t, typename node_t, typename prgkey_t>
inline leaf_t eval(const dpf_key <leaf_t, node_t, prgkey_t> & dpfkey, const size_t input);

template<typename leaf_t, typename node_t, typename prgkey_t>
inline void __evalinterval(const dpf_key<leaf_t, node_t, prgkey_t> & dpfkey, const size_t from, const size_t to, leaf_t * output, uint8_t * _t);

template<typename node_t, typename prgkey_t>
static inline void expand(const prgkey_t & prgkey, const node_t & seed, node_t s[2], uint8_t t[2])
{
	dpf::PRG(prgkey, clear_lsb(seed, 0b11), s, 2);
	t[L] = get_lsb(s[L]);
	s[L] = clear_lsb(s[L], 0b11);
	t[R] = get_lsb(s[R]);
	s[R] = clear_lsb(s[R], 0b11);
} // dpf::expand

 
template< typename prgkey_t>
static inline void traverse_(const prgkey_t & prgkey, __m128i & seed, const bool direction,
	const uint8_t cw_t, const __m128i & cw, const uint8_t prev_t, __m128i & s, uint8_t & t)
{
	dpf::PRG(prgkey, clear_lsb(seed, 0b11), &s, 1, direction);
	t = get_lsb(s) ^ (cw_t & prev_t);
	s = clear_lsb(xor_if(s, cw, prev_t), 0b11);
} // dpf::traverse

template<typename node_t, typename prgkey_t>
static inline void traverse2(const prgkey_t & prgkey, const __m128i & seed,
	const uint8_t cw_t[2], const node_t & cw, const uint8_t prev_t,
	__m128i s[2], uint8_t t[2], int lsbmask = 0b11)
{
	dpf::PRG(prgkey, clear_lsb(seed, 0b11), s, 2);
	t[L] = get_lsb(s[L]) ^ (cw_t[L] & prev_t);;
	s[L] = clear_lsb(xor_if(s[L], cw, prev_t), lsbmask);
	t[R] = get_lsb(s[R]) ^ (cw_t[R] & prev_t);;
	s[R] = clear_lsb(xor_if(s[R], cw, prev_t), lsbmask);
} // dpf::expand

template<typename node_t, typename prgkey_t, size_t nodes_per_leaf>
static inline void stretch_leaf(const prgkey_t & prgkey, const node_t & seed, std::array<node_t, nodes_per_leaf> & s)
{
	dpf::PRG(prgkey, clear_lsb(seed), &s, nodes_per_leaf);
} // dpf::stretch_leaf



template<typename leaf_t, typename node_t, typename prgkey_t>
struct dpf_key final
{
  public:
	static constexpr size_t bits_per_leaf = std::is_same<leaf_t, bool>::value ? 1 : sizeof(leaf_t) * CHAR_BIT;
	static constexpr bool is_packed = (sizeof(leaf_t) < sizeof(node_t));
	static constexpr size_t leaves_per_node = dpf_key::is_packed ? sizeof(node_t) * CHAR_BIT / bits_per_leaf : 1;
	static constexpr size_t nodes_per_leaf = dpf_key::is_packed ? 1 : std::ceil(static_cast<double>(bits_per_leaf) / (sizeof(node_t) * CHAR_BIT));
	inline static constexpr size_t nodes_in_interval(const size_t from, const size_t to) { return (to < from) ? 0 : std::max(1.0, std::ceil(static_cast<double>(to+1) / leaves_per_node) - std::floor(static_cast<double>(from) / leaves_per_node)); }
	static_assert(leaves_per_node * bits_per_leaf == sizeof(node_t) * CHAR_BIT
	    || nodes_per_leaf * sizeof(node_t) == sizeof(leaf_t));

	using block_t = node_t;
	using finalizer_t = std::array<block_t, nodes_per_leaf>;
	typedef std::pair<finalizer_t, finalizer_t> (*finalizer_callback)(const prgkey_t &, const size_t, const leaf_t &, const block_t[2], const uint8_t[2]);

	inline static constexpr size_t depth(const size_t nitems) { return std::ceil(std::log2(std::ceil(static_cast<double>(nitems) / dpf_key::leaves_per_node))); }
	inline constexpr size_t depth() const { return dpf_key::depth(nitems); }

	inline static constexpr size_t input_bits(const size_t nitems) { return std::ceil(std::log2(nitems)); }
	inline constexpr size_t input_bits() const { return dpf_key::input_bits(nitems); }

	inline dpf_key(dpf_key &&) = default;
	inline dpf_key & operator=(dpf_key &&) = default;
	inline dpf_key(const dpf_key &) = default;
	inline dpf_key & operator=(const dpf_key &) = default;

	inline bool operator==(const dpf_key & rhs) const { return nitems == rhs.nitems && root == rhs.root && cw == rhs.cw && finalizer == rhs.finalizer; }
	inline bool operator!=(const dpf_key & rhs) const { return !(*this == rhs); }

 
 
	static auto default_make_finalizer(const prgkey_t & prgkey, const size_t target,
									const leaf_t & val, const block_t s[2], const uint8_t t[2])
	{
		std::cout << "default_make_finalizer" << std::endl;
		finalizer_t finalizer;

		finalizer_t stretched[2];
		stretch_leaf(prgkey, s[L], stretched[L]);
		stretch_leaf(prgkey, s[R], stretched[R]);

		if constexpr (dpf_key::is_packed) {
			std::cout << "it is indeed packed!" << std::endl;
			auto finalizer0 = reinterpret_cast<node_t *>(&finalizer[0]);
			if constexpr (std::numeric_limits<leaf_t>::is_integer) {
				if constexpr (std::is_same<leaf_t, bool>::value) {
					*finalizer0 = val ? 1 : 0;
				} else {
					using unsigned_leaf_t = std::make_unsigned_t<leaf_t>;
					*finalizer0 = static_cast<unsigned_leaf_t>(val);
				}
				finalizer0->shiftl(bits_per_leaf * (target % leaves_per_node));
			} else {
				*finalizer0 = val;
			}
		} else {
			// generic memory copy
			static_assert(sizeof(finalizer_t) == sizeof(val),
						"finalizer_t and leaf_t must have same size");
			std::memcpy(&finalizer, &val, sizeof(finalizer_t));
		}

		for (size_t j = 0; j < nodes_per_leaf; ++j) {
			finalizer[j] ^= stretched[L][j] ^ stretched[R][j];
		}
		return std::make_pair(finalizer, finalizer);
	}





	static auto gen(const prgkey_t & prgkey, size_t nitems, size_t target, const leaf_t & val = 1, const finalizer_callback make_finalizer = default_make_finalizer)
	{
		if (nitems <= target)
		{
			throw std::runtime_error("target point out of range");
		}

		block_t root[2];
		arc4random_buf(root, sizeof(root));
		uint8_t t[2] = { get_lsb(root[0]), !t[0] };
		root[1] = set_lsb(root[1], t[1]);
		block_t s[2] = { root[0], root[1] };
		
		#ifdef VERBOSE
		std::cout << "root0 = " << std::hex << root[0][0] << (root[0][1] & 0xffffffffULL) << std::endl;
		std::cout << "root1 = " << std::hex << root[1][0] << (root[1][1] & 0xffffffffULL) << std::endl;
		std::cout << "t0 = " << (int) t[0] << std::endl;
		std::cout << "t1 = " << (int) t[1] << std::endl;
		#endif

		const size_t depth = dpf_key::depth(nitems);
		std::vector<block_t> cw;
		cw.reserve(depth);

		block_t s0[2], s1[2];
		uint8_t t0[2], t1[2];
		const size_t nbits = input_bits(nitems);
		for (size_t layer = 0; layer < depth; ++layer)
		{
			#ifdef VERBOSE
			std::cout << std::endl << " ------------ " << std::endl;
			#endif

			const uint8_t bit = (target >> (nbits - layer - 1)) & 1U;
		 	
			#ifdef VERBOSE
			std::cout << "bit = " << (int) bit << std::endl;
			#endif

			expand(prgkey, s[0], s0, t0);
			expand(prgkey, s[1], s1, t1);

			#ifdef VERBOSE
			std::cout << "P0: left child: "  << (s0[0][0] & 0xffffffffULL) << (s0[0][1] & 0xffffffffULL) << std::endl;
			std::cout << "P0: right child: " << (s0[1][0] & 0xffffffffULL) << (s0[1][1] & 0xffffffffULL) << std::endl << std::endl;

			std::cout << "P1: left child: "  << (s1[0][0] & 0xffffffffULL) << (s1[0][1] & 0xffffffffULL) << std::endl;
			std::cout << "P1: right child: " << (s1[1][0] & 0xffffffffULL) << (s1[1][1] & 0xffffffffULL)<< std::endl;
			
			std::cout << "P0 L advice bit = " << (int) t0[L] << std::endl; 
			std::cout << "P0 R advice bit = " << (int) t0[R] << std::endl;

			std::cout << "P1 L advice bit = " << (int) t1[L] << std::endl; 
			std::cout << "P1 R advice bit = " << (int) t1[R] << std::endl;
			#endif

			const uint8_t keep = (bit == 0) ? L : R, lose = 1 - keep;
			bool cwt[2] = {
			    cwt[L] = t0[L] ^ t1[L] ^ bit ^ 1,
			    cwt[R] = t0[R] ^ t1[R] ^ bit
			};
			auto nextcw = s0[lose] ^ s1[lose];

			#ifdef VERBOSE
			std::cout << "cw = " << (nextcw[0]& 0xffffffffULL) << (nextcw[1]& 0xffffffffULL) << std::endl;
			std::cout << cwt[L] << ", " << cwt[R] << std::endl;
			#endif

			s[L] = xor_if(s0[keep], nextcw, t[L]);
			t[L] = t0[keep] ^ (t[L] & cwt[keep]);

			s[R] = xor_if(s1[keep], nextcw, t[R]);
			t[R] = t1[keep] ^ (t[R] & cwt[keep]);

			#ifdef VERBOSE
			std::cout << "after the correction words are applied: " << std::endl;
			std::cout << "SL = " << (s[L][0] & 0xffffffffULL) << (s[L][1] & 0xffffffffULL) << std::endl; 
			std::cout << "SR = " << (s[R][0] & 0xffffffffULL) << (s[R][1] & 0xffffffffULL) << std::endl;
			std::cout << "tL = " << (int) t[L] << std::endl; 
			std::cout << "tR = " << (int) t[R] << std::endl;
			#endif
			
			cw.emplace_back(set_lsbs(nextcw, cwt));
		}
		auto [finalizer0, finalizer1] = make_finalizer(prgkey, target, val, s, t);

		return std::make_pair(
		    std::forward<dpf_key>(dpf_key(nitems, root[0], cw, finalizer0, prgkey)),
		    std::forward<dpf_key>(dpf_key(nitems, root[1], cw, finalizer1, prgkey)));
	} // dpf_key::gen

	inline leaf_t eval(const size_t input) const { return std::forward<leaf_t>(dpf::eval(*this, input)); }
	
	const size_t nitems;
	const block_t root;
	const std::vector<block_t> cw;
	const finalizer_t finalizer;
	const prgkey_t prgkey;

  public:
	dpf_key(size_t nitems_, const block_t & root_, const std::vector<block_t> cw_,
		const finalizer_t & finalizer_, const prgkey_t & prgkey_)
	  : nitems(nitems_),
	    root(root_),
	    cw(cw_),
	    finalizer(finalizer_),
	    prgkey(prgkey_) { }
		
}; // struct dpf::dpf_key

 

	template<typename leaf_t, typename node_t>
	inline leaf_t getword(const node_t & S, const size_t input)
	{
		auto S_ = reinterpret_cast<const leaf_t *>(&S);
		if constexpr(sizeof(leaf_t) >= sizeof(node_t)) return *S_;

		return S_[input];
	} // dpf::getword

	template<>
	inline bool getword(const block<__m128i> & S, const size_t input)
	{
		const __m128i mask = bool128_mask[input / 64];
		__m128i vcmp = _mm_xor_si128(_mm_and_si128(S >> (input % 64), mask), mask);

		return static_cast<bool>(_mm_testz_si128(vcmp, vcmp));
	} // dpf::getword<__m128i,bool>

	template<>
	inline bool getword(const block<__m256i> & S, const size_t input)
	{
		const __m256i mask = bool256_mask[input / 64];
		__m256i vcmp = _mm256_xor_si256(_mm256_and_si256(S >> (input % 64), mask), mask);

		return static_cast<bool>(_mm256_testz_si256(vcmp, vcmp));
	} // dpf::getword<__m256i,bool>

	template<typename leaf_t, typename node_t, typename prgkey_t>
	inline void finalize(const prgkey_t & prgkey, std::array<node_t, dpf_key<leaf_t, node_t, prgkey_t>::nodes_per_leaf> finalizer, leaf_t * output, node_t * s, size_t nnodes, uint8_t * t)
	{
		auto output_ = reinterpret_cast<std::array<node_t, dpf_key<leaf_t, node_t, prgkey_t>::nodes_per_leaf> *>(output);
		for (size_t i = 0; i < nnodes; ++i)
		{
			stretch_leaf(prgkey, s[i], output_[i]);
			for (size_t j = 0; j < dpf_key<leaf_t, node_t, prgkey_t>::nodes_per_leaf; ++j)
			{
				output_[i][j] = xor_if(output_[i][j], finalizer[j], t[i]);
			}
		}
	} // dpf::finalize

 

	template<typename leaf_t, typename node_t, typename prgkey_t>
	inline leaf_t eval(const dpf_key<leaf_t, node_t, prgkey_t> & dpfkey, const size_t input)
	{
		auto prgkey = dpfkey.prgkey;
		auto root = dpfkey.root;
		auto depth = dpfkey.depth();
		auto nbits = dpfkey.input_bits();

		node_t S = root;
		uint8_t T = get_lsb(root, 0b01);

		for (size_t layer = 0; layer < depth; ++layer)
		{
			auto & cw = dpfkey.cw[layer];
			const uint8_t nextbit = (input >> (nbits-layer-1)) & 1;
			traverse_(prgkey, S, nextbit, get_lsb(cw, nextbit ? 0b10 : 0b01), cw, T, S, T); 
		}

		std::array<node_t, dpf_key<leaf_t, node_t, prgkey_t>::nodes_per_leaf> final;
		finalize(prgkey, dpfkey.finalizer, &final, &S, 1, &T);

		if constexpr(dpfkey.is_packed)
		{
			auto S_ = reinterpret_cast<node_t *>(&final);
			return std::forward<leaf_t>(getword<leaf_t>(*S_, input % dpfkey.leaves_per_node));
		}
		else
		{
			auto ret = reinterpret_cast<leaf_t *>(&final);
			return *ret;
		}
	} // dpf::eval


	template<typename leaf_t, typename node_t, typename prgkey_t>
	inline void __evalinterval(const dpf_key<leaf_t, node_t, prgkey_t> & dpfkey, const size_t from, const size_t to, leaf_t * output, uint8_t * _t)
	{
		auto nodes_per_leaf = dpfkey.nodes_per_leaf;
		std::cout << "nodes_per_leaf = " << nodes_per_leaf << std::endl;
		auto depth = dpfkey.depth();
		auto nbits = dpfkey.input_bits();
		auto nodes_in_interval = dpfkey.nodes_in_interval(from, to);
		auto root = dpfkey.root;
		auto prgkey = dpfkey.prgkey;

		const size_t from_node = std::floor(static_cast<double>(from) / nodes_per_leaf);

		node_t * s[2] = {
			reinterpret_cast<node_t *>(output) + nodes_in_interval * (nodes_per_leaf - 1),
			s[0] + nodes_in_interval / 2
		};
		uint8_t * t[2] = { _t, _t + nodes_in_interval / 2};

		int curlayer = depth % 2;

		s[curlayer][0] = root;
		t[curlayer][0] = get_lsb(root, 0b01);

	

		for (size_t layer = 0; layer < depth; ++layer)
		{
			auto & cw = dpfkey.cw[layer];
			uint8_t cw_t[2] = { get_lsb(cw, 0b01), get_lsb(cw, 0b10) };
			curlayer = 1-curlayer;

			size_t i=0, j=0;
			auto nextbit = (from_node >> (nbits-layer-1)) & 1;
			size_t nodes_in_prev_layer = std::ceil(static_cast<double>(nodes_in_interval) / (1ULL << (depth-layer)));
			size_t nodes_in_cur_layer = std::ceil(static_cast<double>(nodes_in_interval) / (1ULL << (depth-layer-1)));

		// 	//printf("nextbit = %u\n",  (from_node >> (nbits-layer-1)));

			if (nextbit == 1) traverse_(prgkey, s[1-curlayer][0], R, cw_t[R], cw, t[1-curlayer][j], s[curlayer][0], t[curlayer][0]); // these will not be called in evalfull
			for (i = nextbit, j = nextbit; j < nodes_in_prev_layer-1; ++j, i+=2)
			{
				//printf("j = %u\n", j );
				traverse2(prgkey, s[1-curlayer][j], cw_t, cw, t[1-curlayer][j], &s[curlayer][i], &t[curlayer][i]);
			}
			if (nodes_in_prev_layer > j)
			{
				//printf("jj' = %u\n", j );
				if (i < nodes_in_cur_layer - 1) 
				{
					traverse2(prgkey, s[1-curlayer][j], cw_t, cw, t[1-curlayer][j], &s[curlayer][i], &t[curlayer][i]);
					//printf("If\n");
				}
				else
				{
					traverse_(prgkey, s[1-curlayer][j], L, cw_t[L], cw, t[1-curlayer][j], s[curlayer][i], t[curlayer][i]); // will not be called in evalfull
					printf("else\n");
				} 
			}
		}   
		finalize(prgkey, dpfkey.finalizer, output, s[0], nodes_in_interval, t[0]);
	} // dpf::__evalinterval
		// ---------------------------------------------------------
	// Serialization
 // --- Helpers to serialize/deserialize __m128i safely ---
// --- Serialize and deserialize __m128i safely ---
inline void serialize_m128i(std::vector<uint8_t>& buf, const __m128i& val) {
    alignas(16) uint8_t tmp[16];
    _mm_storeu_si128(reinterpret_cast<__m128i*>(tmp), val);
    buf.insert(buf.end(), tmp, tmp + 16);
}

inline __m128i deserialize_m128i(const uint8_t*& ptr) {
    __m128i val = _mm_loadu_si128(reinterpret_cast<const __m128i*>(ptr));
    ptr += 16;
    return val;
}

// --- Generic byte serialization ---
template<typename T>
inline void serialize_bytes(std::vector<uint8_t>& buf, const T& val) {
    static_assert(std::is_trivially_copyable_v<T>, "serialize_bytes requires trivially copyable type");
    const uint8_t* p = reinterpret_cast<const uint8_t*>(&val);
    buf.insert(buf.end(), p, p + sizeof(T));
}

template<typename T>
inline T deserialize_bytes(const uint8_t*& ptr) {
    static_assert(std::is_trivially_copyable_v<T>, "deserialize_bytes requires trivially copyable type");
    T val;
    std::memcpy(&val, ptr, sizeof(T));
    ptr += sizeof(T);
    return val;
}

// --- Serialize DPF key ---
template<typename leaf_t, typename node_t, typename prgkey_t>
std::vector<uint8_t> serialize_dpf_key(const dpf_key<leaf_t, node_t, prgkey_t>& key) {
    std::vector<uint8_t> buf;

    // nitems
    serialize_bytes(buf, key.nitems);

    // root
    if constexpr (std::is_same_v<node_t, __m128i>)
        serialize_m128i(buf, key.root);
    else
        serialize_bytes(buf, key.root);

    // correction words
    uint64_t cw_size = key.cw.size();
    serialize_bytes(buf, cw_size);
    for (const auto& cw : key.cw) {
        if constexpr (std::is_same_v<node_t, __m128i>)
            serialize_m128i(buf, cw);
        else
            serialize_bytes(buf, cw);
    }

    // finalizer (array<node_t, N>)
    for (const auto& f : key.finalizer) {
        if constexpr (std::is_same_v<node_t, __m128i>)
            serialize_m128i(buf, f);
        else
            serialize_bytes(buf, f);
    }

    // prgkey
    serialize_bytes(buf, key.prgkey);
    return buf;
}

// --- Deserialize DPF key ---
template<typename leaf_t, typename node_t, typename prgkey_t>
dpf_key<leaf_t, node_t, prgkey_t> deserialize_dpf_key(const uint8_t* data, size_t len) {
    const uint8_t* ptr = data;

    // nitems
    uint64_t nitems = deserialize_bytes<uint64_t>(ptr);

    // root
    node_t root;
    if constexpr (std::is_same_v<node_t, __m128i>)
        root = deserialize_m128i(ptr);
    else
        root = deserialize_bytes<node_t>(ptr);

    // correction words
    uint64_t cw_size = deserialize_bytes<uint64_t>(ptr);
    std::vector<node_t> cw(cw_size);
    for (size_t i = 0; i < cw_size; ++i) {
        if constexpr (std::is_same_v<node_t, __m128i>)
            cw[i] = deserialize_m128i(ptr);
        else
            cw[i] = deserialize_bytes<node_t>(ptr);
    }

    // finalizer
    typename dpf_key<leaf_t, node_t, prgkey_t>::finalizer_t finalizer;
    for (auto& f : finalizer) {
        if constexpr (std::is_same_v<node_t, __m128i>)
            f = deserialize_m128i(ptr);
        else
            f = deserialize_bytes<node_t>(ptr);
    }

    // prgkey
    prgkey_t prgkey = deserialize_bytes<prgkey_t>(ptr);

    return dpf_key<leaf_t, node_t, prgkey_t>(nitems, root, cw, finalizer, prgkey);
}


} // namespace dpf

#endif
