#pragma once

#include <immintrin.h>
#include <cstdint>
#include <vector>
#include <type_traits>
#include <stdexcept>

#include <array>
#include <type_traits>
#include <cstring>
#include <sstream>
#include <iomanip>

// -----------------------------------------------------------------------------
// Trait: detect std::array<__m128i, N>
// -----------------------------------------------------------------------------
template<typename T> struct is_array_of_m128i : std::false_type { static constexpr size_t N = 0; };
template<size_t N> struct is_array_of_m128i<std::array<__m128i, N>> : std::true_type { static constexpr size_t N_v = N; };

// -----------------------------------------------------------------------------
// leaf_equal - compare two leaves for equality
// Supports: integral types (uint8_t/uint32_t/uint64_t...), __m128i, and std::array<__m128i,N>
// -----------------------------------------------------------------------------
template<typename Leaf>
bool leaf_equal(const Leaf &a, const Leaf &b)
{
    if constexpr (std::is_integral_v<Leaf>) {
        return a == b;
    } else if constexpr (std::is_same_v<Leaf, __m128i>) {
        // safe: compare raw bytes
        return std::memcmp(&a, &b, sizeof(__m128i)) == 0;
    } else if constexpr (is_array_of_m128i<Leaf>::value) {
        constexpr size_t N = is_array_of_m128i<Leaf>::N_v;
        for (size_t i = 0; i < N; ++i) {
            if (std::memcmp(&a[i], &b[i], sizeof(__m128i)) != 0) return false;
        }
        return true;
    } else {
        static_assert(!std::is_same_v<Leaf, Leaf>, "leaf_equal: unsupported Leaf type");
    }
}

// -----------------------------------------------------------------------------
// leaf_xor - elementwise XOR of two leaves, returns same Leaf type
// -----------------------------------------------------------------------------
template<typename Leaf>
Leaf leaf_xor(const Leaf &a, const Leaf &b)
{
    if constexpr (std::is_integral_v<Leaf>) {
        return static_cast<Leaf>(a ^ b);
    } else if constexpr (std::is_same_v<Leaf, __m128i>) {
        return _mm_xor_si128(a, b);
    } else if constexpr (is_array_of_m128i<Leaf>::value) {
        constexpr size_t N = is_array_of_m128i<Leaf>::N_v;
        Leaf out;
        for (size_t i = 0; i < N; ++i) out[i] = _mm_xor_si128(a[i], b[i]);
        return out;
    } else {
        static_assert(!std::is_same_v<Leaf, Leaf>, "leaf_xor: unsupported Leaf type");
    }
}

// -----------------------------------------------------------------------------
// leaf_to_hex_string - small helper to format a leaf as hex for debugging.
// For integral types it prints the integer. For __m128i/array it prints bytes
// in hex. This is only for debugging (verbose printing).
// -----------------------------------------------------------------------------
template<typename Leaf>
std::string leaf_to_hex_string(const Leaf &x)
{
    std::ostringstream ss;
    ss << std::hex << std::setfill('0');

    if constexpr (std::is_integral_v<Leaf>) {
        ss << "0x" << +x; // + to avoid char printing for small ints
        return ss.str();
    } else if constexpr (std::is_same_v<Leaf, __m128i>) {
        alignas(16) uint8_t bytes[16];
        _mm_storeu_si128((__m128i*)bytes, x);
        for (int i = 0; i < 16; ++i) ss << std::setw(2) << (int)bytes[i];
        return ss.str();
    } else if constexpr (is_array_of_m128i<Leaf>::value) {
        constexpr size_t N = is_array_of_m128i<Leaf>::N_v;
        for (size_t k = 0; k < N; ++k) {
            alignas(16) uint8_t bytes[16];
            _mm_storeu_si128((__m128i*)bytes, x[k]);
            ss << "{" ;
            for (int i = 0; i < 16; ++i) ss << std::setw(2) << (int)bytes[i];
            ss << "}";
            if (k+1 < N) ss << ",";
        }
        return ss.str();
    } else {
        static_assert(!std::is_same_v<Leaf, Leaf>, "leaf_to_hex_string: unsupported Leaf type");
    }
}



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


static const __m128i bool128_mask[2] = {
	_mm_set_epi64x(0,1),                                        // 0b00...0001
	_mm_set_epi64x(1,0)                                         // 0b00...0001 << 64
};

static const __m128i lsb128_mask[4] = {
	_mm_setzero_si128(),                                        // 0b00...0000
	_mm_set_epi64x(0,1),                                        // 0b00...0001
	_mm_set_epi64x(0,2),                                        // 0b00...0010
	_mm_set_epi64x(0,3)                                         // 0b00...0011
};
static const __m128i lsb128_mask_inv[4] = {
	_mm_set1_epi8(-1),                                          // 0b11...1111
	_mm_set_epi64x(-1,-2),                                      // 0b11...1110
	_mm_set_epi64x(-1,-3),                                      // 0b11...1101
	_mm_set_epi64x(-1,-4)                                       // 0b11...1100
};
static const __m128i if128_mask[2] = {
	_mm_setzero_si128(),                                        // 0b00...0000
	_mm_set1_epi8(-1)                                           // 0b11...1111
};

static const __m256i bool256_mask[4] = {
	_mm256_set_epi64x(0,0,0,1),                                 // 0b00...0001
	_mm256_set_epi64x(0,0,1,0),                                 // 0b00...0001 << 64
	_mm256_set_epi64x(0,1,0,0),                                 // 0b00...0001 << 128
	_mm256_set_epi64x(1,0,0,0)                                  // 0b00...0001 << 192
};

static const __m256i lsb256_mask[4] = {
	_mm256_setzero_si256(),                                     // 0b00...0000
	_mm256_set_epi64x(0,0,0,1),                                 // 0b00...0001
	_mm256_set_epi64x(0,0,0,2),                                 // 0b00...0010
	_mm256_set_epi64x(0,0,0,3)                                  // 0b00...0011
};
static const __m256i lsb256_mask_inv[4] = {
	_mm256_set1_epi8(-1),                                       // 0b11...1111
	_mm256_set_epi64x(-1,-1,-1,-2),                             // 0b11...1110
	_mm256_set_epi64x(-1,-1,-1,-3),                             // 0b11...1101
	_mm256_set_epi64x(-1,-1,-1,-4)                              // 0b11...1100
};
static const __m256i if256_mask[2] = {
	_mm256_setzero_si256(),                                     // 0b00...0000
	_mm256_set1_epi8(-1)                                        // 0b11...1111
};

 

inline __m128i xor_if(const __m128i & block1, const __m128i & block2, bool flag)
{
	return  _mm_xor_si128(block1, _mm_and_si128(block2, if128_mask[flag ? 1 : 0]));
}
inline __m256i xor_if(const __m256i & block1, const __m256i & block2, bool flag)
{
	return _mm256_xor_si256(block1, _mm256_and_si256(block2, if256_mask[flag ? 1 : 0]));
}

inline uint8_t get_lsb(const __m128i & block, uint8_t bits = 0b01)
{
	__m128i vcmp = _mm_xor_si128(_mm_and_si128(block, lsb128_mask[bits]), lsb128_mask[bits]);
	return static_cast<uint8_t>(_mm_testz_si128(vcmp, vcmp));
}
inline uint8_t get_lsb(const __m256i & block, uint8_t bits = 0b01)
{
	__m256i vcmp = _mm256_xor_si256(_mm256_and_si256(block, lsb256_mask[bits]), lsb256_mask[bits]);
	return static_cast<uint8_t>(_mm256_testz_si256(vcmp, vcmp));
}

inline __m128i clear_lsb(const __m128i & block, uint8_t bits = 0b01)
{
	return _mm_and_si128(block, lsb128_mask_inv[bits]);
}
inline __m256i clear_lsb(const __m256i & block, uint8_t bits = 0b01)
{
	return _mm256_and_si256(block, lsb256_mask_inv[bits]);
}

__m128i set_lsb(const __m128i & block, const bool val = true);
inline __m128i set_lsb(const __m128i & block, const bool val)
{
	return _mm_or_si128(clear_lsb(block, 0b01), lsb128_mask[val ? 0b01 : 0b00]);
}
__m256i set_lsb(const __m256i & block, const bool val = true);
inline __m256i set_lsb(const __m256i & block, const bool val)
{
	return _mm256_or_si256(clear_lsb(block, 0b01), lsb256_mask[val ? 0b01 : 0b00]);;
}

inline __m128i set_lsbs(const __m128i & block, const bool bits[2])
{
	int i = (bits[0] ? 1 : 0) + 2 * (bits[1] ? 1 : 0);
	return _mm_or_si128(clear_lsb(block, 0b11), lsb128_mask[i]);
}
inline __m256i set_lsbs(const __m256i & block, const bool bits[2])
{
	int i = (bits[0] ? 1 : 0) + 2 * (bits[1] ? 1 : 0);
	return _mm256_or_si256(clear_lsb(block, 0b11), lsb256_mask[i]);
}


// ===============================================================
// 128-bit modular arithmetic (mod 2^128) using __m128i
// ===============================================================
//
// These functions perform modular arithmetic on 128-bit integers,
// represented as __m128i, interpreted as unsigned 128-bit values.
//
// The arithmetic is done mod 2^128 — i.e., overflow wraps around.
// No carry/borrow propagation beyond 128 bits is performed.
//

// 128-bit multiplication mod 2^128
static inline __m128i mul128_mod2_128(__m128i a, __m128i b) {
    // Extract 64-bit halves
    uint64_t a_lo = _mm_cvtsi128_si64(a);
    uint64_t a_hi = _mm_extract_epi64(a, 1);
    uint64_t b_lo = _mm_cvtsi128_si64(b);
    uint64_t b_hi = _mm_extract_epi64(b, 1);

    __extension__ typedef unsigned __int128 u128;

    // Compute cross-products, ignoring high*high term (mod 2^128)
    u128 lo   = (u128)a_lo * b_lo;
    u128 mid1 = (u128)a_lo * b_hi;
    u128 mid2 = (u128)a_hi * b_lo;
    u128 r = lo + ((mid1 + mid2) << 64);

    // Return as __m128i
    uint64_t r_lo = (uint64_t)r;
    uint64_t r_hi = (uint64_t)(r >> 64);
    return _mm_set_epi64x(r_hi, r_lo);
}

// 128-bit addition mod 2^128
static inline __m128i add128_mod2_128(__m128i a, __m128i b) {
    uint64_t a_lo = _mm_cvtsi128_si64(a);
    uint64_t a_hi = _mm_extract_epi64(a, 1);
    uint64_t b_lo = _mm_cvtsi128_si64(b);
    uint64_t b_hi = _mm_extract_epi64(b, 1);

    uint64_t lo = a_lo + b_lo;
    uint64_t carry = (lo < a_lo);
    uint64_t hi = a_hi + b_hi + carry;

    return _mm_set_epi64x(hi, lo);
}

// 128-bit subtraction mod 2^128
static inline __m128i sub128_mod2_128(__m128i a, __m128i b) {
    uint64_t a_lo = _mm_cvtsi128_si64(a);
    uint64_t a_hi = _mm_extract_epi64(a, 1);
    uint64_t b_lo = _mm_cvtsi128_si64(b);
    uint64_t b_hi = _mm_extract_epi64(b, 1);

    uint64_t lo = a_lo - b_lo;
    uint64_t borrow = (a_lo < b_lo);
    uint64_t hi = a_hi - b_hi - borrow;

    return _mm_set_epi64x(hi, lo);
}

// ===============================================================
// mX: 128-bit modular arithmetic element
// ===============================================================
//
// mX wraps __m128i to provide operator overloading for arithmetic
// and bitwise operations mod 2^128. This makes __m128i behave like
// a high-level numeric type.
//
// Typical use:
//     mX a(hi, lo), b(hi, lo);
//     mX c = a + b;
//     mX d = a * b;
//     mX e = a ^ b;  // bitwise XOR
//

struct mX {
    __m128i v;

    // Constructors
    mX() = default;
    explicit mX(__m128i x) : v(x) {}
    mX(uint64_t hi, uint64_t lo) { v = _mm_set_epi64x(hi, lo); }

    // Accessors (split into hi/lo 64-bit halves)
    uint64_t lo() const { return _mm_cvtsi128_si64(v); }
    uint64_t hi() const { return _mm_extract_epi64(v, 1); }

    // ------------------------------------------------------------------
    // Arithmetic operators (+, -, *) mod 2^128
    // ------------------------------------------------------------------
    inline mX operator+(const mX& other) const { return mX(add128_mod2_128(v, other.v)); }
    inline mX operator-(const mX& other) const { return mX(sub128_mod2_128(v, other.v)); }
    inline mX operator*(const mX& other) const { return mX(mul128_mod2_128(v, other.v)); }

    // Compound assignment versions
    inline mX& operator+=(const mX& other) { v = add128_mod2_128(v, other.v); return *this; }
    inline mX& operator-=(const mX& other) { v = sub128_mod2_128(v, other.v); return *this; }
    inline mX& operator*=(const mX& other) { v = mul128_mod2_128(v, other.v); return *this; }

    // ------------------------------------------------------------------
    // Scalar (uint8_t) operations (for convenience)
    // ------------------------------------------------------------------
    inline mX operator*(uint8_t s) const {
        __m128i b = _mm_set_epi64x(0, (uint64_t)s);
        return mX(mul128_mod2_128(v, b));
    }
    inline mX operator+(uint8_t s) const {
        __m128i b = _mm_set_epi64x(0, (uint64_t)s);
        return mX(add128_mod2_128(v, b));
    }
    inline mX operator-(uint8_t s) const {
        __m128i b = _mm_set_epi64x(0, (uint64_t)s);
        return mX(sub128_mod2_128(v, b));
    }

    // Symmetric scalar operations (uint8_t op mX)
    friend inline mX operator*(uint8_t s, const mX& x) { return x * s; }
    friend inline mX operator+(uint8_t s, const mX& x) { return x + s; }
    friend inline mX operator-(uint8_t s, const mX& x) {
        __m128i a = _mm_set_epi64x(0, (uint64_t)s);
        return mX(sub128_mod2_128(a, x.v));
    }

    // ------------------------------------------------------------------
    // Bitwise operators (&, |, ^, ~)
    // ------------------------------------------------------------------
    inline mX operator&(const mX& other) const { return mX(_mm_and_si128(v, other.v)); }
    inline mX operator|(const mX& other) const { return mX(_mm_or_si128(v, other.v)); }
    inline mX operator^(const mX& other) const { return mX(_mm_xor_si128(v, other.v)); }
    inline mX operator~() const { return mX(_mm_xor_si128(v, _mm_set1_epi32(-1))); }

    // In-place versions
    inline mX& operator&=(const mX& other) { v = _mm_and_si128(v, other.v); return *this; }
    inline mX& operator|=(const mX& other) { v = _mm_or_si128(v, other.v); return *this; }
    inline mX& operator^=(const mX& other) { v = _mm_xor_si128(v, other.v); return *this; }

    // ------------------------------------------------------------------
    // Equality comparison
    // ------------------------------------------------------------------
    inline bool operator==(const mX& other) const {
        __m128i cmp = _mm_xor_si128(v, other.v);
        return _mm_testz_si128(cmp, cmp); // true if all bits are zero
    }
    inline bool operator!=(const mX& other) const { return !(*this == other); }
};

// ===============================================================
// set_zero: generic helper for arbitrary numeric types
// ===============================================================
//
// Works for scalars, __m128i, mX, and containers like std::vector<T>.
//
template<typename T>
inline void set_zero(T& x) {
    if constexpr (std::is_arithmetic_v<T>) {
        x = 0;
    } else if constexpr (std::is_same_v<T, __m128i>) {
        x = _mm_setzero_si128();
    } else if constexpr (std::is_same_v<T, mX>) {
        x = mX(_mm_setzero_si128());
    } else if constexpr (requires { x.size(); }) {
        for (auto& e : x) set_zero(e);
    } else {
        static_assert(!sizeof(T*), "set_zero: unsupported type");
    }
}

// ===============================================================
// Vector operations: elementwise arithmetic and dot products
// ===============================================================

// Elementwise vector addition
template<typename T>
static inline std::vector<T> operator^(const std::vector<T>& a,
                                       const std::vector<T>& b) {
    if (a.size() != b.size())
        throw std::runtime_error("operator+: size mismatch");
    std::vector<T> r(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        r[i] = a[i] ^ b[i];
    return r;
}

// Elementwise vector addition
template<typename T>
static inline std::vector<T> operator+(const std::vector<T>& a,
                                       const std::vector<T>& b) {
    if (a.size() != b.size())
        throw std::runtime_error("operator+: size mismatch");
    std::vector<T> r(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        r[i] = a[i] + b[i];
    return r;
}

// Elementwise vector subtraction
template<typename T>
static inline std::vector<T> operator-(const std::vector<T>& a,
                                       const std::vector<T>& b) {
    if (a.size() != b.size())
        throw std::runtime_error("operator-: size mismatch");
    std::vector<T> r(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        r[i] = a[i] - b[i];
    return r;
}

// Homogeneous dot product (vector<T> ⋅ vector<T>)
template<typename T>
static inline T operator*(const std::vector<T>& a,
                          const std::vector<T>& b) {

    if (a.size() != b.size())
        throw std::runtime_error("dot: size mismatch");
    T acc;
    set_zero(acc);
    for (size_t i = 0; i < a.size(); ++i)
        acc = acc + (a[i] * b[i]);
    return acc;
}

// Homogeneous dot product (vector<T> ⋅ vector<T>)
template<typename T>
static inline T operator&(const std::vector<T>& a,
                          const std::vector<T>& b)
{
    if (a.size() != b.size())
        throw std::runtime_error("dot: size mismatch");

    // Ensure acc is initialized to the additive identity
    T acc{};
    set_zero(acc);

    for (size_t i = 0; i < a.size(); ++i) {
        // Optional but strongly recommended sanity checks (remove in production)
        // assert(a[i] & one == a[i] when b[i] == one);
       // std::cout << a[i] << " & " << b[i] << " = " << (a[i] & b[i]) << std::endl;
       // acc ^= (a[i] * b[i]);
        acc ^= (a[i] * b[i]);
    }
    return acc;
}



// Mixed dot product: vector<mX> ⋅ vector<Scalar>
template<typename Scalar, typename = std::enable_if_t<std::is_arithmetic_v<Scalar>>>
static inline mX operator*(const std::vector<mX>& a, const std::vector<Scalar>& b) {
    if (a.size() != b.size())
        throw std::runtime_error("dot (mX,Scalar): size mismatch");
    mX acc(_mm_setzero_si128());
    for (size_t i = 0; i < a.size(); ++i)
        acc = acc + (a[i] * static_cast<uint8_t>(b[i]));
    return acc;
}

#include <immintrin.h>
#include <ostream>
#include <iomanip>

inline std::ostream& operator<<(std::ostream& os, const __m128i& v)
{
    alignas(16) uint8_t bytes[16];
    _mm_storeu_si128(reinterpret_cast<__m128i*>(bytes), v);

    os << "0x";
    for (int i = 15; i >= 0; --i) {  // big-endian display
        os << std::hex << std::setw(2) << std::setfill('0')
           << static_cast<int>(bytes[i]);
    }
    os << std::dec;
    return os;
}


// Symmetric version: vector<Scalar> ⋅ vector<mX>
template<typename Scalar, typename = std::enable_if_t<std::is_arithmetic_v<Scalar>>>
static inline mX operator*(const std::vector<Scalar>& a, const std::vector<mX>& b) {
    return b * a; // commutative
}

// ===============================================================
// Tag struct (used elsewhere for DPF key typing)
// ===============================================================
template<typename leaf_t, typename node_t, typename prgkey_t>
struct dpf_key_tag {};
