#pragma once

#include "dpf/dpf.h"

#include <boost/asio.hpp>
#include <boost/asio/awaitable.hpp>
#include <boost/asio/use_awaitable.hpp>
#include <boost/asio/experimental/awaitable_operators.hpp>
#include <boost/asio/steady_timer.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/detached.hpp>

#include <vector>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <utility>

using boost::asio::awaitable;
using boost::asio::use_awaitable;
using boost::asio::async_write;
using boost::asio::async_read;
using boost::asio::buffer;
using boost::asio::ip::tcp;
#include <immintrin.h>
 
#pragma once

#include <immintrin.h>
#include <cstdint>
#include <vector>
#include <type_traits>
#include <stdexcept>

// ===============================================================
// 128-bit modular arithmetic (mod 2^128) over __m128i
// ===============================================================

static inline __m128i mul128_mod2_128(__m128i a, __m128i b) {
    uint64_t a_lo = _mm_cvtsi128_si64(a);
    uint64_t a_hi = _mm_extract_epi64(a, 1);
    uint64_t b_lo = _mm_cvtsi128_si64(b);
    uint64_t b_hi = _mm_extract_epi64(b, 1);

    __extension__ typedef unsigned __int128 u128;

    u128 lo   = (u128)a_lo * b_lo;
    u128 mid1 = (u128)a_lo * b_hi;
    u128 mid2 = (u128)a_hi * b_lo;

    // ignore high * high because mod 2^128
    u128 r = lo + ((mid1 + mid2) << 64);

    uint64_t r_lo = (uint64_t)r;
    uint64_t r_hi = (uint64_t)(r >> 64);

    return _mm_set_epi64x(r_hi, r_lo);
}

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
// mX: 128-bit arithmetic element with operator overloads
// ===============================================================

struct mX {
    __m128i v;

    mX() = default;
    explicit mX(__m128i x) : v(x) {}
    mX(uint64_t hi, uint64_t lo) { v = _mm_set_epi64x(hi, lo); }

    // Accessors
    uint64_t lo() const { return _mm_cvtsi128_si64(v); }
    uint64_t hi() const { return _mm_extract_epi64(v, 1); }

    // Basic arithmetic
    inline mX operator*(const mX& other) const {
        return mX(mul128_mod2_128(v, other.v));
    }
    inline mX operator+(const mX& other) const {
        return mX(add128_mod2_128(v, other.v));
    }
    inline mX operator-(const mX& other) const {
        return mX(sub128_mod2_128(v, other.v));
    }

    // ===========================================================
    // Mixed scalar (uint8_t) operations — correctly promote to 128-bit
    // ===========================================================
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

    // Friend overloads for commutativity (scalar on left)
    friend inline mX operator*(uint8_t s, const mX& x) { return x * s; }
    friend inline mX operator+(uint8_t s, const mX& x) { return x + s; }
    friend inline mX operator-(uint8_t s, const mX& x) {
        __m128i a = _mm_set_epi64x(0, (uint64_t)s);
        return mX(sub128_mod2_128(a, x.v));
    }
};

// ===============================================================
// set_zero: generic template utility for arithmetic types
// ===============================================================

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
// Vector operations (elementwise arithmetic and dot product)
// ===============================================================

template<typename T>
static inline std::vector<T> operator+(const std::vector<T>& a,
                                       const std::vector<T>& b) {
    size_t n = a.size();
    std::vector<T> r(n);
    for (size_t i = 0; i < n; ++i)
        r[i] = a[i] + b[i];
    return r;
}

template<typename T>
static inline std::vector<T> operator-(const std::vector<T>& a,
                                       const std::vector<T>& b) {
    size_t n = a.size();
    std::vector<T> r(n);
    for (size_t i = 0; i < n; ++i)
        r[i] = a[i] - b[i];
    return r;
}

// dot product
template<typename T>
static inline T operator*(const std::vector<T>& a,
                          const std::vector<T>& b) {
    size_t n = a.size();
    if (n != b.size())
        throw std::runtime_error("dot: size mismatch");

    T acc;
    set_zero(acc);

    for (size_t i = 0; i < n; ++i)
        acc = acc + (a[i] * b[i]);

    return acc;
}

// explicit dotproduct function (same as operator*)
static inline mX dotproduct(const std::vector<mX>& a,
                            const std::vector<mX>& b) {
    return a * b;
}

// -----------------------------------------------------------------------------
// Tag type to select DPF operator overloads
// -----------------------------------------------------------------------------
// Tag struct for DPF key operations
template<typename leaf_t, typename node_t, typename prgkey_t>
struct dpf_key_tag {};


// -----------------------------------------------------------------------------
// Role and NetPeer definitions
// -----------------------------------------------------------------------------
enum class Role { P0, P1, P2 };

struct NetPeer {
    Role role;
    tcp::socket sock;

    NetPeer(Role r, tcp::socket&& s) : role(r), sock(std::move(s)) {}

    // -------------------------
    // Async send for primitive types
    // -------------------------
    template<typename T>
    awaitable<void> send(const T& value) {
        co_await async_write(sock, buffer(&value, sizeof(T)), use_awaitable);
        co_return;
    }

    // Async send for vectors
    template<typename T>
    awaitable<void> send(const std::vector<T>& vec) {
        uint64_t size = vec.size();
        co_await send(size); // send size first
        if (!vec.empty()) {
            co_await async_write(sock, buffer(vec.data(), sizeof(T) * vec.size()), use_awaitable);
        }
        co_return;
    }

    // -------------------------
    // Async receive for primitive types
    // -------------------------
    template<typename T>
    awaitable<T> recv() {
        T value;
        co_await async_read(sock, buffer(&value, sizeof(T)), use_awaitable);
        co_return value;
    }

    // Async receive for vectors
    template<typename T>
    awaitable<std::vector<T>> recv_vector() {
        uint64_t size = co_await recv<uint64_t>();
        std::vector<T> vec(size);
        if (size > 0) {
            co_await async_read(sock, buffer(vec.data(), sizeof(T) * size), use_awaitable);
        }
        co_return vec;
    }

    // -------------------------
    // Async send/recv for DPF keys
    // -------------------------
    template<typename leaf_t, typename node_t, typename prgkey_t>
    awaitable<void> send(const dpf::dpf_key<leaf_t, node_t, prgkey_t>& key) {
        auto buf = dpf::serialize_dpf_key<leaf_t, node_t, prgkey_t>(key);
        uint64_t len = buf.size();
        co_await send(len); // send size first
        if (!buf.empty()) {
            co_await async_write(sock, buffer(buf.data(), buf.size()), use_awaitable);
        }
        co_return;
    }

// -------------------------
// Async receive for DPF keys (fixed for non-trivial dpf_key)
// -------------------------
template<typename leaf_t, typename node_t, typename prgkey_t>
awaitable<dpf::dpf_key<leaf_t, node_t, prgkey_t>> recv_dpf_key() {
    // First, receive the size of the serialized buffer
    uint64_t len = co_await recv<uint64_t>();
    std::vector<uint8_t> buf(len);
    if (len > 0) {
        co_await async_read(sock, buffer(buf.data(), len), use_awaitable);
    }

    // Properly deserialize without default construction
    const uint8_t* data = buf.data();
    size_t offset = 0;

    // Deserialize nitems
    size_t nitems = *reinterpret_cast<const size_t*>(data + offset);
    offset += sizeof(size_t);

    // Deserialize root
    node_t root;
    std::memcpy(&root, data + offset, sizeof(node_t));
    offset += sizeof(node_t);

    // Deserialize cw vector
    size_t cw_size = *reinterpret_cast<const size_t*>(data + offset);
    offset += sizeof(size_t);
    std::vector<node_t> cw(cw_size);
    if (cw_size > 0) {
        std::memcpy(cw.data(), data + offset, cw_size * sizeof(node_t));
        offset += cw_size * sizeof(node_t);
    }

    // Deserialize finalizer
    leaf_t finalizer;
    std::memcpy(&finalizer, data + offset, sizeof(leaf_t));
    offset += sizeof(leaf_t);

    // Deserialize PRG key
    prgkey_t prgkey;
    std::memcpy(&prgkey, data + offset, sizeof(prgkey_t));

    // Call the proper constructor
    co_return dpf::dpf_key<leaf_t, node_t, prgkey_t>(nitems, root, cw, finalizer, prgkey);
}

// -------------------------
// Operator overload for dpf_key using tag
// -------------------------



    // -------------------------
    // Operator overloads for primitives and vectors
    // -------------------------
    template<typename T>
    awaitable<void> operator<<(const T& val) { co_await send(val); }

    template<typename T>
    awaitable<void> operator>>(T& out) { out = co_await recv<T>(); }

    template<typename T>
    awaitable<void> operator<<(const std::vector<T>& vec) { co_await send(vec); }

    template<typename T>
    awaitable<void> operator>>(std::vector<T>& out) { out = co_await recv_vector<T>(); }

    // =====================================================
    // Overload << and >> for DPF keys
    // =====================================================


};
