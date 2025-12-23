#pragma once
#include <cstdint>
#include <vector>
#include <random>
#include <iostream>
#include <type_traits>
#include <utility>
#include <boost/asio/awaitable.hpp>
#include "common.hpp"   // assumes Role, NetPeer and async operator<< / operator>> are defined here

// --------------------- MPC Context ---------------------
struct MPCContext {
    Role role;
    NetPeer& self;
    NetPeer* peer0 = nullptr;
    NetPeer* peer1 = nullptr;
};

template<typename T>
struct XShare {
    static_assert(std::is_integral<T>::value, "AShare requires an integral type");

    T val{};
    MPCContext* ctx = nullptr;  // optional pointer to MPC context

    // Constructors
    XShare() = default;
    explicit XShare(T v, MPCContext* c = nullptr) : val(v), ctx(c) {}

    // Local XOR operations
    XShare operator^(const XShare& other) const {
        return XShare(static_cast<T>(val ^ other.val), ctx);
    }

    XShare& operator^=(const XShare& other) {
        val = static_cast<T>(val ^ other.val);
        return *this;
    }

    inline XShare operator+(const XShare& other) const {
        return XShare(val ^ other.val);   // XOR addition
    }

    inline XShare operator-(const XShare& other) const {
        return XShare(val ^ other.val);   // same as XOR
    }

    inline XShare operator*(const XShare& other) const {
        return XShare(val & other.val);   // bitwise AND as multiplication
    }

    // Output stream for debugging
    friend std::ostream& operator<<(std::ostream& os, const XShare& s) {
        os << s.val;
        return os;
    }
};


template<typename T>
struct AShare {
    static_assert(std::is_integral<T>::value, "AShare requires an integral type");

    T val{};
    MPCContext* ctx = nullptr;  // optional pointer to MPC context

    // Constructors
    AShare() = default;
    explicit AShare(T v, MPCContext* c = nullptr) : val(v), ctx(c) {}

    // Local XOR operations
    AShare operator+(const AShare& other) const {
        return AShare(static_cast<T>(val + other.val), ctx);
    }

    AShare& operator+=(const AShare& other) {
        val = static_cast<T>(val + other.val);
        return *this;
    }

    // Output stream for debugging
    friend std::ostream& operator<<(std::ostream& os, const AShare& s) {
        os << s.val;
        return os;
    }
};

// --------------------- Vector XOR Share ---------------------
template<typename T>
struct XorShareVector {
    static_assert(std::is_integral<T>::value, "XorShareVector requires an integral type");

    std::vector<T> vals;

    XorShareVector() = default;
    explicit XorShareVector(std::vector<T> v) : vals(std::move(v)) {}

    XorShareVector operator^(const XorShareVector& other) const {
        if (vals.size() != other.vals.size())
            throw std::runtime_error("XorShareVector: size mismatch in ^ operator");
        std::vector<T> out(vals.size());
        for (size_t i = 0; i < vals.size(); ++i)
            out[i] = static_cast<T>(vals[i] ^ other.vals[i]);
        return XorShareVector(std::move(out));
    }

    friend std::ostream& operator<<(std::ostream& os, const XorShareVector& s) {
        os << "[";
        for (size_t i = 0; i < s.vals.size(); ++i) {
            os << s.vals[i];
            if (i + 1 < s.vals.size()) os << ", ";
        }
        os << "]";
        return os;
    }
};

// --------------------- Random Share Generation ---------------------
inline uint64_t random_uint64() {
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    static std::uniform_int_distribution<uint64_t> dist;
    return dist(gen);
}

template<typename T>
inline std::pair<AShare<T>, AShare<T>> share_secret(T secret) {
    static_assert(std::is_integral<T>::value, "share_secret requires an integral type");
    // generate a random value of T (cast from uint64_t)
    uint64_t r = random_uint64();
    T s0 = static_cast<T>(r);
    T s1 = static_cast<T>(s0 ^ secret);
    return { AShare<T>(s0), AShare<T>(s1) };
}

// --------------------- Network-based Reconstruction ---------------------
template<typename T>
boost::asio::awaitable<T> reconstruct_remote(NetPeer& peer, const AShare<T>& my_share) {
    static_assert(std::is_integral<T>::value, "reconstruct_remote requires an integral type");

    // Send your share to the other NetPeer
    co_await (peer << my_share.val);

    // Receive their share
    T other_val{};
    co_await (peer >> other_val);

    T secret = static_cast<T>(my_share.val + other_val);
    co_return secret;
}

template<typename T>
boost::asio::awaitable<T> reconstruct_remote(NetPeer& peer, const XShare<T>& my_share) {
    static_assert(std::is_integral<T>::value, "reconstruct_remote requires an integral type");

    // Send your share to the other NetPeer
    co_await (peer << my_share.val);

    // Receive their share
    T other_val{};
    co_await (peer >> other_val);

    T secret = static_cast<T>(my_share.val ^ other_val);
    co_return secret;
}


// Reconstruct for vector XOR shares (send vector and receive peer's vector)
template<typename T>
boost::asio::awaitable<std::vector<T>> reconstruct_remote_vector(NetPeer& peer, const XorShareVector<T>& my_share) {
    static_assert(std::is_integral<T>::value, "reconstruct_remote_vector requires an integral type");

    // Send your vector share
    co_await (peer << my_share.vals);

    // Receive the peer's vector share
    std::vector<T> other_vals;
    co_await (peer >> other_vals);

    if (other_vals.size() != my_share.vals.size())
        throw std::runtime_error("reconstruct_remote_vector: size mismatch");

    // XOR element-wise
    std::vector<T> secret(other_vals.size());
    for (size_t i = 0; i < other_vals.size(); ++i)
        secret[i] = static_cast<T>(my_share.vals[i] ^ other_vals[i]);

    co_return secret;
}
