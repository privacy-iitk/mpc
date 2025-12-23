#pragma once
#include<fstream>
#include <cstdint>
#include <vector>
#include <random>
#include <iostream>
#include <type_traits>
#include <utility>
#include <boost/asio/awaitable.hpp>

#include "types.hpp"
#include "network.hpp"

using boost::asio::awaitable;

enum class RandomnessMode {
    Online,        // arc4random_buf
    Record,        // preprocessing (-p)
    Replay         // online MPC (replay from disk)
};


// MPCContext binds a role with communication peers
struct MPCContext {

    Role role;
    NetPeer& self;       // this NetPeer
    NetPeer* peer0;      // one peer
    NetPeer* peer1;      // optional second peer (for P2)

    // ---------------- NEW ----------------
    RandomnessMode rand_mode = RandomnessMode::Online;

    // Used only when rand_mode == Record
    std::ofstream rand_out;

    // Used only when rand_mode == Replay
    std::ifstream rand_in;
    // -------------------------------------

    std::ofstream mul_out_p0;
    std::ofstream mul_out_p1;
    std::ifstream mul_in_p0;
    std::ifstream mul_in_p1;
        std::ofstream and_out_p0;
    std::ofstream and_out_p1;
    std::ifstream and_in_p0;
    std::ifstream and_in_p1;
    uint64_t       mul_ctr, and_ctr;

    MPCContext(Role r,
               NetPeer& s,
               NetPeer* p0 = nullptr,
               NetPeer* p1 = nullptr)
        : role(r), self(s), peer0(p0), peer1(p1) {}
};


template<typename T>
struct XShare {
    static_assert(std::is_integral_v<T>, "XShare requires integral type");
    T val{};
    MPCContext* ctx = nullptr;
    XShare(T v, MPCContext* c = nullptr) : val(v), ctx(c) {}
    inline XShare operator^(const XShare& o) const { return XShare(static_cast<T>(val ^ o.val), ctx); }
    inline XShare& operator^=(const XShare& o) { val ^= o.val; return *this; }
    inline XShare operator+(const XShare& o) const { return *this ^ o; }
    inline XShare operator-(const XShare& o) const { return *this ^ o; }
     friend std::ostream& operator<<(std::ostream& os, const XShare& s) { os << s.val; return os; }
};

// Additive share
template<typename T>
struct AShare {
    static_assert(std::is_integral_v<T>, "AShare requires integral type");
    T val{};
    MPCContext* ctx = nullptr;
    AShare(T v, MPCContext* c = nullptr) : val(v), ctx(c) {}
    //explicit AShare(T v, MPCContext* c = nullptr) : val(v), ctx(c) {}
    inline AShare operator+(const AShare& o) const { return AShare(static_cast<T>(val + o.val), ctx); }
    inline AShare operator-(const AShare& o) const { return AShare(static_cast<T>(val - o.val), ctx); }
    inline AShare& operator+=(const AShare& o) { val = static_cast<T>(val + o.val); return *this; }
    friend std::ostream& operator<<(std::ostream& os, const AShare& s) { os << s.val; return os; }
};


// XOR share
// template<typename T>
// struct XShare {
//     static_assert(std::is_integral_v<T>, "XShare requires integral type");
//     T val{};
//     MPCContext* ctx = nullptr;
//     XShare() = default;
//     explicit XShare(T v, MPCContext* c = nullptr) : val(v), ctx(c) {}
//     inline XShare operator^(const XShare& o) const { return XShare(static_cast<T>(val ^ o.val), ctx); }
//     inline XShare& operator^=(const XShare& o) { val ^= o.val; return *this; }
//     inline XShare operator+(const XShare& o) const { return *this ^ o; }
//     inline XShare operator-(const XShare& o) const { return *this ^ o; }
//     inline XShare operator*(const XShare& o) const { return XShare(static_cast<T>(val & o.val), ctx); }
//     friend std::ostream& operator<<(std::ostream& os, const XShare& s) { os << s.val; return os; }
// };

// // Additive share
// template<typename T>
// struct AShare {
//     static_assert(std::is_integral_v<T>, "AShare requires integral type");
//     T val{};
//     MPCContext* ctx = nullptr;
//     AShare() = default;
//     explicit AShare(T v, MPCContext* c = nullptr) : val(v), ctx(c) {}
//     inline AShare operator+(const AShare& o) const { return AShare(static_cast<T>(val + o.val), ctx); }
//     inline AShare operator-(const AShare& o) const { return AShare(static_cast<T>(val - o.val), ctx); }
//     inline AShare& operator+=(const AShare& o) { val = static_cast<T>(val + o.val); return *this; }
//     friend std::ostream& operator<<(std::ostream& os, const AShare& s) { os << s.val; return os; }
// };

// vector XOR shares
template<typename T>
struct XorShareVector {
    static_assert(std::is_integral_v<T>, "XorShareVector requires integral type");
    std::vector<T> vals;
    XorShareVector() = default;
    explicit XorShareVector(std::vector<T> v) : vals(std::move(v)) {}
    inline XorShareVector operator^(const XorShareVector& o) const {
        if (vals.size() != o.vals.size()) throw std::runtime_error("XorShareVector size mismatch");
        std::vector<T> out(vals.size());
        for (size_t i = 0; i < vals.size(); ++i) out[i] = static_cast<T>(vals[i] ^ o.vals[i]);
        return XorShareVector(std::move(out));
    }
    friend std::ostream& operator<<(std::ostream& os, const XorShareVector& s) {
        os << "["; for (size_t i = 0; i < s.vals.size(); ++i) { os << s.vals[i]; if (i+1<s.vals.size()) os << ", "; } os << "]"; return os;
    }
};

// template<typename T>
// struct AdditiveShareVector {
//     static_assert(std::is_integral_v<T>, "AdditiveShareVector requires integral type");

//     std::vector<T> vals;

//     AdditiveShareVector() = default;
//     //explicit AdditiveShareVector(std::vector<T> v) : vals(std::move(v)) {}

//     // Element-wise addition of additive shares
//     inline AdditiveShareVector operator+(const AdditiveShareVector& o) const {
//         if (vals.size() != o.vals.size())
//             throw std::runtime_error("AdditiveShareVector size mismatch");
//         std::vector<T> out(vals.size());
//         for (size_t i = 0; i < vals.size(); ++i)
//             out[i] = static_cast<T>(vals[i] + o.vals[i]);
//         return AdditiveShareVector(std::move(out));
//     }

//     // Element-wise subtraction
//     inline AdditiveShareVector operator-(const AdditiveShareVector& o) const {
//         if (vals.size() != o.vals.size())
//             throw std::runtime_error("AdditiveShareVector size mismatch");
//         std::vector<T> out(vals.size());
//         for (size_t i = 0; i < vals.size(); ++i)
//             out[i] = static_cast<T>(vals[i] - o.vals[i]);
//         return AdditiveShareVector(std::move(out));
//     }

//     // In-place addition
//     inline AdditiveShareVector& operator+=(const AdditiveShareVector& o) {
//         if (vals.size() != o.vals.size())
//             throw std::runtime_error("AdditiveShareVector size mismatch");
//         for (size_t i = 0; i < vals.size(); ++i)
//             vals[i] = static_cast<T>(vals[i] + o.vals[i]);
//         return *this;
//     }

//     // In-place subtraction
//     inline AdditiveShareVector& operator-=(const AdditiveShareVector& o) {
//         if (vals.size() != o.vals.size())
//             throw std::runtime_error("AdditiveShareVector size mismatch");
//         for (size_t i = 0; i < vals.size(); ++i)
//             vals[i] = static_cast<T>(vals[i] - o.vals[i]);
//         return *this;
//     }

//     // Print nicely
//     friend std::ostream& operator<<(std::ostream& os, const AdditiveShareVector& s) {
//         os << "[";
//         for (size_t i = 0; i < s.vals.size(); ++i) {
//             os << s.vals[i];
//             if (i + 1 < s.vals.size()) os << ", ";
//         }
//         os << "]";
//         return os;
//     }
// };
template<typename T>
struct AdditiveShareVector {
    static_assert(std::is_integral_v<T>, "AdditiveShareVector requires integral type");

    std::vector<T> vals;
    MPCContext* ctx = nullptr;

    AdditiveShareVector() = default;
    explicit AdditiveShareVector(std::vector<T> v, MPCContext* c = nullptr)
        : vals(std::move(v)), ctx(c) {}

    inline AdditiveShareVector operator+(const AdditiveShareVector& o) const {
        if (vals.size() != o.vals.size())
            throw std::runtime_error("AdditiveShareVector size mismatch");
        std::vector<T> out(vals.size());
        for (size_t i = 0; i < vals.size(); ++i)
            out[i] = static_cast<T>(vals[i] + o.vals[i]);
        return AdditiveShareVector(std::move(out), ctx);
    }

    inline AdditiveShareVector operator-(const AdditiveShareVector& o) const {
        if (vals.size() != o.vals.size())
            throw std::runtime_error("AdditiveShareVector size mismatch");
        std::vector<T> out(vals.size());
        for (size_t i = 0; i < vals.size(); ++i)
            out[i] = static_cast<T>(vals[i] - o.vals[i]);
        return AdditiveShareVector(std::move(out), ctx);
    }

    friend std::ostream& operator<<(std::ostream& os, const AdditiveShareVector& s) {
        os << "[";
        for (size_t i = 0; i < s.vals.size(); ++i) {
            os << s.vals[i];
            if (i + 1 < s.vals.size()) os << ", ";
        }
        os << "]";
        return os;
    }
};



// portable RNG
inline uint64_t random_u64() {
    static thread_local std::mt19937_64 g(std::random_device{}());
    static std::uniform_int_distribution<uint64_t> dist;
    return dist(g);
}

// additive sharing
template<typename T>
inline std::pair<AShare<T>, AShare<T>> share_secret_additive(T secret) {
    static_assert(std::is_integral_v<T>, "share_secret_additive requires integral type");
    T s0 = static_cast<T>(random_u64());
    T s1 = static_cast<T>(secret - s0);
    return { AShare<T>(s0), AShare<T>(s1) };
}

// xor sharing
template<typename T>
inline std::pair<XShare<T>, XShare<T>> share_secret_xor(T secret) {
    static_assert(std::is_integral_v<T>, "share_secret_xor requires integral type");
    T s0 = static_cast<T>(random_u64());
    T s1 = static_cast<T>(s0 ^ secret);
    return { XShare<T>(s0), XShare<T>(s1) };
}

// network-based reconstruction (AShare)
template<typename T>
awaitable<T> reconstruct_remote(NetPeer& peer, const AShare<T>& my_share) {
    //std::cout << "reconstruct_remote" << std::endl;
    static_assert(std::is_integral_v<T>, "reconstruct_remote requires integral type");
    co_await (peer << my_share.val);
    T other = 0;
    co_await (peer >> other);
    co_return static_cast<T>(my_share.val + other);
}

// reconstruction (XShare)
template<typename T>
awaitable<T> reconstruct_remote(NetPeer& peer, const XShare<T>& my_share) {
  //  std::cout << "reconstruct remote (AND Shares!)" << std::endl;
    static_assert(std::is_integral_v<T>, "reconstruct_remote requires integral type");
    co_await (peer << my_share.val);
    T other = 0;
    co_await (peer >> other);
    co_return static_cast<T>(my_share.val ^ other);
}

// vector reconstruct (XorShareVector)
template<typename T>
awaitable<std::vector<T>> reconstruct_remote_vector(NetPeer& peer, const XorShareVector<T>& my_share) {
    static_assert(std::is_integral_v<T>, "reconstruct_remote_vector requires integral type");
    co_await (peer << my_share.vals);
    std::vector<T> other = co_await peer.recv_vector<T>();
    if (other.size() != my_share.vals.size()) throw std::runtime_error("reconstruct_remote_vector: size mismatch");
    std::vector<T> out(other.size());
    for (size_t i = 0; i < out.size(); ++i) out[i] = static_cast<T>(my_share.vals[i] ^ other[i]);
    co_return out;
}

template<typename T>
inline void mpc_random(MPCContext* ctx, T& out) {
    static_assert(std::is_trivially_copyable_v<T>,
                  "mpc_random requires trivially copyable types");

    if (ctx->rand_mode == RandomnessMode::Online) {
        arc4random_buf(&out, sizeof(T));
    }
    else if (ctx->rand_mode == RandomnessMode::Record) {
        arc4random_buf(&out, sizeof(T));
        ctx->rand_out.write(reinterpret_cast<const char*>(&out), sizeof(T));
        if (!ctx->rand_out)
            throw std::runtime_error("Failed to record randomness");
    }
    else { // Replay (will be used later)
        ctx->rand_in.read(reinterpret_cast<char*>(&out), sizeof(T));
        if (!ctx->rand_in)
            throw std::runtime_error("Out of preprocessing randomness");
    }
}
