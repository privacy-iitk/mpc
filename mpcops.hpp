
#ifndef MPCOPS_HPP
#define MPCOPS_HPP
#pragma once
#include <bsd/stdlib.h>   // must come before any template using arc4random_buf
#include <cstdint>
#include <boost/asio.hpp>


#include <boost/asio.hpp>
//#include <boost/asio/experimental/awaitable_operators.hpp>
#include <iostream>
#include <random>
#include <stdexcept>
 
#include "prg.hpp"
using boost::asio::awaitable;
using boost::asio::use_awaitable;
//using namespace boost::asio::experimental::awaitable_operators;


#include "shares.hpp"
 
// Random buffer generator
// ---------------------------------------------------------------------
inline void arc4random_buf(void* buf, size_t n) {
    static thread_local std::mt19937_64 gen(std::random_device{}());
    uint8_t* p = static_cast<uint8_t*>(buf);
    for (size_t i = 0; i < n; ++i)
        p[i] = static_cast<uint8_t>(gen() & 0xFF);
}

// ---------------------------------------------------------------------
// Mock secure multiplication for additive shares (templated)
// ---------------------------------------------------------------------


template<typename T>
inline boost::asio::awaitable<T> mpc_dotproduct(
    const AdditiveShareVector<uint64_t> x,
    const AdditiveShareVector<uint64_t> y,
    Role role,
    NetPeer& self,
    NetPeer* peer_ptr = nullptr,
    NetPeer* second_peer = nullptr
) {
    static_assert(std::is_integral_v<T>, "mpc_dotproduct requires an integral type");
    T out = 0;
 
    // -------------------------------------------------------------------------
    // P2: Correlated randomness generation
    // -------------------------------------------------------------------------
    if (role == Role::P2) {
        if (!peer_ptr || !second_peer) co_return 0;
        NetPeer& p0 = *peer_ptr;
        NetPeer& p1 = *second_peer;

        AES_KEY key;
      //  std::cout << " ----->>> " << x.vals.size() << std::endl;

        std::vector<T> X0(x.vals.size()), X1(x.vals.size());
        std::vector<T> Y0(y.vals.size()), Y1(y.vals.size());

        __m128i seed_X0, seed_X1, seed_Y0, seed_Y1;
        arc4random_buf(&seed_X0, sizeof(__m128i));
        arc4random_buf(&seed_X1, sizeof(__m128i));
        arc4random_buf(&seed_Y0, sizeof(__m128i));
        arc4random_buf(&seed_Y1, sizeof(__m128i));

        crypto::fill_vector_with_prg(X0, key, seed_X0);
        crypto::fill_vector_with_prg(X1, key, seed_X1);
        crypto::fill_vector_with_prg(Y0, key, seed_Y0);
        crypto::fill_vector_with_prg(Y1, key, seed_Y1);

        T gamma0 = X0 * Y1; // vector-vector dotproduct 
        T gamma1 = X1 * Y0; // vector-vector dotproduct

        std::cout << "finished dotproduct" << std::endl;

        // co_await (
        //     (p0 << seed_X0) && (p0 << seed_Y0) && (p0 << gamma0) &&
        //     (p1 << seed_X1) && (p1 << seed_Y1) && (p1 << gamma1)
        // );

        co_await (p0 << seed_X0);
        co_await (p0 << seed_Y0);
        co_await (p0 << gamma0);

        co_await (p1 << seed_X1);
        co_await (p1 << seed_Y1);
        co_await (p1 << gamma1);

    }

    // -------------------------------------------------------------------------
    // P0: Receives seeds and computes blinded vectors
    // -------------------------------------------------------------------------
    else if (role == Role::P0) {
        if (!peer_ptr) co_return 0;
        NetPeer& peer = *peer_ptr;
        AES_KEY key;

        __m128i seed_X0, seed_Y0;
        T gamma0;
        //co_await ((self >> seed_X0) && (self >> seed_Y0) && (self >> gamma0));
        co_await (self >> seed_X0);
        co_await (self >> seed_Y0);
        co_await (self >> gamma0);


        std::vector<T> X_tilde(x.vals.size()), Y_tilde(x.vals.size());
        std::vector<T> X_blind(x.vals.size()), Y_blind(x.vals.size());

        crypto::fill_vector_with_prg(X_blind, key, seed_X0);
        crypto::fill_vector_with_prg(Y_blind, key, seed_Y0);

        // co_await (
        //     (peer << (x.vals + X_blind)) &&
        //     (peer << (y.vals + Y_blind)) &&
        //     (peer >> X_tilde) &&
        //     (peer >> Y_tilde)
        // );


            co_await (peer << (x.vals + X_blind));
            co_await (peer << (y.vals + Y_blind));
            co_await (peer >> X_tilde);
            co_await (peer >> Y_tilde);
        T z0 =  x.vals * (y.vals + Y_tilde) - (Y_blind * X_tilde) + gamma0;
        out = z0;
    }
    // -------------------------------------------------------------------------
    // P1: Receives seeds and computes blinded vectors
    // -------------------------------------------------------------------------
    else if (role == Role::P1) {
        if (!peer_ptr) co_return 0;
        NetPeer& peer = *peer_ptr;
        AES_KEY key;

        __m128i seed_X1, seed_Y1;
        T gamma1;
        
        //co_await ((self >> seed_X1) && (self >> seed_Y1) && (self >> gamma1));
        co_await (self >> seed_X1);
         co_await (self >> seed_Y1);
          co_await (self >> gamma1);

        std::vector<T> X_tilde(x.vals.size()), Y_tilde(x.vals.size());
        std::vector<T> X_blind(x.vals.size()), Y_blind(x.vals.size());

        crypto::fill_vector_with_prg(X_blind, key, seed_X1);
        crypto::fill_vector_with_prg(Y_blind, key, seed_Y1);

        // co_await (
        //     (peer << (x.vals + X_blind)) &&
        //     (peer << (y.vals + Y_blind)) &&
        //     (peer >> X_tilde) &&
        //     (peer >> Y_tilde)
        // );

  
        co_await (peer << (x.vals + X_blind));
        co_await (peer << (y.vals + Y_blind));
        co_await (peer >> X_tilde); 
        co_await (peer >> Y_tilde);


        T z1 = x.vals * (y.vals + Y_tilde) - (Y_blind * X_tilde) + gamma1;
 
        out = z1;
    }

    co_return out;
}

// -----------------------------------------------------------------------------
// Operator* for AdditiveShareVector<T> — triggers MPC dot product
// -----------------------------------------------------------------------------
template<typename T>
inline boost::asio::awaitable<T> operator*(
    const AdditiveShareVector<T>& x,
    const AdditiveShareVector<T>& y)
{
    if (!x.ctx)
        throw std::runtime_error("AdditiveShareVector missing MPC context");
    MPCContext* ctx = x.ctx;
    // Role role = x.ctx->role;
    // NetPeer& self = x.ctx->self;
    // NetPeer* peer0 = x.ctx->peer0;
    // NetPeer* peer1 = x.ctx->peer1;

    co_return co_await mpc_dotproduct<T>(x, y, ctx->role, ctx->self, ctx->peer0, ctx->peer1);
}



template<typename T>
inline boost::asio::awaitable<T> mpc_mul(
    const AShare<T>& x,
    const AShare<T>& y,
    Role role,
    NetPeer& self,
    NetPeer* peer_ptr = nullptr,
    NetPeer* second_peer = nullptr
) {
    static_assert(std::is_integral_v<T>, "mpc_mul requires an integral type");

    T out = 0;
    if (role == Role::P2) {
        if (!peer_ptr || !second_peer) co_return 0;
        NetPeer& p0 = *peer_ptr;
        NetPeer& p1 = *second_peer;

        T X0, Y0, X1, Y1, gamma0, gamma1;
        // arc4random_buf(&X0, sizeof(T));
        // arc4random_buf(&Y0, sizeof(T));
        // arc4random_buf(&X1, sizeof(T));
        // arc4random_buf(&Y1, sizeof(T));

        mpc_random(x.ctx, X0);
        mpc_random(x.ctx, Y0);
        mpc_random(x.ctx, X1);
        mpc_random(x.ctx, Y1);

        gamma0 = X0 * Y1;
        gamma1 = X1 * Y0;

        if (x.ctx->rand_mode == RandomnessMode::Online) {
       
        // Send concurrently
            // co_await (
            //     (p0 << X0) && (p1 << Y0) &&
            //     (p0 << X1) && (p1 << Y1) &&
            //     (p0 << gamma0) && (p1 << gamma1)
            // );
        co_await (p0 << X0);
        co_await (p1 << Y0);

        co_await (p0 << X1);
        co_await (p1 << Y1);

        co_await (p0 << gamma0);
        co_await (p1 << gamma1);

         std::cout << std::endl << "RandomnessMode::Online" << std::endl << std::endl;
       
        }

        if (x.ctx->rand_mode == RandomnessMode::Record) {

            std::cout << std::endl << 
                                "RandomnessMode::Record" << std::endl << std::endl;

            uint64_t idx = x.ctx->mul_ctr++;

            // Write P0 material
            x.ctx->mul_out_p0.write(reinterpret_cast<char*>(&idx), sizeof(idx));
            x.ctx->mul_out_p0.write(reinterpret_cast<char*>(&X0), sizeof(T));
            x.ctx->mul_out_p0.write(reinterpret_cast<char*>(&Y0), sizeof(T));
            x.ctx->mul_out_p0.write(reinterpret_cast<char*>(&gamma0), sizeof(T));

            // Write P1 material
            x.ctx->mul_out_p1.write(reinterpret_cast<char*>(&idx), sizeof(idx));
            x.ctx->mul_out_p1.write(reinterpret_cast<char*>(&X1), sizeof(T));
            x.ctx->mul_out_p1.write(reinterpret_cast<char*>(&Y1), sizeof(T));
            x.ctx->mul_out_p1.write(reinterpret_cast<char*>(&gamma1), sizeof(T));

            #ifdef VERBOSE 
                std::cout << "X0 = " << X0 << std::endl;
                std::cout << "Y0 = " << Y0 << std::endl;
                std::cout << "gamma0 = " << gamma0 << std::endl << std::endl;

                std::cout << "X1 = " << X1 << std::endl;
                std::cout << "Y1 = " << Y1 << std::endl;
                std::cout << "gamma1 = " << gamma1 << std::endl;
            #endif
        }        
    }

    else if (role == Role::P0 || role == Role::P1) {

        if (!peer_ptr) co_return 0;
        NetPeer& peer = *peer_ptr;

        const bool is_p0 = (role == Role::P0);

        T X, Y, gamma;
        uint64_t idx_file;

        // --------------------------------
        // Read preprocessing material
        // --------------------------------
        if (x.ctx->rand_mode == RandomnessMode::Online) {
            //co_await ((self >> X) && (self >> Y) && (self >> gamma));

            co_await (self >> X); 
            co_await (self >> Y);
            co_await  (self >> gamma);
        }
        else if (x.ctx->rand_mode == RandomnessMode::Replay) {

            std::ifstream& tape =
                is_p0 ? x.ctx->mul_in_p0 : x.ctx->mul_in_p1;

            tape.read(reinterpret_cast<char*>(&idx_file), sizeof(idx_file));

            uint64_t expected = x.ctx->mul_ctr++;

            if (idx_file != expected)
                throw std::runtime_error(
                    std::string("Replay desync in mpc_mul at ")
                    + (is_p0 ? "P0" : "P1")
                );

            tape.read(reinterpret_cast<char*>(&X), sizeof(T));
            tape.read(reinterpret_cast<char*>(&Y), sizeof(T));
            tape.read(reinterpret_cast<char*>(&gamma), sizeof(T));
            
            #ifdef VERBOSE
                std::cout << "index: " << idx_file << " =? " << expected << std::endl;    
                std::cout << "X = " << X << std::endl;
                std::cout << "Y = " << Y << std::endl;
                std::cout << "gamma = " << gamma << std::endl << std::endl;
            #endif
        }

        // --------------------------------
        // Local shares
        // --------------------------------
        const T x_i = x.val;
        const T y_i = y.val;

        T X_tilde, Y_tilde;

        // --------------------------------
        // Exchange masked values
        // --------------------------------
        // if (is_p0) {
        //     co_await (
        //         (peer << (x_i + X)) &&
        //         (peer << (y_i + Y)) &&
        //         (peer >> X_tilde) &&
        //         (peer >> Y_tilde)
        //     );
        // } else {
        //     co_await (
        //         (peer >> X_tilde) &&
        //         (peer >> Y_tilde) &&
        //         (peer << (x_i + X)) &&
        //         (peer << (y_i + Y))
        //     );
        // }

        if (is_p0) {
    co_await (peer << (x_i + X));
    co_await (peer << (y_i + Y));
    co_await (peer >> X_tilde);
    co_await (peer >> Y_tilde);
} else {
    co_await (peer >> X_tilde);
    co_await (peer >> Y_tilde);
    co_await (peer << (x_i + X));
    co_await (peer << (y_i + Y));
}


        // --------------------------------
        // Beaver reconstruction
        // --------------------------------
        out = x_i * (y_i + Y_tilde)
            - Y * X_tilde
            + gamma;
    }

    co_return out;
}

template<typename T>
inline boost::asio::awaitable<T> mpc_or(
    const XShare<T>& x,
    const XShare<T>& y,
    Role role,
    NetPeer& self,
    NetPeer* peer_ptr = nullptr,
    NetPeer* second_peer = nullptr
) {
    static_assert(std::is_integral_v<T>, "mpc_or requires an integral type");

    // Free XOR on shares
    XShare<T> xor_xy = x ^ y;

    // One MPC AND (uses operator*)
    XShare<T> and_xy = co_await (x * y);

    // OR = XOR ⊕ AND
    XShare<T> or_xy = xor_xy ^ and_xy;

    co_return or_xy.val;
}



// template<typename T>
// inline boost::asio::awaitable<T> mpc_and(
//     const XShare<T>& x, 
//     const XShare<T>& y,
//     Role role,    
//     NetPeer& self,
//     NetPeer* peer_ptr = nullptr,
//     NetPeer* second_peer = nullptr
// ) {
//     static_assert(std::is_integral_v<T>, "mpc_and requires an integral type");

//     T out = 0;

//     if (role == Role::P2) {
//         if (!peer_ptr || !second_peer) co_return 0;
//         NetPeer& p0 = *peer_ptr;
//         NetPeer& p1 = *second_peer;

//         T X0, Y0, X1, Y1, gamma0, gamma1;
//         arc4random_buf(&X0, sizeof(T));
//         arc4random_buf(&Y0, sizeof(T));
//         arc4random_buf(&X1, sizeof(T));
//         arc4random_buf(&Y1, sizeof(T));

//         gamma0 = X0 & Y1;
//         gamma1 = X1 & Y0;

//         // Parallel sends
//         co_await (
//             (p0 << X0) && (p0 << Y0) && (p0 << gamma0) &&
//             (p1 << X1) && (p1 << Y1) && (p1 << gamma1)
//         );
//     }

//     else if (role == Role::P0 || role == Role::P1) {
//         if (!peer_ptr) co_return 0;
//         NetPeer& peer = *peer_ptr;  // P1

//         T X0, Y0, gamma0;
//         co_await ((self >> X0) && (self >> Y0) && (self >> gamma0));

//         T x0 = x.val;
//         T y0 = y.val; 

//         T a0 = x0 ^ X0;
//         T b0 = y0 ^ Y0;

//         T a1, b1;
//         // Exchange masked values
//         co_await ((peer << a0) && (peer << b0) && (peer >> a1) && (peer >> b1));

//         T z0 = (x0 & (y0 ^ b1)) ^ (Y0 & a1) ^ gamma0;
//         out = z0;
//     }

//     // else if (role == Role::P1) {
//     //     if (!peer_ptr) co_return 0;
//     //     NetPeer& peer = *peer_ptr;  // P0

//     //     T X1, Y1, gamma1;
//     //     co_await ((self >> X1) && (self >> Y1) && (self >> gamma1));

//     //     T x1 = x.val;
//     //     T y1 = y.val;

//     //     T a0, b0;
        
//     //     T a1 = x1 ^ X1;
//     //     T b1 = y1 ^ Y1;
        
//     //     co_await ((peer >> a0) && (peer >> b0) && (peer << a1) && (peer << b1));

//     //     T z1 = (x1 & (y1 ^ b0)) ^ (Y1 & a0) ^ gamma1;
//     //     out = z1;
//     // }

//     co_return out;
// }

template<typename T>
inline boost::asio::awaitable<T> mpc_and(
    const XShare<T>& x,
    const XShare<T>& y,
    Role role,
    NetPeer& self,
    NetPeer* peer_ptr = nullptr,
    NetPeer* second_peer = nullptr
) {
    static_assert(std::is_integral_v<T>, "mpc_and requires an integral type");

    T out = 0;

    // --------------------------------------------------
    // P2: preprocessing / online randomness
    // --------------------------------------------------
    if (role == Role::P2) {
        if (!peer_ptr || !second_peer) co_return 0;

        NetPeer& p0 = *peer_ptr;
        NetPeer& p1 = *second_peer;

        T X0, Y0, X1, Y1, gamma0, gamma1;

        if (x.ctx->rand_mode == RandomnessMode::Online ||
            x.ctx->rand_mode == RandomnessMode::Record) {

            mpc_random(x.ctx, X0);
            mpc_random(x.ctx, Y0);
            mpc_random(x.ctx, X1);
            mpc_random(x.ctx, Y1);

            gamma0 = X0 & Y1;
            gamma1 = X1 & Y0;
        }

        if (x.ctx->rand_mode == RandomnessMode::Online) {
            // co_await (
            //     (p0 << X0) && (p0 << Y0) && (p0 << gamma0) &&
            //     (p1 << X1) && (p1 << Y1) && (p1 << gamma1)
            // );

            co_await (p0 << X0);
            co_await (p0 << Y0);
            co_await (p0 << gamma0);

            co_await (p1 << X1);
            co_await (p1 << Y1);
            co_await (p1 << gamma1);

        }
        else if (x.ctx->rand_mode == RandomnessMode::Record) {
            uint64_t idx = x.ctx->and_ctr++;

            x.ctx->and_out_p0.write(reinterpret_cast<char*>(&idx), sizeof(idx));
            x.ctx->and_out_p0.write(reinterpret_cast<char*>(&X0), sizeof(T));
            x.ctx->and_out_p0.write(reinterpret_cast<char*>(&Y0), sizeof(T));
            x.ctx->and_out_p0.write(reinterpret_cast<char*>(&gamma0), sizeof(T));

            x.ctx->and_out_p1.write(reinterpret_cast<char*>(&idx), sizeof(idx));
            x.ctx->and_out_p1.write(reinterpret_cast<char*>(&X1), sizeof(T));
            x.ctx->and_out_p1.write(reinterpret_cast<char*>(&Y1), sizeof(T));
            x.ctx->and_out_p1.write(reinterpret_cast<char*>(&gamma1), sizeof(T));
        }

        co_return 0;
    }

    // --------------------------------------------------
    // P0 / P1: unified branch
    // --------------------------------------------------
    else if (role == Role::P0 || role == Role::P1) {

        if (!peer_ptr) co_return 0;
        NetPeer& peer = *peer_ptr;

        const bool is_p0 = (role == Role::P0);

        T X, Y, gamma;
        uint64_t idx_file;

        // -------------------------------
        // Receive / replay preprocessing
        // -------------------------------
        if (x.ctx->rand_mode == RandomnessMode::Online) {
            
           // co_await ((self >> X) && (self >> Y) && (self >> gamma));

            co_await (self >> X);
            co_await (self >> Y);
            co_await (self >> gamma);


        }
        else if (x.ctx->rand_mode == RandomnessMode::Replay) {

            std::ifstream& tape =
                is_p0 ? x.ctx->and_in_p0 : x.ctx->and_in_p1;

            tape.read(reinterpret_cast<char*>(&idx_file), sizeof(idx_file));

            uint64_t expected = x.ctx->and_ctr++;
            if (idx_file != expected)
                throw std::runtime_error(
                    std::string("Replay desync in mpc_and at ")
                    + (is_p0 ? "P0" : "P1")
                );

            tape.read(reinterpret_cast<char*>(&X), sizeof(T));
            tape.read(reinterpret_cast<char*>(&Y), sizeof(T));
            tape.read(reinterpret_cast<char*>(&gamma), sizeof(T));


            std::cout << "X = " << X << std::endl 
                      << "Y = " << Y << std::endl 
                      << "gamma = " << gamma << std::endl;
        }

        // -------------------------------
        // Local shares
        // -------------------------------
        const T xi = x.val;
        const T yi = y.val;

        const T ai = xi ^ X;
        const T bi = yi ^ Y;

        T aj, bj;

        // -------------------------------
        // Exchange masked values
        // -------------------------------
        // if (is_p0) {
        //     co_await (
        //         (peer << ai) &&
        //         (peer << bi) &&
        //         (peer >> aj) &&
        //         (peer >> bj)
        //     );
        // } else {
        //     co_await (
        //         (peer >> aj) &&
        //         (peer >> bj) &&
        //         (peer << ai) &&
        //         (peer << bi)
        //     );
        // }
        if (is_p0) {
            co_await (peer << ai);
            co_await (peer << bi);
            co_await (peer >> aj);
            co_await (peer >> bj);
        } else {
            co_await (peer >> aj);
            co_await (peer >> bj);
            co_await (peer << ai);
            co_await (peer << bi);
        }

        // -------------------------------
        // AND reconstruction
        // -------------------------------
        out = (xi & (yi ^ bj)) ^ (Y & aj) ^ gamma;
    }

    co_return out;
}



// ---------------------------------------------------------------------
// operator* for AShare<T>
// ---------------------------------------------------------------------
template<typename T>
awaitable<AShare<T>> operator*(const AShare<T>& x, const AShare<T>& y) {
    if (!x.ctx || !y.ctx)
        throw std::runtime_error("AShare::operator* missing ctx");

    MPCContext* ctx = x.ctx;
    uint64_t result = co_await mpc_mul(
        static_cast<const AShare<uint64_t>&>(x),
        static_cast<const AShare<uint64_t>&>(y),
        ctx->role, ctx->self, ctx->peer0, ctx->peer1);
    co_return AShare<T>(static_cast<T>(result), ctx);
}

// ---------------------------------------------------------------------
// operator* for XShare<T>
// ---------------------------------------------------------------------
template<typename T>
awaitable<XShare<T>> operator*(const XShare<T>& x, const XShare<T>& y) {
    if (!x.ctx || !y.ctx)
        throw std::runtime_error("XShare::operator* missing ctx");

    MPCContext* ctx = x.ctx;
    T result = co_await mpc_and(x, y, ctx->role, ctx->self, ctx->peer0, ctx->peer1);
    co_return XShare<T>(result, ctx);
}

// ---------------------------------------------------------------------
// XOR-like arithmetic for XShare<T>
// ---------------------------------------------------------------------
template<typename T>
inline XShare<T> operator+(const XShare<T>& a, const XShare<T>& b) {
    if (!a.ctx || !b.ctx) throw std::runtime_error("XShare::operator+ missing ctx");
    return XShare<T>(a.val ^ b.val, a.ctx);
}
template<typename T>
inline XShare<T> operator-(const XShare<T>& a, const XShare<T>& b) {
    if (!a.ctx || !b.ctx) throw std::runtime_error("XShare::operator- missing ctx");
    return XShare<T>(a.val ^ b.val, a.ctx);
}
template<typename T>
inline XShare<T> operator^(const XShare<T>& a, const XShare<T>& b) {
    if (!a.ctx || !b.ctx) throw std::runtime_error("XShare::operator^ missing ctx");
    return XShare<T>(a.val ^ b.val, a.ctx);
}

template<typename T>
inline boost::asio::awaitable<XShare<T>>
operator|(const XShare<T>& x, const XShare<T>& y) {
    if (!x.ctx || !y.ctx)
        throw std::runtime_error("XShare::operator| missing ctx");

    MPCContext* ctx = x.ctx;

    T result = co_await mpc_or(
        x, y,
        ctx->role,
        ctx->self,
        ctx->peer0,
        ctx->peer1
    );

    co_return XShare<T>(result, ctx);
}

 
 
 inline boost::asio::awaitable<XShare<uint64_t>> mpc_eqz(
    const AShare<uint64_t>& x,
    Role role,
    NetPeer& self,
    NetPeer* peer_ptr,
    NetPeer* second_peer
) {
    

    // --------------------------------------------------
    // P2
    // --------------------------------------------------
    if (role == Role::P2) {
        if (!peer_ptr || !second_peer)
            co_return XShare<uint64_t>{0};

        NetPeer& p0 = *peer_ptr;
        NetPeer& p1 = *second_peer;

        uint64_t r, r0AS, r1AS, r0XS, r1XS;
        arc4random_buf(&r, sizeof(r));
        arc4random_buf(&r0AS, sizeof(r0AS));
        arc4random_buf(&r0XS, sizeof(r0XS));

        r1AS = r - r0AS;
        r1XS = r ^ r0XS;

        // Send masks
        // co_await (
        //     (p0 << r0AS) && (p1 << r1AS) &&
        //     (p0 << r0XS) && (p1 << r1XS)
        // );
        co_await (p0 << r0AS);
        co_await (p1 << r1AS);

        co_await (p0 << r0XS);
        co_await (p1 << r1XS);

        // --------------------------------------------------
        // Dummy OR-reduction to match P0/P1 control flow
        // --------------------------------------------------
        std::vector<XShare<uint64_t>> level;
        level.reserve(64);

        // Dummy leaves
        for (int i = 0; i < 64; i++)
            level.emplace_back(0, x.ctx);

        // Tree reduction (identical structure)
        while (level.size() > 1) {
            std::vector<XShare<uint64_t>> next;
            next.reserve((level.size() + 1) / 2);

            for (size_t i = 0; i < level.size(); i += 2) {
                if (i + 1 < level.size()) {
                    XShare<uint64_t> t =
                        co_await (level[i] | level[i + 1]);
                    next.push_back(t);
                } else {
                    next.push_back(level[i]);
                }
            }

            level.swap(next);
        }

        // Result meaningless for P2
 
        co_return XShare<uint64_t>{0, x.ctx};
    }

    // --------------------------------------------------
    // P0
    // --------------------------------------------------
    else if (role == Role::P0 || role == Role::P1) {
        if (!peer_ptr)
            co_return XShare<uint64_t>{0};

        NetPeer& peer = *peer_ptr;

        uint64_t rAS, rXS;
        //co_await ((self >> rAS) && (self >> rXS));
        co_await (self >> rAS);
        co_await (self >> rXS);
        
        uint64_t c = x.val + rAS;

        uint64_t c_recv;
        //co_await ((peer << c) && (peer >> c_recv));

        co_await (peer << c);
        co_await (peer >> c_recv);
        
        uint64_t c_rec = c + c_recv;

        uint64_t tmps[64];
        for (int i = 0; i < 64; i++) {
            uint64_t rXS_i = (rXS >> i) & 1ULL;
            tmps[i] = rXS_i;
            if(role == Role::P1) 
            {
                uint64_t c_rec_i = (c_rec >> i) & 1ULL; 
                tmps[i] = c_rec_i ^ rXS_i;
            }
        }

        std::vector<XShare<uint64_t>> level;
        level.reserve(64);

        // Initialize leaves
        for (int i = 0; i < 64; i++)
            level.emplace_back(tmps[i], x.ctx);

        // Tree reduction
        while (level.size() > 1) {
            std::vector<XShare<uint64_t>> next;
            next.reserve((level.size() + 1) / 2);

            for (size_t i = 0; i < level.size(); i += 2) {
                if (i + 1 < level.size()) {
                    XShare<uint64_t> t =
                        co_await (level[i] | level[i + 1]);
                    next.push_back(t);
                } else {
                    next.push_back(level[i]);
                }
            }

            level.swap(next);
        }

        // Final OR result
        XShare<uint64_t> acc = level[0];
        co_return acc;
    }

    // --------------------------------------------------
    // P1
    // --------------------------------------------------
    // else {
    //     if (!peer_ptr)
    //         co_return XShare<uint64_t>{0};

    //     NetPeer& peer = *peer_ptr;

    //     uint64_t rAS, rXS;
    //     co_await ((self >> rAS) && (self >> rXS));

    //     uint64_t c = x.val + rAS;

    //     uint64_t c_recv;
    //     co_await ((peer << c) && (peer >> c_recv));

    //     uint64_t c_rec = c + c_recv;

    //     uint64_t tmps[64];
    //     for (int i = 0; i < 64; i++) {
    //         uint64_t c_rec_i = (c_rec >> i) & 1ULL;
    //         uint64_t rXS_i  = (rXS >> i) & 1ULL;
    //         tmps[i] = c_rec_i ^ rXS_i;
    //     }

    //     std::vector<XShare<uint64_t>> level;
    //     level.reserve(64);

    //     // Initialize leaves
    //     for (int i = 0; i < 64; i++)
    //         level.emplace_back(tmps[i], x.ctx);

    //     // Tree reduction
    //     while (level.size() > 1) {
    //         std::vector<XShare<uint64_t>> next;
    //         next.reserve((level.size() + 1) / 2);

    //         for (size_t i = 0; i < level.size(); i += 2) {
    //             if (i + 1 < level.size()) {
    //                 XShare<uint64_t> t =
    //                     co_await (level[i] | level[i + 1]);
    //                 next.push_back(t);
    //             } else {
    //                 next.push_back(level[i]);
    //             }
    //         }

    //         level.swap(next);
    //     }

    //     // Final OR result
    //     XShare<uint64_t> acc = level[0];
    //     co_return acc;
    // }
}

inline boost::asio::awaitable<XShare<uint64_t>>
operator==(const AShare<uint64_t>& x, uint64_t zero) {
    if (zero != 0)
        throw std::runtime_error("operator== only supports == 0 for AShare");

    if (!x.ctx)
        throw std::runtime_error("AShare::operator== missing ctx");

    MPCContext* ctx = x.ctx;

    co_return co_await mpc_eqz(
        x,
        ctx->role,
        ctx->self,
        ctx->peer0,
        ctx->peer1
    );
}


#endif