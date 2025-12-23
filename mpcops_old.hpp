#pragma once
#include <cstdint>
#include <iostream>
#include <boost/asio/awaitable.hpp>
#include "shares.hpp"   // contains AShare<T>, MPCContext, Role, NetPeer, etc.
#include "common.hpp"   // if needed for arc4random_buf, operator<< / >> on NetPeer, etc.

template<typename leaf_t, typename node_t, typename prgkey_t>
std::vector<uint8_t> serialize_dpf_key(const dpf_key<leaf_t, node_t, prgkey_t>& key)
{
    static_assert(std::is_trivially_copyable_v<leaf_t>, "leaf_t must be trivially copyable for this serializer.");
    static_assert(std::is_trivially_copyable_v<node_t>, "node_t must be trivially copyable for this serializer.");
    std::vector<uint8_t> buf;
    auto append = [&](const void* data, size_t len) {
        const uint8_t* p = reinterpret_cast<const uint8_t*>(data);
        buf.insert(buf.end(), p, p + len);
    };

    // nitems and root
    uint64_t nitems_net = static_cast<uint64_t>(key.nitems);
    append(&nitems_net, sizeof(nitems_net));
    append(&key.root, sizeof(node_t));

    // cw length and contents
    uint64_t cw_size = static_cast<uint64_t>(key.cw.size());
    append(&cw_size, sizeof(cw_size));
    if (cw_size > 0) {
        append(key.cw.data(), sizeof(std::array<node_t,2>) * cw_size);
    }

    // finalizer array
    append(&key.finalizer, sizeof(typename dpf_key<leaf_t,node_t,prgkey_t>::finalizer_t));

    return buf;
}

template<typename leaf_t, typename node_t, typename prgkey_t>
dpf_key<leaf_t, node_t, prgkey_t> deserialize_dpf_key(const uint8_t* data, size_t len)
{
    dpf_key<leaf_t, node_t, prgkey_t> key;
    size_t offset = 0;
    auto read = [&](void* dst, size_t n) {
        if (offset + n > len) throw std::runtime_error("deserialize buffer overrun");
        std::memcpy(dst, data + offset, n);
        offset += n;
    };

    uint64_t nitems_net = 0;
    read(&nitems_net, sizeof(nitems_net));
    key.nitems = static_cast<size_t>(nitems_net);

    read(&key.root, sizeof(node_t));

    uint64_t cw_size = 0;
    read(&cw_size, sizeof(cw_size));
    key.cw.resize(static_cast<size_t>(cw_size));
    if (cw_size > 0) {
        read(key.cw.data(), sizeof(std::array<node_t,2>) * cw_size);
    }

    read(&key.finalizer, sizeof(typename dpf_key<leaf_t,node_t,prgkey_t>::finalizer_t));

    return key;
}






// --------------------- mpc_and (unchanged logic) ---------------------
// keep as you already wrote it; slightly changed signature to accept AShare<T> by const ref
// and kept the same behavior. (If you already have it in another .cpp/.h, keep single definition.)
inline boost::asio::awaitable<uint64_t> mpc_mul(
    const AShare<uint64_t>& x,
    const AShare<uint64_t>& y,
    Role role,
    NetPeer& self,
    NetPeer* peer_ptr = nullptr,
    NetPeer* second_peer = nullptr
) {
    // std::cout << "x.val = " << x.val << std::endl;
    // std::cout << "y.val = " << y.val << std::endl;
    uint64_t out = 0;
    if (role == Role::P2) {
        if (!peer_ptr || !second_peer) co_return 0;
        NetPeer& p0 = *peer_ptr;
        NetPeer& p1 = *second_peer;

        uint64_t X0, Y0, X1, Y1, gamma0, gamma1;
        arc4random_buf(&X0, sizeof(X0));
        arc4random_buf(&Y0, sizeof(Y0));
        arc4random_buf(&X1, sizeof(X1));
        arc4random_buf(&Y1, sizeof(Y1));
        gamma0 = X0 * Y1;
        gamma1 = X1 * Y0;

        // Send concurrently
        co_await (
            (p0 << X0) && (p1 << Y0) &&
            (p0 << X1) && (p1 << Y1) &&
            (p0 << gamma0) && (p1 << gamma1)
        );
    }
    else if (role == Role::P0) {
        if (!peer_ptr) co_return 0;
        NetPeer& peer = *peer_ptr;  // P1

        uint64_t X0, Y0, gamma0;
        co_await ((self >> X0) && (self >> Y0) && (self >> gamma0));

        uint64_t x0 = x.val;
        uint64_t y0 = y.val;
        uint64_t X_tilde, Y_tilde;
        // Send to P1 (no yield needed)
        co_await ((peer << (x0 + X0)) && (peer << (y0 + Y0)));
        co_await ((peer >> X_tilde) && (peer >> Y_tilde));

        uint64_t z0 = x0 * (y0 + Y_tilde) - (Y0 * X_tilde) + gamma0;

        // std::cout << "X_tilde = " << X_tilde << std::endl
        //           << "Y_tilde = " << Y_tilde << std::endl;
        // std::cout << "z0 = " << z0 << std::endl;

        out = z0;
    }
    else if (role == Role::P1) {
        if (!peer_ptr) co_return 0;
        NetPeer& peer = *peer_ptr;  // P0

        uint64_t X1, Y1, gamma1;
        co_await ((self >> X1) && (self >> Y1) && (self >> gamma1));

        uint64_t x1 = x.val;
        uint64_t y1 = y.val;

        uint64_t X_tilde, Y_tilde;
        co_await ((peer >> X_tilde) && (peer >> Y_tilde));
        co_await ((peer << (x1 + X1)) && (peer << (y1 + Y1)));

        uint64_t z1 = x1 * (y1 + Y_tilde) - (Y1 * X_tilde) + gamma1;

        // std::cout << "X_tilde = " << X_tilde << std::endl
        //           << "Y_tilde = " << Y_tilde << std::endl;
        // std::cout << "z1 = " << z1 << std::endl;

        out = z1;
    }

    co_return out;
}


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

    // std::cout << "x.val = " << +x.val << std::endl;
    // std::cout << "y.val = " << +y.val << std::endl;

    T out = 0;

    // ------------------------------------------------------------------------
    // P2: generates correlated randomness and sends to P0 and P1
    // ------------------------------------------------------------------------
    if (role == Role::P2) {
        if (!peer_ptr || !second_peer)
            co_return 0;

        NetPeer& p0 = *peer_ptr;
        NetPeer& p1 = *second_peer;

        T X0, Y0, X1, Y1, gamma0, gamma1;
        arc4random_buf(&X0, sizeof(T));
        arc4random_buf(&Y0, sizeof(T));
        arc4random_buf(&X1, sizeof(T));
        arc4random_buf(&Y1, sizeof(T));

        gamma0 = X0 & Y1;
        gamma1 = X1 & Y0;

    //    std::cout << "[P2] Generated correlated randomness" << std::endl;

        co_await (p0 << X0);
        co_await (p0 << Y0);
        co_await (p0 << gamma0);

        co_await (p1 << X1);
        co_await (p1 << Y1);
        co_await (p1 << gamma1);

      //  std::cout << "[P2] Sent correlated triples" << std::endl;
        co_return 0;
    }

    // ------------------------------------------------------------------------
    // P0: receives triple from P2, engages in AND protocol with P1
    // ------------------------------------------------------------------------
    else if (role == Role::P0) {
        if (!peer_ptr)
            co_return 0;

        NetPeer& peer = *peer_ptr; // P1

        T X0, Y0, gamma0;
        co_await (self >> X0);
        co_await (self >> Y0);
        co_await (self >> gamma0);

        T x0 = x.val;
        T y0 = y.val;

        // Send masked shares to P1
        T a0 = x0 ^ X0;
        T b0 = y0 ^ Y0;

        co_await (peer << a0);
        co_await (peer << b0);

        // Receive masked shares from P1
        T a1, b1;
        co_await (peer >> a1);
        co_await (peer >> b1);

        T z0 = (x0 & (y0 ^ b1)) ^ (Y0 & a1) ^ gamma0;

        // std::cout << "[P0] X_tilde=" << +a1
        //           << " Y_tilde=" << +b1
        //           << " z0=" << +z0 << std::endl;

        out = z0;
    }

    // ------------------------------------------------------------------------
    // P1: receives triple from P2, engages in AND protocol with P0
    // ------------------------------------------------------------------------
    else if (role == Role::P1) {
        if (!peer_ptr)
            co_return 0;

        NetPeer& peer = *peer_ptr; // P0

        T X1, Y1, gamma1;
        co_await (self >> X1);
        co_await (self >> Y1);
        co_await (self >> gamma1);

        T x1 = x.val;
        T y1 = y.val;

        // Receive masked shares from P0
        T a0, b0;
        co_await (peer >> a0);
        co_await (peer >> b0);

        // Send masked shares to P0
        T a1 = x1 ^ X1;
        T b1 = y1 ^ Y1;

        co_await (peer << a1);
        co_await (peer << b1);

        T z1 = (x1 & (y1 ^ b0)) ^ (Y1 & a0) ^ gamma1;

        // std::cout << "[P1] X_tilde=" << +a0
        //           << " Y_tilde=" << +b0
        //           << " z1=" << +z1 << std::endl;

        out = z1;
    }

    co_return out;
}



// --------------------- operator* using MPCContext (no globals) ---------------------
// Returns an awaitable<AShare<T>> so `auto a = co_await (x * y);` yields AShare<T>.
template<typename T>
boost::asio::awaitable<AShare<T>> operator*(const AShare<T>& x, const AShare<T>& y) {
    // Both shares should have a context pointer set (and ideally the same context).
    if (!x.ctx || !y.ctx) {
        throw std::runtime_error("AShare::operator* requires both operands to have MPCContext set");
    }
    MPCContext* ctx = x.ctx;
    if (ctx != y.ctx) {
        throw std::runtime_error("AShare::operator*: mismatched MPCContext pointers");
    }

    // Call the existing mpc_mul implementation with concrete Role and NetPeer references from ctx.
    // Note: mpc_mul currently returns uint64_t for the result value. Adapt as necessary for T != uint64_t.
    if constexpr (!std::is_same<T, uint64_t>::value) {
        // If you plan to use types other than uint64_t, implement type conversion or templated mpc_mul.
        static_assert(std::is_same<T, uint64_t>::value, "mpc_mul currently supports uint64_t only");
    }

    uint64_t result = co_await mpc_mul(
        static_cast<const AShare<uint64_t>&>(x),
        static_cast<const AShare<uint64_t>&>(y),
        ctx->role,
        ctx->self,
        ctx->peer0,
        ctx->peer1
    );

    co_return AShare<T>(static_cast<T>(result), ctx);
}


    // --------------------- mpc_eqz (unchanged logic, minor const/ref cleanup) ---------------------
    inline boost::asio::awaitable<uint64_t> mpc_eqz(
        const AShare<uint64_t>& x,
        Role role,
        NetPeer& self,
        NetPeer* peer_ptr = nullptr,
        NetPeer* second_peer = nullptr
    ) {
        uint64_t out = 0;
         uint64_t eqz = 1;
        if (role == Role::P2) {
            if (!peer_ptr || !second_peer) co_return 0;
            NetPeer& p0 = *peer_ptr;
            NetPeer& p1 = *second_peer;

            uint64_t r, r0AS, r1AS, r0XS, r1XS;
            arc4random_buf(&r, sizeof(r));
            arc4random_buf(&r0AS, sizeof(r0AS));
            arc4random_buf(&r0XS, sizeof(r0XS));
          
            r1AS = r - r0AS;
            r1XS = r ^ r0XS;
                  XShare<uint64_t> yy{1}; 
                XShare<uint64_t> xx{1};
         co_await mpc_and(xx, yy, Role::P2, p0, &p0, &p1);
            // Send concurrently
            co_await (
                (p0 << r0AS) && (p1 << r1AS) &&
                (p0 << r0XS) && (p1 << r1XS) 
            );
        }
        else if (role == Role::P0) {
            if (!peer_ptr) co_return 0;
            NetPeer& peer = *peer_ptr;  // P1

            uint64_t rAS, rXS;
            co_await ((self >> rAS) && (self >> rXS));

            uint64_t c = x.val + rAS;
            
            uint64_t c_recv;
            co_await ((peer << (c)));
            co_await ((peer >> c_recv));
            //

            uint64_t c_rec = c + c_recv;
               for (int i = 0; i < 64; i++) { 
                uint64_t c_rec_i = (c_rec >> i) & 1ULL; 
                uint64_t rXS_i = (rXS >> i) & 1ULL; 
                uint64_t diff = c_rec_i ^ rXS_i; // XOR //std::cout << "diff = " << diff << std::endl; 
                
                XShare<uint64_t> yy{diff}; 
                XShare<uint64_t> xx{eqz};
                
                auto aaa = co_await mpc_and(xx, yy, Role::P0, self, &peer);
                std::cout << "aaa = " << aaa << std::endl;

                uint64_t diff_recv;
                co_await ((peer << (1ULL - diff)));
                co_await ((peer >> (diff_recv)));
                diff_recv = diff_recv ^ (1ULL - diff);
                //std::cout << "diff_recv = " << diff_recv << std::endl;
                eqz &= (1ULL - diff_recv); 
            }

            out = eqz;
        }
        else if (role == Role::P1) {
            if (!peer_ptr) co_return 0;
            NetPeer& peer = *peer_ptr;  // P0

            uint64_t rAS, rXS;
            co_await ((self >> rAS) && (self >> rXS));
            uint64_t c = x.val + rAS;
 

            uint64_t c_recv;
            co_await ((peer << (c)));
            co_await ((peer >> c_recv));
            // 

            uint64_t c_rec = c + c_recv;
            //std::cout << "c_rec = " << c_rec << std::endl;
            uint64_t diff_recv = 0;
            for (int i = 0; i < 64; i++) { 
                //uint64_t c_rec_i = (c_rec >> i) & 1ULL; 
                uint64_t rXS_i = (rXS >> i) & 1ULL; 
                uint64_t diff =  rXS_i; // XOR //std::cout << "diff = " << diff << std::endl; 
                                
                XShare<uint64_t> yy{diff}; 
                XShare<uint64_t> xx{eqz};
                
                auto aaa = co_await mpc_and(xx, yy, Role::P0, self, &peer);
                //std::cout << "aaa = " << aaa << std::endl;
                co_await ((peer << (1ULL - diff)));
                co_await ((peer >> (diff_recv)));
                diff_recv = diff_recv ^ (1ULL - diff);
                eqz &= (1ULL - diff_recv);  

                //std::cout << "eqz = " << eqz << std::endl;
            }



            out = eqz;
        }

        co_return out;
}
