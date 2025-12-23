#ifndef LOCORAM_HPP
#define LOCORAM_HPP

#include <vector>
#include <stdexcept>
#include <cstdint>
#include <iostream>
#include <type_traits>
#include <boost/asio/awaitable.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/detached.hpp>

using boost::asio::awaitable;

#include "shares.hpp"   // defines AShare, MPCContext, Role, NetPeer, etc.

template<typename leaf_t, typename node_t = __m128i, typename prgkey_t = AES_KEY>
struct Locoram {
    std::vector<leaf_t> secret_share;
    std::vector<leaf_t> blind;
    std::vector<leaf_t> peer_blinded;
    std::vector<leaf_t> sbv_secret_share;
    std::vector<leaf_t> sbv_blind;
    std::vector<leaf_t> sbv_peer_blinded;
    std::vector<XShare<uint64_t>> update_inds;
    std::vector<XShare<uint64_t>> update_vals;
    std::vector<XShare<uint64_t>> read_inds;
    
    leaf_t cancellation_term;
    size_t nitems = 0;
    prgkey_t prgkey;

    Locoram();
    Locoram(size_t n, const prgkey_t& key);

    void set_share(const std::vector<leaf_t>& share);
    void dump(const std::string& label = "secret_share") const;
    boost::asio::awaitable<void> reconstruct_with_ctx(MPCContext* ctx,
                                                      const std::string& label = "reconstruction");

    struct Ref {
        Locoram& parent;
        size_t idx;
        Ref(Locoram& p, size_t i);

        boost::asio::awaitable<leaf_t> read(const AShare<uint64_t>& shared_index);
        boost::asio::awaitable<leaf_t> read(const XShare<uint64_t>& shared_index);
        boost::asio::awaitable<leaf_t> update(const XShare<uint64_t>& shared_index,
                                              const XShare<uint64_t>& value_share);
    };

    Ref operator[](size_t idx);
    boost::asio::awaitable<leaf_t> operator[](const AShare<uint64_t>& shared_index);
    boost::asio::awaitable<leaf_t> operator[](const XShare<uint64_t>& shared_index);

    boost::asio::awaitable<void> update(const XShare<uint64_t>& shared_index,
                                        const XShare<uint64_t>& value_share);
};

#include "locoram.tpp"

#endif // LOCORAM_HPP
