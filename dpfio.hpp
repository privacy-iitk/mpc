#pragma once
 
#include "dpf.hpp"
#include <boost/asio.hpp>
#include <vector>

using boost::asio::awaitable;
using boost::asio::use_awaitable;
using boost::asio::async_write;
using boost::asio::async_read;
using boost::asio::buffer;

namespace dpf_io {

// =====================================================
// Send DPF key
// =====================================================
template<typename leaf_t, typename node_t, typename prgkey_t>
awaitable<void> send_dpf_key(NetPeer& p, const dpf::dpf_key<leaf_t, node_t, prgkey_t>& key, bool debug = false) {
    // Serialize key
    auto buf = dpf::serialize_dpf_key<leaf_t, node_t, prgkey_t>(key);
    uint64_t len = static_cast<uint64_t>(buf.size());

    // Send length and buffer
    co_await p.send<uint64_t>(len);
    if (len > 0) {
        co_await async_write(p.sock, buffer(buf.data(), len), use_awaitable);
    }

    if (debug) std::cout << "[send_dpf_key] Sent " << len << " bytes\n";
    co_return;
}

// =====================================================
// Receive DPF key
// =====================================================
template<typename leaf_t, typename node_t, typename prgkey_t>
awaitable<dpf::dpf_key<leaf_t, node_t, prgkey_t>> recv_dpf_key(NetPeer& p, bool debug = false) {
    // Receive length
    uint64_t len = co_await p.recv<uint64_t>();
    std::vector<uint8_t> buf(static_cast<size_t>(len));

    // Receive buffer
    if (len > 0) {
        co_await async_read(p.sock, buffer(buf.data(), len), use_awaitable);
    }

    if (debug) std::cout << "[recv_dpf_key] Received " << len << " bytes\n";

    // Deserialize and return
    co_return dpf::deserialize_dpf_key<leaf_t, node_t, prgkey_t>(buf.data(), buf.size());
}

    template<typename leaf_t, typename node_t, typename prgkey_t>
    awaitable<void> operator<<(NetPeer& p, const dpf::dpf_key<leaf_t, node_t, prgkey_t>& key) {
        co_await send_dpf_key<leaf_t, node_t, prgkey_t>(p, key);
    }

    template<typename leaf_t, typename node_t, typename prgkey_t>
    awaitable<void> operator>>(NetPeer& p, dpf::dpf_key<leaf_t, node_t, prgkey_t>& key) {
        key = co_await recv_dpf_key<leaf_t, node_t, prgkey_t>(p);
    }

} // namespace dpf_io




 