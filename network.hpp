#pragma once

#include <boost/asio.hpp>
#include <boost/asio/awaitable.hpp>
#include <boost/asio/use_awaitable.hpp>
#include <boost/asio/experimental/parallel_group.hpp>

#include <vector>
#include <cstdint>
#include <string>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <type_traits>

using boost::asio::awaitable;
using boost::asio::use_awaitable;
using boost::asio::async_write;
using boost::asio::async_read;
using boost::asio::buffer;
using boost::asio::ip::tcp;

// Role enum
enum class Role { P0, P1, P2 };

// NetPeer: lightweight connected peer wrapper
struct NetPeer {
    Role role;
    tcp::socket sock;

    explicit NetPeer(Role r, tcp::socket&& s) : role(r), sock(std::move(s)) {}

    // primitive send/recv (trivially copyable types)
    template<typename T>
    awaitable<void> send(const T& value) {
        static_assert(std::is_trivially_copyable_v<T>);
        co_await async_write(sock, buffer(&value, sizeof(T)), use_awaitable);
        co_return;
    }

    template<typename T>
    awaitable<T> recv() {
        static_assert(std::is_trivially_copyable_v<T>);
        T value;
        co_await async_read(sock, buffer(&value, sizeof(T)), use_awaitable);
        co_return value;
    }

    // vector send/recv (prefix length + raw data)
    template<typename T>
    awaitable<void> send(const std::vector<T>& vec) {
        static_assert(std::is_trivially_copyable_v<T>);
        uint64_t sz = vec.size();
        co_await send(sz);
        if (sz > 0) co_await async_write(sock, buffer(vec.data(), sizeof(T) * sz), use_awaitable);
        co_return;
    }

    template<typename T>
    awaitable<std::vector<T>> recv_vector() {
        static_assert(std::is_trivially_copyable_v<T>);
        uint64_t sz = co_await recv<uint64_t>();
        std::vector<T> v(sz);
        if (sz > 0) co_await async_read(sock, buffer(v.data(), sizeof(T) * sz), use_awaitable);
        co_return v;
    }

    // legacy operator<< / >> that return awaitable<void> to preserve old syntax
    template<typename T>
    awaitable<void> operator<<(const T& val) {
        co_await send(val);
        co_return;
    }

    template<typename T>
    awaitable<void> operator>>(T& out) {
        out = co_await recv<T>();
        co_return;
    }

    template<typename T>
    awaitable<void> operator<<(const std::vector<T>& vec) {
        co_await send(vec);
        co_return;
    }

    template<typename T>
    awaitable<void> operator>>(std::vector<T>& out) {
        out = co_await recv_vector<T>();
        co_return;
    }
};

awaitable<tcp::socket> connect_with_retry(boost::asio::io_context& io,
                                          const std::string& host,
                                          uint16_t port,
                                          int max_tries = -1) {
    tcp::resolver resolver(io);
    tcp::socket sock(io);
    int attempt = 0;

    for (;;) {
        try {
            auto endpoints = co_await resolver.async_resolve(host, std::to_string(port), use_awaitable);
            co_await boost::asio::async_connect(sock, endpoints, use_awaitable);
            co_return std::move(sock);
        } catch (std::exception& e) {
            attempt++;
            if (max_tries > 0 && attempt >= max_tries)
                throw;
            std::cerr << "[connect_with_retry] " << host << ":" << port
                      << " failed (" << e.what() << "), retrying..." << std::endl;
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
        }
    }
}

// helpers to create sockets
inline awaitable<tcp::socket> make_server(boost::asio::io_context& io, uint16_t port) {
    tcp::acceptor acceptor(io, tcp::endpoint(tcp::v4(), port));
    tcp::socket sock(io);
    co_await acceptor.async_accept(sock, use_awaitable);
    co_return std::move(sock);
}

inline awaitable<tcp::socket> make_client(boost::asio::io_context& io, const std::string& host, uint16_t port) {
    tcp::resolver resolver(io);
    auto results = co_await resolver.async_resolve(host, std::to_string(port), use_awaitable);
    tcp::socket sock(io);
    co_await boost::asio::async_connect(sock, results, use_awaitable);
    co_return std::move(sock);
}

// NetContext: container for multiple peers
struct NetContext {
    Role self_role;
    std::vector<std::unique_ptr<NetPeer>> peers;
    explicit NetContext(Role r) : self_role(r) {}

    void add_peer(Role r, tcp::socket&& sock) {
        peers.emplace_back(std::make_unique<NetPeer>(r, std::move(sock)));
    }

    NetPeer& peer(Role r) {
        for (auto& p : peers) if (p->role == r) return *p;
        throw std::runtime_error("peer not found");
    }
};
