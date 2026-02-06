 
 
#pragma once

#include <boost/asio.hpp>
#include <boost/asio/awaitable.hpp>
#include <boost/asio/use_awaitable.hpp>
#include <boost/asio/experimental/awaitable_operators.hpp>

#include <vector>
#include <cstdint>
#include <string>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <type_traits>
#include <tuple>
#include <utility>
#include <thread>

using boost::asio::awaitable;
using boost::asio::use_awaitable;
using boost::asio::async_write;
using boost::asio::async_read;
using boost::asio::buffer;
using boost::asio::ip::tcp;
using namespace boost::asio::experimental::awaitable_operators;

// -----------------------------------------------------------------------------
// Role enum
// -----------------------------------------------------------------------------
enum class Role { P0, P1, P2 };

// -----------------------------------------------------------------------------
// Bundle
// -----------------------------------------------------------------------------
template<typename... Ts>
struct Bundle {
    std::tuple<Ts...> data;
    Bundle() = default;
    Bundle(Ts... xs) : data(std::move(xs)...) {}
};

// -----------------------------------------------------------------------------
// Type trait: detect std::vector
// -----------------------------------------------------------------------------
template<typename T>
struct is_std_vector : std::false_type {};

template<typename T, typename A>
struct is_std_vector<std::vector<T, A>> : std::true_type {};

template<typename T>
inline constexpr bool is_std_vector_v = is_std_vector<T>::value;

// -----------------------------------------------------------------------------
// Coroutine-safe tuple iteration
// -----------------------------------------------------------------------------
namespace net_detail {

template<std::size_t I = 0, typename Tuple, typename F>
awaitable<void> tuple_for_each(Tuple& t, F&& f) {
    if constexpr (I < std::tuple_size_v<Tuple>) {
        co_await f(std::get<I>(t));
        co_await tuple_for_each<I + 1>(t, std::forward<F>(f));
    } else {
        co_return;
    }
}

} // namespace net_detail

// -----------------------------------------------------------------------------
// NetPeer
// -----------------------------------------------------------------------------
struct NetPeer {
    Role role;
    tcp::socket sock;

    // explicit NetPeer(Role r, tcp::socket&& s)
    //     : role(r), sock(std::move(s)) {}

explicit NetPeer(Role r, tcp::socket&& s)
    : role(r), sock(std::move(s))
{
    sock.set_option(tcp::no_delay(true));
}


    // -----------------------------------------------------------------
    // Primitive send / recv (UNCHANGED API)
    // -----------------------------------------------------------------
    template<typename T>
    awaitable<void> send(const T& value) {
        static_assert(std::is_trivially_copyable_v<T>);
        co_await async_write(sock, buffer(&value, sizeof(T)), use_awaitable);
        co_return;
    }

    template<typename T>
    awaitable<T> recv() {
        static_assert(std::is_trivially_copyable_v<T>);
        T value{};
        co_await async_read(sock, buffer(&value, sizeof(T)), use_awaitable);
        co_return value;
    }

    // -----------------------------------------------------------------
    // Vector send / recv
    // -----------------------------------------------------------------
    template<typename T>
    awaitable<void> send(const std::vector<T>& vec) {
        static_assert(std::is_trivially_copyable_v<T>);
        uint64_t sz = vec.size();
        co_await send(sz);
        if (sz > 0)
            co_await async_write(
                sock, buffer(vec.data(), sizeof(T) * sz), use_awaitable);
        co_return;
    }

    template<typename T>
    awaitable<std::vector<T>> recv_vector() {
        static_assert(std::is_trivially_copyable_v<T>);
        uint64_t sz = co_await recv<uint64_t>();
        std::vector<T> v(sz);
        if (sz > 0)
            co_await async_read(
                sock, buffer(v.data(), sizeof(T) * sz), use_awaitable);
        co_return v;
    }

    // -----------------------------------------------------------------
    // Bundle send
    // -----------------------------------------------------------------
    template<typename... Ts>
    awaitable<void> send_bundle(const Bundle<Ts...>& b) {
        co_await net_detail::tuple_for_each(
            b.data,
            [&](auto const& x) -> awaitable<void> {
                co_await send(x);
                co_return;
            });
        co_return;
    }

    // -----------------------------------------------------------------
    // Bundle recv (FIXED, backward-compatible)
    // -----------------------------------------------------------------
    template<typename... Ts>
    awaitable<Bundle<Ts...>> recv_bundle() {
        Bundle<Ts...> b;
        co_await net_detail::tuple_for_each(
            b.data,
            [&](auto& x) -> awaitable<void> {
                using T = std::decay_t<decltype(x)>;
                if constexpr (is_std_vector_v<T>) {
                    using Elem = typename T::value_type;
                    x = co_await recv_vector<Elem>();
                } else {
                    x = co_await recv<T>();
                }
                co_return;
            });
        co_return b;
    }

    // -----------------------------------------------------------------
    // Operator overloads (UNCHANGED)
    // -----------------------------------------------------------------
    template<typename T>
    awaitable<void> operator<<(const T& v) {
        co_await send(v);
        co_return;
    }

    template<typename T>
    awaitable<void> operator>>(T& v) {
        v = co_await recv<T>();
        co_return;
    }

    template<typename T>
    awaitable<void> operator<<(const std::vector<T>& v) {
        co_await send(v);
        co_return;
    }

    template<typename T>
    awaitable<void> operator>>(std::vector<T>& v) {
        v = co_await recv_vector<T>();
        co_return;
    }

    template<typename... Ts>
    awaitable<void> operator<<(const Bundle<Ts...>& b) {
        co_await send_bundle(b);
        co_return;
    }

    template<typename... Ts>
    awaitable<void> operator>>(Bundle<Ts...>& b) {
        b = co_await recv_bundle<Ts...>();
        co_return;
    }

    
};

inline boost::asio::awaitable<tcp::socket>
connect_with_retry(boost::asio::io_context& io,
                   const std::string& host,
                   uint16_t port,
                   int max_tries = -1)
{ 
    tcp::resolver resolver(io);
    tcp::socket sock(io);
    int attempt = 0; 

    for (;;) {
        try {
            auto endpoints =
                co_await resolver.async_resolve(
                    host, std::to_string(port), use_awaitable);

            co_await boost::asio::async_connect(
                sock, endpoints, use_awaitable);

            co_return std::move(sock);
        } catch (const std::exception& e) {
            ++attempt;
            if (max_tries > 0 && attempt >= max_tries)
                throw;

            std::cerr << "[connect_with_retry] "
                      << host << ":" << port
                      << " failed (" << e.what()
                      << "), retrying...\n";

            std::this_thread::sleep_for(
                std::chrono::milliseconds(500));
        }
    }
}


// -----------------------------------------------------------------------------
// Connection helpers (RESTORED)
// -----------------------------------------------------------------------------
inline awaitable<tcp::socket>
make_server(boost::asio::io_context& io, uint16_t port) {
    tcp::acceptor acc(io, tcp::endpoint(tcp::v4(), port));
    tcp::socket sock(io);
    co_await acc.async_accept(sock, use_awaitable);
    co_return std::move(sock);
}

inline awaitable<tcp::socket>
make_client(boost::asio::io_context& io,
            const std::string& host,
            uint16_t port) {
    tcp::resolver resolver(io);
    auto eps = co_await resolver.async_resolve(
        host, std::to_string(port), use_awaitable);
    tcp::socket sock(io);
    co_await boost::asio::async_connect(sock, eps, use_awaitable);
    co_return std::move(sock);
}

// -----------------------------------------------------------------------------
// NetContext
// -----------------------------------------------------------------------------
struct NetContext {
    Role self_role;
    std::vector<std::unique_ptr<NetPeer>> peers;

    explicit NetContext(Role r) : self_role(r) {}

    void add_peer(Role r, tcp::socket&& sock) {
        peers.emplace_back(
            std::make_unique<NetPeer>(r, std::move(sock)));
    }

    NetPeer& peer(Role r) {
        for (auto& p : peers)
            if (p->role == r) return *p;
        throw std::runtime_error("peer not found");
    }
};
