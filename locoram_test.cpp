#include <iostream>
#include <chrono>
#include <random>
#include <thread>

// Boost.Asio core
#include <boost/asio.hpp>
#include <boost/asio/steady_timer.hpp>
#include <boost/asio/awaitable.hpp>
#include <boost/asio/use_awaitable.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/detached.hpp>
#include <boost/asio/experimental/awaitable_operators.hpp>

#include "mpcops.hpp"
#include "dpfio.hpp"
#include "locoram.hpp"

using boost::asio::awaitable;
using boost::asio::use_awaitable;
using boost::asio::ip::tcp;
using namespace boost::asio;
using namespace std::chrono_literals;
using namespace boost::asio::experimental::awaitable_operators;

using leaf_t = uint64_t;    // <======== FIXED: define leaf_t globally

// -----------------------------------------------------------------------------
// Helper: create Locoram with test data
// -----------------------------------------------------------------------------
template<typename T>
Locoram<T, __m128i, AES_KEY> make_test_locoram(size_t n, const AES_KEY& key) {
    Locoram<T, __m128i, AES_KEY> LL(n, key);
    std::vector<T> vals(n);
    for (size_t i = 0; i < n; ++i)
       // vals[i] = static_cast<T>(i * 10 + 1); // predictable data
    LL.set_share(vals);
    return LL;
}

 awaitable<void> run_p2(boost::asio::io_context& io) {
    std::cout << "[P2] Listening for P0 and P1..." << std::endl;
    auto s0 = co_await make_server(io, 9000);
    auto s1 = co_await make_server(io, 9001);

    NetPeer p0(Role::P0, std::move(s0));
    NetPeer p1(Role::P1, std::move(s1));

    MPCContext ctx(Role::P2, p0, &p0, &p1);
    AES_KEY key;

    Locoram<leaf_t, __m128i, AES_KEY> loc(8, key);

    std::cout << "[P2] === READ TESTS ===" << std::endl;
    AShare<uint64_t> idx_a{3, &ctx};
    XShare<uint64_t> idx_x{2, &ctx};

    //c//o_await loc[idx_a];
    //co_await loc[idx_x];

    std::cout << "[P2] === UPDATE TEST ===" << std::endl;
    XShare<uint64_t> upd_idx{5, &ctx};
    XShare<uint64_t> upd_val{0, &ctx}; // value ignored by current update logic

    co_await loc.update(upd_idx, upd_val);
    co_await loc[idx_x];
    std::cout << "[P2] Done with Locoram tests." << std::endl;
    co_return;
}


awaitable<void> run_p0(boost::asio::io_context& io) {
    
    auto s_self = co_await connect_with_retry(io, "127.0.0.1", 9000);
    auto s_peer = co_await make_server(io, 9100);

    NetPeer self(Role::P0, std::move(s_self));
    NetPeer peer(Role::P1, std::move(s_peer));
    MPCContext ctx(Role::P0, self, &peer);

    AES_KEY key;
    auto loc = make_test_locoram<leaf_t>(8, key);

    std::cout << "[P0] === BEFORE UPDATE ===" << std::endl;

    AShare<uint64_t> idx_a{0, &ctx};
    XShare<uint64_t> idx_x{0, &ctx};
    //co_await loc[idx_a];

    std::cout << "[P0] === UPDATE ===" << std::endl;
    XShare<uint64_t> upd_idx{5, &ctx};
    XShare<uint64_t> upd_val{0, &ctx};
    co_await loc.update(upd_idx, upd_val);

    std::cout << "[P0] === AFTER UPDATE (RECONSTRUCT) ===" << std::endl;
    co_await loc.reconstruct_with_ctx(&ctx, "P0 post-update");


    auto dot = co_await loc[idx_x];
    std::cout << "dot = " << dot << std::endl; 
    XShare<uint64_t> my_share{dot, &ctx};
    auto reconstruct_add = co_await reconstruct_remote(peer, my_share);
    std::cout << "reconstruct_add = " << reconstruct_add << std::endl;
    co_return;
}


awaitable<void> run_p1(boost::asio::io_context& io) {
    auto s_self = co_await connect_with_retry(io, "127.0.0.1", 9001);
    auto s_peer = co_await connect_with_retry(io, "127.0.0.1", 9100);

    NetPeer self(Role::P1, std::move(s_self));
    NetPeer peer(Role::P0, std::move(s_peer));
    MPCContext ctx(Role::P1, self, &peer);

    AES_KEY key;
    auto loc = make_test_locoram<leaf_t>(8, key);

    std::cout << "[P1] === BEFORE UPDATE ===" << std::endl;
    AShare<uint64_t> idx_a{5, &ctx};
    XShare<uint64_t> idx_x{3, &ctx};
    //co_await loc[idx_a];

    std::cout << "[P1] === UPDATE ===" << std::endl;
    XShare<uint64_t> upd_idx{0, &ctx};
    XShare<uint64_t> upd_val{100, &ctx};
    co_await loc.update(upd_idx, upd_val);

    std::cout << "[P1] === AFTER UPDATE (RECONSTRUCT) ===" << std::endl;
    co_await loc.reconstruct_with_ctx(&ctx, "P1 post-update");


    std::cout << "[P1] === BEFORE UPDATE ===" << std::endl << std::endl 
                                            << std::endl << std::endl << std::endl; 
 auto dot = co_await loc[idx_x];
 std::cout << "dot = " << dot << std::endl << std::endl << std::endl << std::endl; 
 XShare<uint64_t> my_share{dot, &ctx};
auto reconstruct_add = co_await reconstruct_remote(peer, my_share);
    std::cout << "reconstruct_add = " << reconstruct_add << std::endl;
    co_return;
}

// -----------------------------------------------------------------------------
// main()
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "usage: ./locoram_test [p0|p1|p2]\n";
        return 1;
    }

    boost::asio::io_context io;
    std::string role = argv[1];

    if (role == "p0")
        boost::asio::co_spawn(io, run_p0(io), boost::asio::detached);
    else if (role == "p1")
        boost::asio::co_spawn(io, run_p1(io), boost::asio::detached);
    else if (role == "p2")
        boost::asio::co_spawn(io, run_p2(io), boost::asio::detached);
    else {
        std::cerr << "Invalid role.\n";
        return 1;
    }

    io.run();
}
