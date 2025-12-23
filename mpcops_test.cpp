// #include <boost/asio.hpp>
// #include <boost/asio/co_spawn.hpp>
// #include <boost/asio/detached.hpp>

// #include "mpcops.hpp"

// using boost::asio::awaitable;
// using boost::asio::use_awaitable;
// using boost::asio::ip::tcp;
// using namespace boost::asio::experimental::awaitable_operators;

// // -------------------------------------------------------------
// // P2 – Generates correlated randomness
// // -------------------------------------------------------------
// awaitable<void> run_p2(boost::asio::io_context& io) {
//     std::cout << "[P2] Listening for P0 and P1..." << std::endl;

//     auto s0 = co_await make_server(io, 9000);
//     auto s1 = co_await make_server(io, 9001);

//     NetPeer p0(Role::P0, std::move(s0));
//     NetPeer p1(Role::P1, std::move(s1));

//     MPCContext ctx(Role::P2, p0, &p0, &p1);

//     // ---------------------------------------------------------
//     // Existing multiplication tests
//     // ---------------------------------------------------------
//     AShare<uint64_t> a{0, &ctx};
//     AShare<uint64_t> b{0, &ctx};
//     co_await (a * b);

//     XShare<uint64_t> x{1, &ctx};
//     XShare<uint64_t> y{1, &ctx};

//     co_await (x | y);
//     co_await (x * y);

//     // ---------------------------------------------------------
//     // Dot-product test
//     // ---------------------------------------------------------
//     std::cout << "[P2] Running MPC dot product..." << std::endl;

//     std::vector<uint64_t> vv1 = {10, 20, 30};
//     std::vector<uint64_t> vv2 = {1, 2, 3};

//     AdditiveShareVector<uint64_t> aa(vv1, &ctx);
//     AdditiveShareVector<uint64_t> bb(vv2, &ctx);

//     co_await (aa * bb);
//     std::cout << "[P2] Done with dot product." << std::endl;

//     // ---------------------------------------------------------
//     // EQZ tests
//     // ---------------------------------------------------------
//     std::cout << "[P2] Running MPC EQZ tests..." << std::endl;

//     AShare<uint64_t> x0{0, &ctx};
//     AShare<uint64_t> x1{0, &ctx};

//     co_await mpc_eqz(x0, Role::P2, p0, &p0, &p1);
//     co_await mpc_eqz(x1, Role::P2, p0, &p0, &p1);

//     std::cout << "[P2] EQZ tests dispatched." << std::endl;
//     co_return;
// }

// // -------------------------------------------------------------
// // P0 – connects to P2, waits for P1
// // -------------------------------------------------------------
// awaitable<void> run_p0(boost::asio::io_context& io) {
//     std::cout << "[P0] Connecting to P2..." << std::endl;

//     auto s_self = co_await connect_with_retry(io, "127.0.0.1", 9000);

//     std::cout << "[P0] Waiting for P1 on port 9100..." << std::endl;
//     auto s_peer = co_await make_server(io, 9100);

//     NetPeer self(Role::P0, std::move(s_self));
//     NetPeer peer(Role::P1, std::move(s_peer));
//     MPCContext ctx(Role::P0, self, &peer);

//     // ---------------------------------------------------------
//     // Scalar tests
//     // ---------------------------------------------------------
//     AShare<uint64_t> a{5, &ctx};
//     AShare<uint64_t> b{7, &ctx};

//     AShare<uint64_t> c = co_await (a * b);
//     std::cout << "[P0] a*b = " << c.val << std::endl;

//     auto reconstruct_add = co_await reconstruct_remote(peer, c);
//     std::cout << "[P0] Reconstructed multiplication = "
//               << reconstruct_add << std::endl;

//     XShare<uint64_t> x{1, &ctx};
//     XShare<uint64_t> y{0, &ctx};

//     auto x_or_y = co_await (x | y);
//     auto reconstruct_or = co_await reconstruct_remote(peer, x_or_y);
//     std::cout << "[P0] OR = " << reconstruct_or << std::endl;

//     auto s4 = co_await (x * y);
//     auto reconstruct_and = co_await reconstruct_remote(peer, s4);
//     std::cout << "[P0] AND = " << reconstruct_and << std::endl;

//     // ---------------------------------------------------------
//     // Dot-product test
//     // ---------------------------------------------------------
//     std::cout << "[P0] Running MPC dot product..." << std::endl;

//     std::vector<uint64_t> vv1 = {10, 20, 30};
//     std::vector<uint64_t> vv2 = {1, 2, 3};

//     AdditiveShareVector<uint64_t> aa(vv1, &ctx);
//     AdditiveShareVector<uint64_t> bb(vv2, &ctx);

//     AShare<uint64_t> aaa = co_await (aa * bb);
//     auto reconstruct_dot = co_await reconstruct_remote(peer, aaa);
//     std::cout << "[P0] Dot product = " << reconstruct_dot << std::endl;

//     // ---------------------------------------------------------
//     // EQZ tests
//     // ---------------------------------------------------------
//     std::cout << "[P0] Running MPC EQZ tests..." << std::endl;

//     AShare<uint64_t> x_zero{0, &ctx};
//     AShare<uint64_t> x_nonzero{375, &ctx};

//     XShare<uint64_t> z0 =
//         co_await (x_zero == 0);//      mpc_eqz(x_zero, Role::P0, self, &peer, nullptr);

//     XShare<uint64_t> z1 =
//         co_await  (x_nonzero == 0);   // mpc_eqz(x_nonzero, Role::P0, self, &peer, nullptr);

//     auto r0 = co_await reconstruct_remote(peer, z0);
//     std::cout << "EQZ(0) = " << r0 << std::endl;

//     auto r1 = co_await reconstruct_remote(peer, z1);
//     std::cout << "EQZ(375) = " << r1 << std::endl;

//     co_return;
// }

// // -------------------------------------------------------------
// // P1 – connects to P2 and to P0
// // -------------------------------------------------------------
// awaitable<void> run_p1(boost::asio::io_context& io) {
//     std::cout << "[P1] Connecting to P2 and P0..." << std::endl;

//     auto s_self = co_await connect_with_retry(io, "127.0.0.1", 9001);
//     auto s_peer = co_await connect_with_retry(io, "127.0.0.1", 9100);

//     NetPeer self(Role::P1, std::move(s_self));
//     NetPeer peer(Role::P0, std::move(s_peer));
//     MPCContext ctx(Role::P1, self, &peer);

//     // ---------------------------------------------------------
//     // Scalar tests
//     // ---------------------------------------------------------
//     AShare<uint64_t> a{5, &ctx};
//     AShare<uint64_t> b{7, &ctx};

//     AShare<uint64_t> c = co_await (a * b);
//     auto reconstruct_mul = co_await reconstruct_remote(peer, c);
//     std::cout << "reconstruct_mul = " << reconstruct_mul << std::endl;

//     XShare<uint64_t> x{1, &ctx};
//     XShare<uint64_t> y{1, &ctx};

//     auto x_or_y = co_await (x | y);
//     auto reconstruct_or = co_await reconstruct_remote(peer, x_or_y);
//     std::cout << "[P1] OR = " << reconstruct_or << std::endl;

//     auto s4 = co_await (x * y);
//     auto reconstruct_and = co_await reconstruct_remote(peer, s4);
//     std::cout << "[P1] AND = " << reconstruct_and << std::endl;

//     // ---------------------------------------------------------
//     // Dot-product test
//     // ---------------------------------------------------------
//     std::cout << "[P1] Running MPC dot product..." << std::endl;

//     std::vector<uint64_t> vv1 = {10, 20, 30};
//     std::vector<uint64_t> vv2 = {1, 2, 3};

//     AdditiveShareVector<uint64_t> aa(vv1, &ctx);
//     AdditiveShareVector<uint64_t> bb(vv2, &ctx);

//     AShare<uint64_t> aaa = co_await (aa * bb);
//     auto reconstruct_dot = co_await reconstruct_remote(peer, aaa);
//     std::cout << "[P1] Dot product = " << reconstruct_dot << std::endl;

//     // ---------------------------------------------------------
//     // EQZ tests
//     // ---------------------------------------------------------
//     AShare<uint64_t> x_zero{0, &ctx};
//     AShare<uint64_t> x_nonzero{1330, &ctx};

//     XShare<uint64_t> z0 =
//         co_await (x_zero == 0);//      mpc_eqz(x_zero, Role::P0, self, &peer, nullptr);

//     XShare<uint64_t> z1 =
//         co_await  (x_nonzero == 0);   // mpc_eqz(x_nonzero, Role::P0, self, &peer, nullptr);

//     auto r0 = co_await reconstruct_remote(peer, z0);
//     std::cout << "EQZ(0) = " << r0 << std::endl;

//     auto r1 = co_await reconstruct_remote(peer, z1);
//     std::cout << "EQZ(10) = " << r1 << std::endl;

//     co_return;
// }

// // -------------------------------------------------------------
// // main
// // -------------------------------------------------------------
// int main(int argc, char* argv[]) {
//     if (argc != 2) {
//         std::cerr << "usage: ./mpcops_test [p0|p1|p2]\n";
//         return 1;
//     }

//     boost::asio::io_context io;
//     std::string r = argv[1];

//     if (r == "p0") {
//         boost::asio::co_spawn(io, run_p0(io), boost::asio::detached);
//     } else if (r == "p1") {
//         boost::asio::co_spawn(io, run_p1(io), boost::asio::detached);
//     } else if (r == "p2") {
//         boost::asio::co_spawn(io, run_p2(io), boost::asio::detached);
//     } else {
//         return 1;
//     }

//     io.run();
// }


#include <boost/asio.hpp>
#include <boost/asio/co_spawn.hpp>
#include <boost/asio/detached.hpp>

#include "mpcops.hpp"

using boost::asio::awaitable;
using boost::asio::use_awaitable;
using boost::asio::ip::tcp;
using namespace boost::asio::experimental::awaitable_operators;

// -------------------------------------------------------------
// Helper: init randomness mode + files
// -------------------------------------------------------------
static void init_randomness(MPCContext& ctx, bool preprocess, bool replay) {
    std::string fname = "p" + std::to_string(int(ctx.role)) + ".rand";

    if (preprocess) {
        ctx.rand_mode = RandomnessMode::Record;
        ctx.rand_out.open(fname, std::ios::binary);
        if (!ctx.rand_out)
            throw std::runtime_error("Failed to open " + fname + " for record");
        
        ctx.mul_out_p0.open("mul_p0.rand", std::ios::binary);
        ctx.mul_out_p1.open("mul_p1.rand", std::ios::binary);
        ctx.and_out_p0.open("and_p0.rand", std::ios::binary);
        ctx.and_out_p1.open("and_p1.rand", std::ios::binary);
    }
    else if (replay) {
        ctx.rand_mode = RandomnessMode::Replay;
        ctx.rand_in.open(fname, std::ios::binary);
        if (!ctx.rand_in)
            throw std::runtime_error("Failed to open " + fname + " for replay");

        ctx.mul_in_p0.open("mul_p0.rand", std::ios::binary);
        ctx.mul_in_p1.open("mul_p1.rand", std::ios::binary);
        ctx.and_in_p0.open("and_p0.rand", std::ios::binary);
        ctx.and_in_p1.open("and_p1.rand", std::ios::binary);
    }
    else {
        ctx.rand_mode = RandomnessMode::Online;
    }

    std::cout << "[P" << int(ctx.role) << "] RandomnessMode = ";
    switch (ctx.rand_mode) {
        case RandomnessMode::Record:
            std::cout << "Record";
            break;
        case RandomnessMode::Replay:
            std::cout << "Replay";
            break;
        case RandomnessMode::Online:
            std::cout << "Online";
            break;
    }
    std::cout << std::endl;
}

// -------------------------------------------------------------
// P2
// -------------------------------------------------------------
awaitable<void> run_p2(boost::asio::io_context& io, bool preprocess, bool replay) {
    std::cout << "[P2] Listening for P0 and P1..." << std::endl;

    auto s0 = co_await make_server(io, 9000);
    auto s1 = co_await make_server(io, 9001);

    NetPeer p0(Role::P0, std::move(s0));
    NetPeer p1(Role::P1, std::move(s1));

    MPCContext ctx(Role::P2, p0, &p0, &p1);
    init_randomness(ctx, preprocess, replay);

    AShare<uint64_t> a{0, &ctx};
    AShare<uint64_t> b{0, &ctx};
    co_await (a * b);

    AShare<uint64_t> a1{0, &ctx};
    AShare<uint64_t> b1{0, &ctx};
    co_await (a1 * b1);

    XShare<uint64_t> x{1, &ctx};
    XShare<uint64_t> y{1, &ctx};
    co_await (x | y);
    co_await (x * y);

    std::vector<uint64_t> vv1 = {10, 20, 30};
    std::vector<uint64_t> vv2 = {1, 2, 3};
    AdditiveShareVector<uint64_t> aa(vv1, &ctx);
    AdditiveShareVector<uint64_t> bb(vv2, &ctx);
    co_await (aa * bb);

    AShare<uint64_t> x0{0, &ctx};
    AShare<uint64_t> x1{0, &ctx};
    co_await mpc_eqz(x0, Role::P2, p0, &p0, &p1);
    co_await mpc_eqz(x1, Role::P2, p0, &p0, &p1);

    if (ctx.rand_mode == RandomnessMode::Replay) {
        char extra;
        if (ctx.rand_in.read(&extra, 1))
            throw std::runtime_error("[P2] Replay mismatch: extra randomness");
    }

    // std::cout << "[P2] Random values consumed = "
    //           << ctx.rand_counter << std::endl;
    co_return;
}

// -------------------------------------------------------------
// P0
// -------------------------------------------------------------
awaitable<void> run_p0(boost::asio::io_context& io, bool preprocess, bool replay) {
    auto s_self = co_await connect_with_retry(io, "127.0.0.1", 9000);
    auto s_peer = co_await make_server(io, 9100);

    NetPeer self(Role::P0, std::move(s_self));
    NetPeer peer(Role::P1, std::move(s_peer));
    MPCContext ctx(Role::P0, self, &peer);
    init_randomness(ctx, preprocess, replay);

    AShare<uint64_t> a{5, &ctx};
    AShare<uint64_t> b{7, &ctx};
    AShare<uint64_t> c = co_await (a * b);
    auto rec = co_await reconstruct_remote(peer, c);
    std::cout << "[P0] mul = " << rec << std::endl;

    AShare<uint64_t> a1{15, &ctx};
    AShare<uint64_t> b1{17, &ctx};
    AShare<uint64_t> c1 = co_await (a1 * b1);
    auto rec1 = co_await reconstruct_remote(peer, c1);
    std::cout << "[P0] mul = " << rec1 << std::endl;


    XShare<uint64_t> x{1, &ctx};
    XShare<uint64_t> y{0, &ctx};
    auto r_or = co_await reconstruct_remote(peer, co_await (x | y));
    std::cout << "r_or = " << r_or << std::endl;

    auto r_and = co_await reconstruct_remote(peer, co_await (x * y));
    std::cout << "r_and = " << r_and << std::endl;
    
    std::vector<uint64_t> vv1 = {10, 20, 30};
    std::vector<uint64_t> vv2 = {1, 2, 3};
    AdditiveShareVector<uint64_t> aa(vv1, &ctx);
    AdditiveShareVector<uint64_t> bb(vv2, &ctx);
 
    AShare<uint64_t> dot = co_await (aa * bb);
    auto r_dot = co_await reconstruct_remote(peer, dot);
    std::cout << "[P0] Dot product = " << r_dot << std::endl;



    XShare<uint64_t> z0 = co_await (AShare<uint64_t>{0, &ctx} == 0);
    XShare<uint64_t> z1 = co_await (AShare<uint64_t>{375, &ctx} == 0);

    auto rz0 = co_await reconstruct_remote(peer, z0);
    auto rz1 = co_await reconstruct_remote(peer, z1);

    std::cout << "rz0 = " << rz0 << std::endl;
    std::cout << "rz1 = " << rz1 << std::endl;

    if (ctx.rand_mode == RandomnessMode::Replay) {
        char extra;
        if (ctx.rand_in.read(&extra, 1))
            throw std::runtime_error("[P0] Replay mismatch: extra randomness");
    }

    // std::cout << "[P0] Random values consumed = "
    //           << ctx.rand_counter << std::endl;
    co_return;
}

// -------------------------------------------------------------
// P1
// -------------------------------------------------------------
awaitable<void> run_p1(boost::asio::io_context& io, bool preprocess, bool replay) {
    auto s_self = co_await connect_with_retry(io, "127.0.0.1", 9001);
    auto s_peer = co_await connect_with_retry(io, "127.0.0.1", 9100);

    NetPeer self(Role::P1, std::move(s_self));
    NetPeer peer(Role::P0, std::move(s_peer));
    MPCContext ctx(Role::P1, self, &peer);
    init_randomness(ctx, preprocess, replay);

    AShare<uint64_t> a{5, &ctx};
    AShare<uint64_t> b{7, &ctx};
    auto r_mul = co_await reconstruct_remote(peer, co_await (a * b));
    std::cout << "[P0] mul = " << r_mul << std::endl;

    AShare<uint64_t> a1{15, &ctx};
    AShare<uint64_t> b1{17, &ctx};
    AShare<uint64_t> c1 = co_await (a1 * b1);
    auto rec1 = co_await reconstruct_remote(peer, c1);
    std::cout << "[P0] mul = " << rec1 << std::endl;

    XShare<uint64_t> x{1, &ctx};
    XShare<uint64_t> y{1, &ctx};
    
    auto r_or  = co_await reconstruct_remote(peer, co_await (x | y));
    std::cout << "r_or = " << r_or << std::endl;

    auto r_and = co_await reconstruct_remote(peer, co_await (x * y));
    std::cout << "r_and = " << r_and << std::endl;

    std::vector<uint64_t> vv1 = {10, 20, 30};
    std::vector<uint64_t> vv2 = {1, 2, 3};
    AdditiveShareVector<uint64_t> aa(vv1, &ctx);
    AdditiveShareVector<uint64_t> bb(vv2, &ctx);

    AShare<uint64_t> dot = co_await (aa * bb);  
    auto r_dot = co_await reconstruct_remote(peer, dot);
    std::cout << "[P0] Dot product = " << r_dot << std::endl;

    XShare<uint64_t> z0 = co_await (AShare<uint64_t>{0, &ctx} == 0);
    XShare<uint64_t> z1 = co_await (AShare<uint64_t>{375, &ctx} == 0);

    auto rz0 = co_await reconstruct_remote(peer, z0);
    auto rz1 = co_await reconstruct_remote(peer, z1);

    std::cout << "rz0 = " << rz0 << std::endl;
    std::cout << "rz1 = " << rz1 << std::endl;

    if (ctx.rand_mode == RandomnessMode::Replay) {
        char extra;
        if (ctx.rand_in.read(&extra, 1))
            throw std::runtime_error("[P1] Replay mismatch: extra randomness");
    }

    // std::cout << "[P1] Random values consumed = "
    //           << ctx.rand_counter << std::endl;
    co_return;
}

// -------------------------------------------------------------
// main
// -------------------------------------------------------------
int main(int argc, char* argv[]) {
    bool preprocess = false;
    bool replay = false;
    std::string role;

    if (argc == 2) {

        std::cout << "No Preprocessing Here" << std::endl;
        role = argv[1];
    }
    else if (argc == 3 && std::string(argv[1]) == "-p") {
        preprocess = true;
        role = argv[2];
    }
    
    else if (argc == 3 && std::string(argv[1]) == "-o") {
        replay = true;
        role = argv[2];
    }
    else {
        std::cerr << "usage: ./mpcops_test [-p] [p0|p1|p2]\n";
        return 1;
    }

    boost::asio::io_context io;

    if (role == "p0")
        boost::asio::co_spawn(io, run_p0(io, preprocess, replay), boost::asio::detached);
    else if (role == "p1")
        boost::asio::co_spawn(io, run_p1(io, preprocess, replay), boost::asio::detached);
    else if (role == "p2")
        boost::asio::co_spawn(io, run_p2(io, preprocess, replay), boost::asio::detached);
    else {
        std::cerr << "invalid role\n";
        return 1;
    }

    io.run();
    return 0;
}
