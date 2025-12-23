// locoram.tpp
#pragma once

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <cstdint>
#include <type_traits>

// Implementation for Locoram template

template<typename leaf_t, typename node_t, typename prgkey_t>
Locoram<leaf_t, node_t, prgkey_t>::Locoram() = default;

template<typename leaf_t, typename node_t, typename prgkey_t>
Locoram<leaf_t, node_t, prgkey_t>::Locoram(size_t n, const prgkey_t& key)
    : secret_share(n),
      blind(n),
      peer_blinded(n),
      sbv_secret_share(n),
      sbv_blind(n),
      sbv_peer_blinded(n),
      nitems(n),
      prgkey(key)
{
    for (size_t i = 0; i < n; ++i) {
        set_zero(secret_share[i]);
        set_zero(blind[i]);
        set_zero(peer_blinded[i]);
        set_zero(sbv_secret_share[i]);
        set_zero(sbv_blind[i]);
        set_zero(sbv_peer_blinded[i]);
    }
    set_zero(cancellation_term);
}

template<typename leaf_t, typename node_t, typename prgkey_t>
void Locoram<leaf_t, node_t, prgkey_t>::set_share(const std::vector<leaf_t>& share) {
    if (share.size() != nitems)
        throw std::runtime_error("Size mismatch in Locoram::set_share");
    secret_share = share;
}

template<typename leaf_t, typename node_t, typename prgkey_t>
void Locoram<leaf_t, node_t, prgkey_t>::dump(const std::string& label) const {
    std::cout << "==== Locoram Dump: " << label << " ====" << std::endl;
    for (size_t i = 0; i < secret_share.size(); ++i) {
        if constexpr (std::is_arithmetic_v<leaf_t>) {
            std::cout << "[" << i << "] = " << +secret_share[i] << std::endl;
        } else {
            std::cout << "[" << i << "] = " << secret_share[i][0] << " "
                      << secret_share[i][1] << std::endl;
        }
    }
    std::cout << "====================================" << std::endl;
}

template<typename leaf_t, typename node_t, typename prgkey_t>
boost::asio::awaitable<void>
Locoram<leaf_t, node_t, prgkey_t>::reconstruct_with_ctx(MPCContext* ctx,
                                                        const std::string& label) {
    std::cout << "REC" << std::endl;
    if (ctx == nullptr)
        throw std::runtime_error("Missing MPCContext in Locoram::reconstruct_with_ctx");

    Role role = ctx->role;

    if (role == Role::P2) {
        std::cout << "[Locoram::reconstruct] P2: nothing to do for reconstruction" << std::endl;
        co_return;
    }

    //NetPeer& self = ctx->self;

    NetPeer* peer_ptr = nullptr;
    if (role == Role::P0)
        peer_ptr = ctx->peer1 ? ctx->peer1 : ctx->peer0;
    else if (role == Role::P1)
        peer_ptr = ctx->peer0 ? ctx->peer0 : ctx->peer1;

    if (!peer_ptr)
        throw std::runtime_error("Locoram::reconstruct_with_ctx: missing peer NetPeer pointer (P0↔P1)");

    NetPeer& peer = *peer_ptr;

    size_t n = nitems;
    size_t bytes = n * sizeof(leaf_t);

    static_assert(std::is_trivially_copyable_v<leaf_t>,
                  "leaf_t must be trivially copyable for reconstruct serialization");

    std::vector<uint8_t> buf(bytes);
    if (bytes)
        std::memcpy(buf.data(), reinterpret_cast<const uint8_t*>(secret_share.data()), bytes);

    std::vector<uint8_t> buf_peer;

    if (role == Role::P0) {
        std::cout << "[Locoram::reconstruct] P0: sending " << bytes << " bytes to peer" << std::endl;
        co_await (peer << buf);
        std::cout << "[Locoram::reconstruct] P0: waiting to receive peer buffer" << std::endl;
        co_await (peer >> buf_peer);
    } else {
        std::cout << "[Locoram::reconstruct] P1: waiting to receive peer buffer" << std::endl;
        co_await (peer >> buf_peer);
        std::cout << "[Locoram::reconstruct] P1: sending " << bytes << " bytes to peer" << std::endl;
        co_await (peer << buf);
    }

    std::cout << "------>>> received " << buf_peer.size() << " bytes, expected " << bytes << std::endl;

    if (buf_peer.size() != bytes) {
        std::cout << "------>>> size mismatch: got " << buf_peer.size() << " expected " << bytes << std::endl;
        throw std::runtime_error("Locoram::reconstruct_with_ctx: received buffer size mismatch");
    }

    const leaf_t* peer_data = reinterpret_cast<const leaf_t*>(buf_peer.data());
    std::vector<leaf_t> reconstructed(n);
    for (size_t i = 0; i < n; ++i)
        reconstructed[i] = secret_share[i] ^ peer_data[i];

    std::cout << "==== Locoram Reconstruction: " << label << " ====" << std::endl;
    for (size_t i = 0; i < reconstructed.size(); ++i) {
        if constexpr (std::is_arithmetic_v<leaf_t>) {
             std::cout << "[" << i << "] = " << +reconstructed[i] << std::endl;
        } else {
            std::cout << "[" << i << "] = " << reconstructed[i][0]
                      << " " << reconstructed[i][1] << std::endl;
        }
    }
    std::cout << "====================================" << std::endl;

    co_return;
}

// ------------------------------------------------------------------------
// Ref methods
// ------------------------------------------------------------------------
template<typename leaf_t, typename node_t, typename prgkey_t>
Locoram<leaf_t, node_t, prgkey_t>::Ref::Ref(Locoram& p, size_t i) : parent(p), idx(i) {}


// ============================================================================
// Helper: receive and evaluate DPF key from P2
// ============================================================================
template<typename leaf_t, typename node_t, typename prgkey_t>
boost::asio::awaitable<void>
locoram_recv_and_eval_read(NetPeer& self,
                           Locoram<leaf_t, node_t, prgkey_t>& parent,
                           const char* who)
{
    std::vector<uint8_t> buffer;
    std::cout << "[" << who << "] Waiting to receive DPF key from P2...\n";
    co_await (self >> buffer);
    std::cout << "[" << who << "] Received " << buffer.size() << " bytes from P2\n";

    auto key = dpf::deserialize_dpf_key<leaf_t, __m128i, AES_KEY>(buffer.data(), buffer.size());
    size_t nitems = parent.nitems;

    std::vector<leaf_t> output(nitems);
    std::vector<uint8_t> flags(nitems);
    __evalinterval(key, 0, nitems - 1, output.data(), flags.data());
    
    for (size_t j = 0; j < nitems; ++j)
        std::cout << j << " --->> " << output[j][0] << std::endl; 
      //  parent.sbv_secret_share[j] = parent.sbv_secret_share[j] + flags[j];

    co_return;
}



// ============================================================================
// Ref::read for XShare
// ============================================================================
// template<typename leaf_t, typename node_t, typename prgkey_t>
// boost::asio::awaitable<leaf_t>
// Locoram<leaf_t, node_t, prgkey_t>::Ref::read(const XShare<uint64_t>& shared_index)
// {
//     std::cout << "READ!!! (in the XOR Secret Sharing Scheme!!)" << std::endl; 
//     if (shared_index.ctx == nullptr)
//         throw std::runtime_error("Missing MPCContext in Locoram::Ref::read");

//     MPCContext* ctx = shared_index.ctx;
//     Role role = ctx->role;

//     if (idx >= parent.nitems)
//         throw std::runtime_error("Index out of bounds in Locoram::Ref::read");

//     leaf_t local_value = parent.secret_share[idx];
//     leaf_t result = local_value;

//     switch (role) {
//         case Role::P0:
//         case Role::P1:
//             co_await locoram_recv_and_eval_read<leaf_t, node_t, prgkey_t>(ctx->self, parent, role == Role::P0 ? "P0" : "P1");
//             break;

//         case Role::P2: {
//             uint64_t nitems = parent.nitems;
//             uint64_t alpha;
//             arc4random_buf(&alpha, sizeof(uint64_t));
//             std::cout << "nitems = " << nitems << std::endl;
//             alpha %= nitems;
//             alpha = 5;
//             std::cout << "alpha = " << alpha << std::endl;
//             leaf_t target_value{};
//             target_value[0] = 100;
//             target_value[1] = 200;
//             AES_KEY prgkey;
//             auto [k0, k1] = dpf::dpf_key<leaf_t, __m128i, AES_KEY>::gen(prgkey, nitems, alpha, target_value);

//             if (!ctx->peer0 || !ctx->peer1)
//                 throw std::runtime_error("P2 missing peer connections");

//             NetPeer& p0 = *(ctx->peer0);
//             NetPeer& p1 = *(ctx->peer1);

//             auto buf0 = dpf::serialize_dpf_key(k0);
//             auto buf1 = dpf::serialize_dpf_key(k1);

//             std::cout << "[P2] Sending keys to P0/P1 (alpha=" << alpha << ")\n";
//             co_await ((p0 << buf0) && (p1 << buf1));
//             co_return leaf_t{};
//         }

//         default:
//             throw std::runtime_error("Unknown role in Locoram::Ref::read");
//     }

//     co_return result;
// }



template<typename leaf_t, typename node_t, typename prgkey_t>
boost::asio::awaitable<leaf_t>
Locoram<leaf_t, node_t, prgkey_t>::Ref::read(const AShare<uint64_t>& shared_index)
{


    std::cout << "Read Operation for AShare" << std::endl;
    if (shared_index.ctx == nullptr)
        throw std::runtime_error("Missing MPCContext in Locoram::Ref::read");

    MPCContext* ctx = shared_index.ctx;
    Role role = ctx->role;

    const size_t n = parent.nitems;
 
    // ---------------------------------------------------------------------
    // Step 1: P0 / P1 receive DPF key and evaluate full vector
    // ---------------------------------------------------------------------
    std::vector<uint64_t> w_vals(n, 0);  // DPF output
    std::vector<uint64_t> x_vals(n, 0);  // database shares

    if (role == Role::P0 || role == Role::P1) {
        std::vector<uint8_t> buffer;
        co_await (ctx->self >> buffer);

        auto key =
            dpf::deserialize_dpf_key<uint64_t, __m128i, AES_KEY>(
                buffer.data(), buffer.size());

        std::vector<uint64_t> eval(n);
        std::vector<uint8_t> flags(n);

        __evalinterval(key, 0, n - 1, eval.data(), flags.data()); 

        for (size_t i = 0; i < n; ++i) {
            w_vals[i] = flags[i];

            parent.sbv_secret_share[i] ^= w_vals[i];
            std::cout << "sanity --> " << parent.sbv_secret_share[i] << " + " << (int) w_vals[i] << std::endl;
            parent.sbv_blind[i] ^= w_vals[i];
            // parent.sbv_peer_blinded[i] ^= flags[i];
            // parent.sbv_peer_blinded[i] ^= flags[i];
            // interpret leaf_t as holding additive share in lane 0
            x_vals[i] = reinterpret_cast<const uint64_t*>(&parent.secret_share[i])[0];
        }
        
         parent.cancellation_term ^= (parent.blind & w_vals);
    }

    AdditiveShareVector<uint64_t> w(std::move(w_vals), ctx);
    AdditiveShareVector<uint64_t> x(std::move(x_vals), ctx);

    // ---------------------------------------------------------------------
    // Step 2: P2 generates and distributes DPF keys
    // ---------------------------------------------------------------------
    if (role == Role::P2) {
        uint64_t alpha;
        arc4random_buf(&alpha, sizeof(uint64_t));
        alpha %= n;
        alpha = 5;
        std::cout << "alpha = " << alpha << std::endl;

        AES_KEY prgkey;
        auto [k0, k1] =
            dpf::dpf_key<uint64_t, __m128i, AES_KEY>::gen(
                prgkey, n, alpha, uint64_t(1));

        auto buf0 = dpf::serialize_dpf_key(k0);
        auto buf1 = dpf::serialize_dpf_key(k1);

        if (!ctx->peer0 || !ctx->peer1)
            throw std::runtime_error("P2 missing peer connections");

        co_await (
            (*(ctx->peer0) << buf0) &&
            (*(ctx->peer1) << buf1)
        );
    }

    // ---------------------------------------------------------------------
    // Step 3: MPC dot product
    // ---------------------------------------------------------------------
 
 
leaf_t result{};
   

    auto tmp = (parent.sbv_secret_share ^ parent.sbv_peer_blinded);
    result = (parent.secret_share & (parent.sbv_secret_share ^ parent.sbv_peer_blinded)) 
                ^ (parent.peer_blinded & parent.sbv_blind) ^ parent.cancellation_term;
     leaf_t result_a = (parent.secret_share & (parent.sbv_secret_share ^ parent.sbv_peer_blinded));
    
    for(size_t i = 0; i < 8; ++i) std::cout << "tmp -->> " << tmp[i] << std::endl;
    for(size_t i = 0; i < 8; ++i) std::cout << "parent.secret_share -->> " << parent.secret_share[i] << std::endl;

      for(size_t i = 0; i < 8; ++i) std::cout << "parent.sbv_secret_share -->> " << parent.sbv_secret_share[i] << std::endl;
    std::cout << "result = " << result << std::endl;
    std::cout << "result_a = " << result_a << std::endl;
    std::cout << "parent.cancellation_term = " << parent.cancellation_term << std::endl;
    co_return result;
}



// ============================================================================
// Ref::read for XShare
// ============================================================================
template<typename leaf_t, typename node_t, typename prgkey_t>
boost::asio::awaitable<leaf_t>
Locoram<leaf_t, node_t, prgkey_t>::Ref::read(const XShare<uint64_t>& shared_index)
{
    std::cout << "It is the XOR Read Operation that is being called!" << std::endl;
 
    if (shared_index.ctx == nullptr)
        throw std::runtime_error("Missing MPCContext in Locoram::Ref::read");

    MPCContext* ctx = shared_index.ctx;
    Role role = ctx->role;

    const size_t n = parent.nitems;
    parent.read_inds.push_back(shared_index);
    // ---------------------------------------------------------------------
    // Step 1: P0 / P1 receive DPF key and evaluate full vector
    // ---------------------------------------------------------------------
    std::vector<uint64_t> w_vals(n, 0);  // DPF output
    std::vector<uint64_t> x_vals(n, 0);  // database shares

    if (role == Role::P0 || role == Role::P1) {
        std::vector<uint8_t> buffer;
        co_await (ctx->self >> buffer);

        auto key =
            dpf::deserialize_dpf_key<uint64_t, __m128i, AES_KEY>(
                buffer.data(), buffer.size());

        std::vector<uint64_t> eval(n);
        std::vector<uint8_t> flags(n);

        __evalinterval(key, 0, n - 1, eval.data(), flags.data());

        for (size_t i = 0; i < n; ++i) {
            w_vals[i] = flags[i];

            parent.sbv_secret_share[i] ^= w_vals[i];
            std::cout << "sanity --> " << parent.sbv_secret_share[i] << " + " << (int) w_vals[i] << std::endl;
            parent.sbv_blind[i] ^= w_vals[i];
            // parent.sbv_peer_blinded[i] ^= flags[i];
            // parent.sbv_peer_blinded[i] ^= flags[i];
            // interpret leaf_t as holding additive share in lane 0
            x_vals[i] = reinterpret_cast<const uint64_t*>(&parent.secret_share[i])[0];
        }

        parent.cancellation_term ^= (parent.blind & w_vals) ^ parent.update_vals[0].val;
    }

    AdditiveShareVector<uint64_t> w(std::move(w_vals), ctx);
    AdditiveShareVector<uint64_t> x(std::move(x_vals), ctx);

    // ---------------------------------------------------------------------
    // Step 2: P2 generates and distributes DPF keys
    // ---------------------------------------------------------------------
    if (role == Role::P2) {
        uint64_t alpha;
        arc4random_buf(&alpha, sizeof(uint64_t));
        alpha %= n;
        alpha = 5;
        std::cout << "alpha = " << alpha << std::endl;

        AES_KEY prgkey;
        auto [k0, k1] =
            dpf::dpf_key<uint64_t, __m128i, AES_KEY>::gen(
                prgkey, n, alpha, uint64_t(1));

        auto buf0 = dpf::serialize_dpf_key(k0);
        auto buf1 = dpf::serialize_dpf_key(k1);

        if (!ctx->peer0 || !ctx->peer1)
            throw std::runtime_error("P2 missing peer connections");

        co_await (
            (*(ctx->peer0) << buf0) &&
            (*(ctx->peer1) << buf1)
        );
    }

    // ---------------------------------------------------------------------
    // Step 3: MPC dot product
    // ---------------------------------------------------------------------
 
 
leaf_t result{};
   
  auto tmp = (parent.sbv_secret_share ^ parent.sbv_peer_blinded);
    result = (parent.secret_share & (parent.sbv_secret_share ^ parent.sbv_peer_blinded)) 
                ^ (parent.peer_blinded & parent.sbv_blind) ^ parent.cancellation_term;
     leaf_t result_a = (parent.secret_share & (parent.sbv_secret_share ^ parent.sbv_peer_blinded));
    
    for(size_t i = 0; i < 8; ++i) std::cout << "tmp -->> " << tmp[i] << std::endl;
    for(size_t i = 0; i < 8; ++i) std::cout << "parent.secret_share -->> " << parent.secret_share[i] << std::endl;

      for(size_t i = 0; i < 8; ++i) std::cout << "parent.sbv_secret_share -->> " << parent.sbv_secret_share[i] << std::endl;
    std::cout << "result = " << result << std::endl;
    std::cout << "result_a = " << result_a << std::endl;
    std::cout << "parent.cancellation_term = " << parent.cancellation_term << std::endl;

    co_return result;
}

// ============================================================================
// Helper for update
// ============================================================================
template<typename leaf_t, typename node_t, typename prgkey_t>
boost::asio::awaitable<void>
locoram_recv_and_apply_p0p1(NetPeer& self, Locoram<leaf_t, node_t, prgkey_t>& parent) {
    std::vector<uint8_t> buffer;
    std::cout << "[P?] Waiting to receive DPF key from P2...\n";
    co_await (self >> buffer);
    std::cout << "[P?] Received " << buffer.size() << " bytes from P2\n";

    auto key = dpf::deserialize_dpf_key<leaf_t, __m128i, AES_KEY>(buffer.data(), buffer.size());
    size_t nitems = parent.nitems;
    if (nitems == 0) co_return;

    std::vector<leaf_t> output(nitems);
    std::vector<uint8_t> flags(nitems);
    __evalinterval(key, 0, nitems - 1, output.data(), flags.data());
std::vector<leaf_t> w_vals(nitems);
    for (size_t j = 0; j < nitems; ++j)
    {   
        w_vals[j] = flags[j];
        parent.secret_share[j] = parent.secret_share[j] ^ output[j];
        parent.blind[j] = parent.blind[j] ^ output[j];
    }
   // parent.cancellation_term ^= (output & parent.sbv_blind);
    co_return; 
}

// ============================================================================
// Ref::update
// ============================================================================
template<typename leaf_t, typename node_t, typename prgkey_t>
boost::asio::awaitable<leaf_t>
Locoram<leaf_t, node_t, prgkey_t>::Ref::update(const XShare<uint64_t>& shared_index,
                                               const XShare<uint64_t>& value_share) {



    std::cout << "UPDATE ------ is being called!!! " << std::endl;
    if (shared_index.ctx == nullptr)
        throw std::runtime_error("Missing MPCContext in Locoram::Ref::update");

    parent.update_inds.push_back(shared_index);
    parent.update_vals.push_back(value_share);
    
    MPCContext* ctx = shared_index.ctx;
    Role role = ctx->role;

    leaf_t local_value = parent.secret_share[idx];
    leaf_t result = local_value;

    switch (role) {
        case Role::P0:
        case Role::P1:
            co_await locoram_recv_and_apply_p0p1<leaf_t, node_t, prgkey_t>(ctx->self, parent);
            
            break;

        case Role::P2: {
            uint64_t nitems = parent.nitems;
            uint64_t alpha;
            arc4random_buf(&alpha, sizeof(uint64_t));
            alpha %= nitems;
            alpha = 2;//shared_index.val;

            leaf_t target_value{};
            AES_KEY prgkey;
            arc4random_buf(&target_value, sizeof(leaf_t));
            target_value = 100;
    
            auto [k0, k1] = dpf::dpf_key<leaf_t, __m128i, AES_KEY>::gen(prgkey, nitems, alpha, target_value);

            if (!ctx->peer0 || !ctx->peer1)
                throw std::runtime_error("P2 missing peer connections");

            NetPeer& p0 = *(ctx->peer0);
            NetPeer& p1 = *(ctx->peer1);

            auto buf0 = dpf::serialize_dpf_key(k0);
            auto buf1 = dpf::serialize_dpf_key(k1);

            std::cout << "[P2] Sending keys to P0/P1 (alpha=" << alpha << ")\n";
            co_await ((p0 << buf0) && (p1 << buf1));
            break;
        }

        default:
            throw std::runtime_error("Unknown role in Locoram::Ref::update");
    }

    co_return result;
}

// ------------------------------------------------------------------------
// operator[] and high-level update
// ------------------------------------------------------------------------
template<typename leaf_t, typename node_t, typename prgkey_t>
typename Locoram<leaf_t, node_t, prgkey_t>::Ref
Locoram<leaf_t, node_t, prgkey_t>::operator[](size_t idx) {
    if (idx >= nitems)
        throw std::runtime_error("Index out of bounds in Locoram::operator[]");
    return Ref(*this, idx);
}

template<typename leaf_t, typename node_t, typename prgkey_t>
boost::asio::awaitable<leaf_t>
Locoram<leaf_t, node_t, prgkey_t>::operator[](const AShare<uint64_t>& shared_index) {
    if (shared_index.ctx == nullptr)
        throw std::runtime_error("AShare index missing MPCContext in Locoram::operator[]");
    size_t local_idx = static_cast<size_t>(shared_index.val);
    auto ref = (*this)[local_idx];
    co_return co_await ref.read(shared_index);
}

template<typename leaf_t, typename node_t, typename prgkey_t>
boost::asio::awaitable<leaf_t>
Locoram<leaf_t, node_t, prgkey_t>::operator[](const XShare<uint64_t>& shared_index) {
    if (shared_index.ctx == nullptr)
        throw std::runtime_error("XShare index missing MPCContext in Locoram::operator[]");
    size_t local_idx = static_cast<size_t>(shared_index.val);
    auto ref = (*this)[local_idx];
    co_return co_await ref.read(shared_index);
}

template<typename leaf_t, typename node_t, typename prgkey_t>
boost::asio::awaitable<void>
Locoram<leaf_t, node_t, prgkey_t>::update(const XShare<uint64_t>& shared_index,
                                          const XShare<uint64_t>& value_share) {
    if (shared_index.ctx == nullptr)
        throw std::runtime_error("AShare index missing MPCContext in Locoram::update");

    size_t local_idx = static_cast<size_t>(shared_index.val);
    if (local_idx >= nitems)
        throw std::runtime_error("Index out of bounds in Locoram::update");

    auto ref = (*this)[local_idx];
    co_await ref.update(shared_index, value_share);
}
