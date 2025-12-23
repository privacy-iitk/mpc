#include <assert.h>
#include <bsd/stdlib.h>
#include <iostream>
#include <emmintrin.h>   // SSE2
#include <tmmintrin.h>   // SSSE3
#include <cstring>
#include <chrono>
#include <array>
#include <type_traits>
#include <fstream>
#include <stdexcept>

#include "dpf.h"
#include "shares.h"

using namespace dpf;

__m128i *proofdb;
__m128i *BB0, *BB1;

typedef std::array<__m128i, 4> leaf_t;

// Helper trait: detect .size()
template <typename T, typename = void>
struct has_size_method : std::false_type {};

template <typename T>
struct has_size_method<T, std::void_t<decltype(std::declval<T>().size())>> : std::true_type {};

// Exclude C-strings and char arrays
template<typename T>
using enable_for_non_string = std::enable_if_t<
    !std::is_same_v<std::decay_t<T>, const char*> &&
    !std::is_same_v<std::decay_t<T>, char*> &&
    !std::is_array_v<T>, int
>;

// Compare two __m128i vectors for equality
inline bool m128i_equal(const __m128i &a, const __m128i &b) {
    __m128i cmp = _mm_cmpeq_epi8(a, b);
    return _mm_movemask_epi8(cmp) == 0xFFFF;
}

// Compare two leaf_t objects
inline bool operator==(const leaf_t &a, const leaf_t &b) {
    for (size_t i = 0; i < a.size(); ++i)
        if (!m128i_equal(a[i], b[i])) return false;
    return true;
}

inline bool operator!=(const leaf_t &a, const leaf_t &b) {
    return !(a == b);
}

// Output operator for leaf_t and __m128i
template<typename T, enable_for_non_string<T> = 0>
std::ostream& operator<<(std::ostream& os, const T& leaf) {
    if constexpr (std::is_same_v<T, __m128i>) {
        alignas(16) uint64_t v[2];
        _mm_store_si128(reinterpret_cast<__m128i*>(v), leaf);
        os << "[" << std::hex << "0x" << v[1] << ", 0x" << v[0] << std::dec << "]";
    } else if constexpr (std::is_arithmetic_v<T>) {
        os << leaf;
    } else if constexpr (has_size_method<T>::value) {
        os << "{ ";
        for (size_t i = 0; i < leaf.size(); ++i) {
            os << leaf[i];
            if (i + 1 != leaf.size()) os << ", ";
        }
        os << " }";
    } else {
        static_assert(sizeof(T) == 0, "Unsupported leaf_t type for operator<<");
    }
    return os;
}

// Save and load DPF key + share value
template<typename leaf_t, typename node_t, typename prgkey_t>
void save_dpf_key(const dpf::dpf_key<leaf_t, node_t, prgkey_t> &key,
                  uint64_t share_val,
                  const std::string &filename)
{
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) throw std::runtime_error("Failed to open file for writing");

    ofs.write(reinterpret_cast<const char*>(&key.nitems), sizeof(key.nitems));
    ofs.write(reinterpret_cast<const char*>(&key.root), sizeof(key.root));

    size_t cw_size = key.cw.size();
    ofs.write(reinterpret_cast<const char*>(&cw_size), sizeof(cw_size));
    ofs.write(reinterpret_cast<const char*>(key.cw.data()), cw_size * sizeof(node_t));

    ofs.write(reinterpret_cast<const char*>(key.finalizer.data()), sizeof(key.finalizer));
    ofs.write(reinterpret_cast<const char*>(&key.prgkey), sizeof(key.prgkey));
    ofs.write(reinterpret_cast<const char*>(&share_val), sizeof(share_val));

    ofs.close();
}

template<typename leaf_t, typename node_t, typename prgkey_t>
std::pair<dpf::dpf_key<leaf_t, node_t, prgkey_t>, uint64_t>
load_dpf_key(const std::string &filename)
{
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) throw std::runtime_error("Failed to open file for reading");

    size_t nitems;
    node_t root;
    size_t cw_size;

    ifs.read(reinterpret_cast<char*>(&nitems), sizeof(nitems));
    ifs.read(reinterpret_cast<char*>(&root), sizeof(root));
    ifs.read(reinterpret_cast<char*>(&cw_size), sizeof(cw_size));

    std::vector<node_t> cw(cw_size);
    ifs.read(reinterpret_cast<char*>(cw.data()), cw_size * sizeof(node_t));

    std::array<node_t, dpf::dpf_key<leaf_t, node_t, prgkey_t>::nodes_per_leaf> finalizer;
    ifs.read(reinterpret_cast<char*>(finalizer.data()), sizeof(finalizer));

    prgkey_t prgkey;
    ifs.read(reinterpret_cast<char*>(&prgkey), sizeof(prgkey));

    uint64_t share_val;
    ifs.read(reinterpret_cast<char*>(&share_val), sizeof(share_val));

    return {
        dpf::dpf_key<leaf_t, node_t, prgkey_t>(nitems, root, cw, finalizer, prgkey),
        share_val
    };
}

// Save/load array functions
template<typename T>
void save_array(const std::string &filename, const T *data, size_t nitems) {
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs)
        throw std::runtime_error("Failed to open file for writing: " + filename);
    ofs.write(reinterpret_cast<const char*>(&nitems), sizeof(nitems));
    ofs.write(reinterpret_cast<const char*>(data), nitems * sizeof(T));
    ofs.close();
}

template<typename T>
T* load_array(const std::string &filename, size_t &nitems_out) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs)
        throw std::runtime_error("Failed to open file for reading: " + filename);
    ifs.read(reinterpret_cast<char*>(&nitems_out), sizeof(nitems_out));
    T* data = new T[nitems_out];
    ifs.read(reinterpret_cast<char*>(data), nitems_out * sizeof(T));
    ifs.close();
    return data;
}

// Main
int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <log2(nitems)>\n";
        return 1;
    }

    uint64_t log_nitems = std::stoull(argv[1]);
    uint64_t nitems = 1ULL << log_nitems;
    uint64_t target_ind = 555;

    AES_KEY prgkey;
    leaf_t target_value;

    leaf_t *output0 = new leaf_t[nitems];
    leaf_t *output1 = new leaf_t[nitems];
    uint8_t *_t0 = new uint8_t[nitems];
    uint8_t *_t1 = new uint8_t[nitems];

    auto [dpfkey0, dpfkey1] = dpf_key<leaf_t, __m128i, AES_KEY>::gen(prgkey, nitems, target_ind, target_value);

    auto shares = share<uint64_t>(target_ind, 2);
    uint64_t share0 = shares[0].value;
    uint64_t share1 = shares[1].value;
    print_shares(shares);

    uint64_t rec = reconstruct(shares);
    std::cout << "Reconstructed secret: " << rec << std::endl;

    // Save both key and share
    save_dpf_key(dpfkey0, share0, "dpfkey0.bin");
    save_dpf_key(dpfkey1, share1, "dpfkey1.bin");

    // Load them back
    auto [loaded0, share0_loaded] = load_dpf_key<leaf_t, __m128i, AES_KEY>("dpfkey0.bin");
    auto [loaded1, share1_loaded] = load_dpf_key<leaf_t, __m128i, AES_KEY>("dpfkey1.bin");

    std::cout << "Loaded share0: " << share0_loaded << std::endl;
    std::cout << "Loaded share1: " << share1_loaded << std::endl;

    uint64_t recc = share0_loaded + share1_loaded;
    std::cout << "rec = " << recc << std::endl;

    __evalinterval(loaded0, 0, nitems - 1, output0, _t0);
    __evalinterval(loaded1, 0, nitems - 1, output1, _t1);

    save_array("output0.bin", output0, nitems);
    save_array("output1.bin", output1, nitems);

    leaf_t *loaded_outs0 = load_array<leaf_t>("output0.bin", nitems);
    leaf_t *loaded_outs1 = load_array<leaf_t>("output1.bin", nitems);

    for (size_t j = 0; j < nitems; ++j) {
        if (output0[j] != output1[j]) {
            std::cout << "Inequality at: " << j << std::endl;
            std::cout << output0[j] << " <> " << output1[j] << std::endl;
        }
        if (loaded_outs0[j] != loaded_outs1[j]) {
            std::cout << "Inequality (loaded) at: " << j << std::endl;
            std::cout << output0[j] << " <> " << output1[j] << std::endl;
        }
    }



 
    uint64_t a0 = 335530;
    uint64_t a1 = static_cast<uint64_t>(-335530);
    
    // Random mask r
    uint64_t r;
    arc4random_buf(&r, sizeof(uint64_t));

    // XOR shares of r
    uint64_t r0, r1;
    arc4random_buf(&r0, sizeof(uint64_t));
    r1 = r ^ r0;

    // Additive shares of r
    uint64_t r0_ashare;
    arc4random_buf(&r0_ashare, sizeof(uint64_t));
    uint64_t r1_ashare = r - r0_ashare;

    // Verify shares
    std::cout << "XOR:       " << (r0 ^ r1) << std::endl;
    std::cout << "Additive:  " << (r0_ashare + r1_ashare) << std::endl;

    // Compute masked sum
    uint64_t c = a0 + a1 + r0_ashare + r1_ashare;

    // Bitwise equality check with r
    uint64_t eqz = 1;
    
    for (int i = 0; i < 64; i++) { 
        uint64_t c_i = (c >> i) & 1ULL; 
        uint64_t r0_i = (r0 >> i) & 1ULL; 
        uint64_t r1_i = (r1 >> i) & 1ULL; 
        uint64_t diff = c_i ^ (r0_i ^ r1_i); // XOR //std::cout << "diff = " << diff << std::endl; 
        eqz &= (1ULL - diff); 
    }
    std::cout << "EQZ = " << eqz << std::endl;

    return 0;
}
