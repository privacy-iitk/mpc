#include "shares.hpp"
#include <iostream>

int main() {
    std::cout << "================== SHARES TEST ==================\n";

    using T = uint64_t;
    T secret = 0xDEADBEEFCAFEBABEULL;
    std::cout << "Original secret: 0x" << std::hex << secret << std::dec << "\n";

    // -------------------------------------------------
    // Test additive sharing
    // -------------------------------------------------
    auto [a0, a1] = share_secret_additive(secret);
    T reconstructed_add = static_cast<T>(a0.val + a1.val);


    std::cout << "a =  " << std::hex << a0 << " " << a1 << std::endl;
    std::cout << "\nAdditive shares:\n";
    std::cout << "a0 = 0x" << std::hex << a0.val << "\n";
    std::cout << "a1 = 0x" << std::hex << a1.val << "\n";
    std::cout << "Reconstructed = 0x" << reconstructed_add << std::dec << "\n";

    if (reconstructed_add == secret)
        std::cout << "✅ Additive share reconstruction succeeded\n";
    else
        std::cout << "❌ Additive share reconstruction failed\n";

    // -------------------------------------------------
    // Test XOR sharing
    // -------------------------------------------------
    auto [x0, x1] = share_secret_xor(secret);
    T reconstructed_xor = static_cast<T>(x0.val ^ x1.val);

    std::cout << "\nXOR shares:\n";
    std::cout << "x0 = 0x" << std::hex << x0.val << "\n";
    std::cout << "x1 = 0x" << std::hex << x1.val << "\n";
    std::cout << "Reconstructed = 0x" << reconstructed_xor << std::dec << "\n";

    if (reconstructed_xor == secret)
        std::cout << "✅ XOR share reconstruction succeeded\n";
    else
        std::cout << "❌ XOR share reconstruction failed\n";

    // -------------------------------------------------
    // Test vector XOR shares
    // -------------------------------------------------
    std::vector<uint64_t> vec = {1, 2, 3, 4};
    std::vector<uint64_t> mask = {5, 6, 7, 8};

    XorShareVector<uint64_t> v1(vec);
    XorShareVector<uint64_t> v2(mask);

    auto vxor = v1 ^ v2;

    std::cout << "\nVector XOR test:\n";
    std::cout << "v1 = " << v1 << "\n";
    std::cout << "v2 = " << v2 << "\n";
    std::cout << "v1 ^ v2 = " << vxor << "\n";


    std::vector<uint64_t> vv1 = {10, 20, 30};
    std::vector<uint64_t> vv2 = {1, 2, 3};

    AdditiveShareVector<uint64_t> a(vv1);
    AdditiveShareVector<uint64_t> b(vv2);

    auto c = a + b;
    auto d = a - b;

    std::cout << "a + b = " << c << "\n";
    std::cout << "a - b = " << d << "\n";


    std::cout << "\n================== TEST COMPLETE ==================\n";
    return 0;
}
