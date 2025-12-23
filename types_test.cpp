#include "types.hpp"
#include <iostream>
#include <iomanip>

// Helper: print __m128i as hex
static void print128(const __m128i& v, const char* label = nullptr) {
    uint64_t lo = _mm_cvtsi128_si64(v);
    uint64_t hi = _mm_extract_epi64(v, 1);
    if (label) std::cout << label << ": ";
    std::cout << "hi=0x" << std::hex << hi
              << " lo=0x" << std::hex << lo << std::dec << std::endl;
}

// Helper: print mX
static void print_mX(const mX& x, const char* label = nullptr) {
    uint64_t lo = x.lo();
    uint64_t hi = x.hi();
    if (label) std::cout << label << ": ";
    std::cout << "hi=0x" << std::hex << hi
              << " lo=0x" << std::hex << lo << std::dec << std::endl;
}

int main() {
    std::cout << "================== TYPES TEST ==================\n";

    // -------------------------------------------------------
    // Basic arithmetic
    // -------------------------------------------------------
    mX x(0xAAAABBBBCCCCDDDDULL, 0x1111222233334444ULL);
    mX y(0x99990000FFFF0000ULL, 0x5555666677778888ULL);

    std::cout << "\nBasic arithmetic:\n";
    print_mX(x + y, "x + y");
    print_mX(x - y, "x - y");
    print_mX(x * y, "x * y");

    // -------------------------------------------------------
    // Bitwise operators
    // -------------------------------------------------------
    std::cout << "\nBitwise ops:\n";
    print_mX(x & y, "x & y");
    print_mX(x ^ y, "x ^ y");
    print_mX(~x, "~x");

    // -------------------------------------------------------
    // Dot product between vector<mX> and vector<mX>
    // -------------------------------------------------------
    std::vector<mX> v1 = { mX(0,1), mX(0,2), mX(0,3), mX(0,4) };
    std::vector<mX> v2 = { mX(0,5), mX(0,6), mX(0,7), mX(0,8) };

    mX dp1 = v1 * v2;
    std::cout << "\nDot product (mX ⋅ mX):\n";
    print_mX(dp1, "v1 ⋅ v2");

    // -------------------------------------------------------
    // Dot product between vector<mX> and vector<uint8_t>
    // -------------------------------------------------------
    std::vector<uint8_t> scalars = { 10, 20, 30, 40 };

    mX dp2 = v1 * scalars;   // mX * uint8_t
    mX dp3 = scalars * v1;   // uint8_t * mX (commutative)

    std::cout << "\nDot product (mX ⋅ uint8_t):\n";
    print_mX(dp2, "v1 ⋅ scalars");
    print_mX(dp3, "scalars ⋅ v1");

    // -------------------------------------------------------
    // Scalar mix test
    // -------------------------------------------------------
    std::cout << "\nScalar mix test:\n";
    print_mX(x + (uint8_t)5, "x + 5");
    print_mX(x * (uint8_t)3, "x * 3");
    print_mX((uint8_t)7 + x, "7 + x");
    print_mX((uint8_t)9 - x, "9 - x");

    std::cout << "\n================== TEST COMPLETE ==================\n";
    return 0;
}
