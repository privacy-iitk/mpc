#ifndef SHARES_H
#define SHARES_H

#include <vector>
#include <cstdint>
#include <cassert>
#include <iostream>
#include <bsd/stdlib.h> // for arc4random_buf

// ----------------------------
// Template for additive shares
// ----------------------------
template<typename T>
struct Share {
    T value;

    Share() = default;
    explicit Share(T val) : value(val) {}

    // Arithmetic operators
    Share<T> operator+(const Share<T> &other) const {
        return Share<T>(value + other.value);
    }

    Share<T> operator-(const Share<T> &other) const {
        return Share<T>(value - other.value);
    }

    Share<T>& operator+=(const Share<T> &other) {
        value += other.value;
        return *this;
    }

    Share<T>& operator-=(const Share<T> &other) {
        value -= other.value;
        return *this;
    }
};

// ---------------------------------
// Split a value into n additive shares
// Uses arc4random_buf() for randomness
// ---------------------------------
template<typename T>
std::vector<Share<T>> share(T secret, size_t n) {
    assert(n >= 2);
    std::vector<Share<T>> shares(n);

    T sum = 0;
    for (size_t i = 0; i < n - 1; ++i) {
        T s;
        arc4random_buf(&s, sizeof(T));
        shares[i] = Share<T>(s);
        sum += s;
    }

    // Last share ensures additive property
    shares[n - 1] = Share<T>(secret - sum);
    return shares;
}

// ---------------------------------
// Reconstruct secret from shares
// ---------------------------------
template<typename T>
T reconstruct(const std::vector<Share<T>>& shares) {
    T sum = 0;
    for (const auto &s : shares) {
        sum += s.value;
    }
    return sum;
}

// ---------------------------------
// Utility: print shares
// ---------------------------------
template<typename T>
void print_shares(const std::vector<Share<T>>& shares) {
    std::cout << "Shares: ";
    for (auto &s : shares) {
        std::cout << s.value << " ";
    }
    std::cout << std::endl;
}

#endif // SHARES_H
