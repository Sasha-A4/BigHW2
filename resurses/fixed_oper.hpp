#pragma once

#include "fixed.hpp"
#include "fast_fixed.hpp"

template<size_t M, size_t L>
Fixed<M,L> operator-(const FastFixed<M,L>& x, const Fixed<M,L>& y) {
    return static_cast<Fixed<M,L>>(x) - y;
}
template<size_t M, size_t L>
Fixed<M,L> operator-(const Fixed<M,L>& x, const FastFixed<M,L>& y) {
    return x - static_cast<Fixed<M,L>>(y);
}
template<size_t M, size_t L>
Fixed<M,L> operator*(const Fixed<M,L>& x, const FastFixed<M,L>& y) {
    return x * static_cast<Fixed<M,L>>(y);
}
template<size_t M, size_t L>
Fixed<M,L> operator+(const Fixed<M,L>& x, const FastFixed<M,L>& y) {
    return x + static_cast<Fixed<M,L>>(y);
}
template<size_t M, size_t L>
Fixed<M,L> operator*(const FastFixed<M,L>& x, const Fixed<M,L>& y) {
    return static_cast<Fixed<M,L>>(x) * y;
}
template<size_t M, size_t L>
Fixed<M,L> operator/(const Fixed<M,L>& x, const FastFixed<M,L>& y) {
    return x / static_cast<Fixed<M,L>>(y);
}
template<size_t M, size_t L>
Fixed<M,L> operator/(const FastFixed<M,L>& x, const Fixed<M,L>& y) {
    return static_cast<Fixed<M,L>>(x) / y;
}
template<size_t M, size_t L>
Fixed<M,L> operator+(const FastFixed<M,L>& x, const Fixed<M,L>& y) {
    return static_cast<Fixed<M,L>>(x) + y;
}
template<size_t M, size_t L>
Fixed<M,L>& operator+=(Fixed<M,L>& x, const FastFixed<M,L>& y) {
    return x += static_cast<Fixed<M,L>>(y);
}
template<size_t M, size_t L>
Fixed<M,L>& operator-=(Fixed<M,L>& x, const FastFixed<M,L>& y) {
    return x -= static_cast<Fixed<M,L>>(y);
}
template<size_t M, size_t L>
Fixed<M,L>& operator/=(Fixed<M,L>& x, const FastFixed<M,L>& y) {
    return x /= static_cast<Fixed<M,L>>(y);
}
template<size_t M, size_t L>
Fixed<M,L>& operator*=(Fixed<M,L>& x, const FastFixed<M,L>& y) {
    return x *= static_cast<Fixed<M,L>>(y);
}
template<size_t M, size_t L>
bool operator<(const Fixed<M,L>& x, const FastFixed<M,L>& y) {
    return x < static_cast<Fixed<M,L>>(y);
}
template<size_t M, size_t L>
bool operator<=(const Fixed<M,L>& x, const FastFixed<M,L>& y) {
    return x <= static_cast<Fixed<M,L>>(y);
}
template<size_t M, size_t L>
bool operator<(const FastFixed<M,L>& x, const Fixed<M,L>& y) {
    return static_cast<Fixed<M,L>>(x) < y;
}
template<size_t M, size_t L>
bool operator<=(const FastFixed<M,L>& x, const Fixed<M,L>& y) {
    return static_cast<Fixed<M,L>>(x) <= y;
}

template<size_t M, size_t L>
bool operator>(const Fixed<M,L>& x, const FastFixed<M,L>& y) {
    return x > static_cast<Fixed<M,L>>(y);
}

template<size_t M, size_t L>
bool operator>(const FastFixed<M,L>& x, const Fixed<M,L>& y) {
    return static_cast<Fixed<M,L>>(x) > y;
}

template<size_t M, size_t L>
bool operator>=(const Fixed<M,L>& x, const FastFixed<M,L>& y) {
    return x >= static_cast<Fixed<M,L>>(y);
}

template<size_t M, size_t L>
bool operator==(const Fixed<M,L>& x, const FastFixed<M,L>& y) {
    return x == static_cast<Fixed<M,L>>(y);
}

template<size_t M, size_t L>
bool operator>=(const FastFixed<M,L>& x, const Fixed<M,L>& y) {
    return static_cast<Fixed<M,L>>(x) >= y;
}
template<size_t M, size_t L>
bool operator==(const FastFixed<M,L>& x, const Fixed<M,L>& y) {
    return static_cast<Fixed<M,L>>(x) == y;
}

