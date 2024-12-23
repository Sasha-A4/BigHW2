#pragma once

#include <iostream>
#include <type_traits>
#include <cstdint>
#include <cmath>

template<size_t P, size_t Q>
struct FastFixed;

template<size_t P, size_t Q>
struct Fixed {
    static_assert(P > Q, "P must be large than Q");
    static_assert(P <= 64, "P must be under than or equal to 64");

    static constexpr size_t bits = P;
    static constexpr size_t frac = Q;

    using StorageType = typename std::conditional<P <= 32, int32_t, int64_t>::type;

    static constexpr Fixed from_raw(StorageType x) {
        Fixed res;
        res.value = x;
        return res;
    }

    constexpr Fixed(int n = 0): value(static_cast<StorageType>(n) << Q) {}
    constexpr Fixed(float f): value(f * (StorageType(1) << Q)) {}
    constexpr Fixed(double f): value(f * (StorageType(1) << Q)) {}


    StorageType value;


    explicit operator float() const { return value / float(StorageType(1) << Q); }
    explicit operator double() const { return value / double(StorageType(1) << Q); }
    auto operator<=>(const Fixed&) const = default;
    bool operator==(const Fixed&) const = default;

    friend Fixed operator/(Fixed a, int b) {
        return Fixed::from_raw(a.value / b);
    }
    friend Fixed operator*(Fixed a, int b) {
        return Fixed::from_raw(a.value * b);
    }
    friend Fixed operator*(int a, Fixed b) {
        return b * a;
    }

    template<size_t P2, size_t Q2>
    explicit operator Fixed<P2,Q2>() const {
        if constexpr (Q2 >= Q) {
            return Fixed<P2,Q2>::from_raw(static_cast<typename Fixed<P2,Q2>::StorageType>(value) << (Q2 - Q));
        } else {
            constexpr size_t shift = Q - Q2;
            if constexpr (shift >= P2) {
                auto temp = value >> (shift - P2 + 1);
                return Fixed<P2,Q2>::from_raw(static_cast<typename Fixed<P2,Q2>::StorageType>(temp) >> 1);
            } else {
                return Fixed<P2,Q2>::from_raw(static_cast<typename Fixed<P2,Q2>::StorageType>(value) >> shift);
            }
        }
    }

    template<size_t P2, size_t Q2>
    explicit operator FastFixed<P2,Q2>() const {
        if constexpr (Q2 >= Q) {
            return FastFixed<P2,Q2>::from_raw(static_cast<typename FastFixed<P2,Q2>::StorageType>(value) << (Q2 - Q));
        } else {
            constexpr size_t shift = Q - Q2;
            if constexpr (shift >= P2) {
                auto temp = value >> (shift - P2 + 1);
                return FastFixed<P2,Q2>::from_raw(static_cast<typename FastFixed<P2,Q2>::StorageType>(temp) >> 1);
            } else {
                return FastFixed<P2,Q2>::from_raw(static_cast<typename FastFixed<P2,Q2>::StorageType>(value) >> shift);
            }
        }
    }
};

template<size_t P, size_t Q>
Fixed<P,Q> operator+(Fixed<P,Q> x, Fixed<P,Q> y) {
    return Fixed<P,Q>::from_raw(x.value + y.value);
}

template<size_t P, size_t Q>
Fixed<P,Q> operator-(Fixed<P,Q> x, Fixed<P,Q> y) {
    return Fixed<P,Q>::from_raw(x.value - y.value);
}

template<size_t P, size_t Q>
Fixed<P,Q> operator*(Fixed<P,Q> x, Fixed<P,Q> y) {
    using ST = typename Fixed<P,Q>::StorageType;
    if constexpr (P <= 32) {
        return Fixed<P,Q>::from_raw((static_cast<int64_t>(x.value) * y.value) >> Q);
    } else {
        ST high = (x.value >> Q) * y.value;
        ST low = (x.value & ((ST(1) << Q) - 1)) * y.value >> Q;
        return Fixed<P,Q>::from_raw(high + low);
    }
}

template<size_t P, size_t Q>
Fixed<P,Q> operator/(Fixed<P,Q> x, Fixed<P,Q> y) {
    using ST = typename Fixed<P,Q>::StorageType;
    if constexpr (P <= 32) {
        return Fixed<P,Q>::from_raw((static_cast<int64_t>(x.value) << Q) / y.value);
    } else {
        return Fixed<P,Q>::from_raw((x.value << Q) / y.value);
    }
}
template<size_t P, size_t Q>
Fixed<P,Q> abs(Fixed<P,Q> x) {
    Fixed<P,Q> res = x;
    if (res.value < 0) res.value = -res.value;
    return res;
}

template<size_t P, size_t Q>
Fixed<P,Q> &operator+=(Fixed<P,Q> &x, Fixed<P,Q> y) { return x = x + y; }
template<size_t P, size_t Q>
Fixed<P,Q> &operator-=(Fixed<P,Q> &x, Fixed<P,Q> y) { return x = x - y; }
template<size_t P, size_t Q>
Fixed<P,Q> &operator*=(Fixed<P,Q> &x, Fixed<P,Q> y) { return x = x * y; }
template<size_t P, size_t Q>
Fixed<P,Q> &operator/=(Fixed<P,Q> &x, Fixed<P,Q> y) { return x = x / y; }

template<size_t P, size_t Q>
Fixed<P,Q> operator+(double a, Fixed<P,Q> x) {
    return Fixed<P,Q>(a) + x;
}
template<size_t P, size_t Q>
Fixed<P,Q>& operator-=(Fixed<P,Q>& x, float a) {
    return x = x - Fixed<P,Q>(a);
}
template<size_t P, size_t Q>
Fixed<P,Q>& operator+=(Fixed<P,Q>& x, float a) {
    return x = x + Fixed<P,Q>(a);
}
template<size_t P, size_t Q>
Fixed<P,Q>& operator*=(Fixed<P,Q>& x, float a) {
    return x = x * Fixed<P,Q>(a);
}

template<size_t P, size_t Q>

Fixed<P,Q> operator-(Fixed<P,Q> x) { return Fixed<P,Q>::from_raw(-x.value); }


template<size_t P, size_t Q>
std::ostream &operator<<(std::ostream &out, Fixed<P,Q> x) {
    return out << static_cast<double>(x);
}

template<size_t P, size_t Q>
Fixed<P,Q> operator+(float a, Fixed<P,Q> x) {
    return Fixed<P,Q>(a) + x;
}

template<size_t P, size_t Q>
Fixed<P,Q> operator+(Fixed<P,Q> x, float a) {
    return x + Fixed<P,Q>(a);
}

template<size_t P, size_t Q>
Fixed<P,Q> operator-(float a, Fixed<P,Q> x) {
    return Fixed<P,Q>(a) - x;
}

template<size_t P, size_t Q>
Fixed<P,Q> operator-(Fixed<P,Q> x, float a) {
    return x - Fixed<P,Q>(a);
}

template<size_t P, size_t Q>
Fixed<P,Q> operator*(float a, Fixed<P,Q> x) {
    return Fixed<P,Q>(a) * x;
}

template<size_t P, size_t Q>
Fixed<P,Q> operator*(Fixed<P,Q> x, float a) {
    return x * Fixed<P,Q>(a);
}

template<size_t P, size_t Q>
Fixed<P,Q> operator/(float a, Fixed<P,Q> x) {
    return Fixed<P,Q>(a) / x;
}

template<size_t P, size_t Q>
Fixed<P,Q> operator/(Fixed<P,Q> x, float a) {
    return x / Fixed<P,Q>(a);
}

template<size_t P, size_t Q>
Fixed<P,Q>& operator+=(float& a, Fixed<P,Q> x) {
    a = static_cast<float>(Fixed<P,Q>(a) + x);
    return a;
}

template<size_t P, size_t Q>
float& operator-=(float& a, Fixed<P,Q> x) {
    a = static_cast<float>(Fixed<P,Q>(a) - x);
    return a;
}


