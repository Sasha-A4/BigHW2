#pragma once

#include <iostream>
#include <type_traits>
#include <cstdint>
#include <cmath>

template<size_t M, size_t L>
struct Fixed;

template<size_t M, size_t L>
struct FastFixed {
    static_assert(M > L, "M must be large than L");
    static constexpr size_t bits = M;
    static constexpr size_t frac = L;

    using StorageType = typename std::conditional<(M <= 8), int_fast8_t,
            typename std::conditional<(M <= 16), int_fast16_t,
                    typename std::conditional<(M <= 32), int_fast32_t, int_fast64_t>::type
            >::type
    >::type;

    StorageType value;

    constexpr FastFixed(int n = 0): value(static_cast<StorageType>(n) << L) {}
    constexpr FastFixed(float f): value(f * (StorageType(1) << L)) {}
    constexpr FastFixed(double f): value(f * (StorageType(1) << L)) {}

    static constexpr FastFixed from_raw(StorageType x) {
        FastFixed result;
        result.value = x;
        return result;
    }

    auto operator<=>(const FastFixed&) const = default;
    bool operator==(const FastFixed&) const = default;

    explicit operator float() const { return value / float(StorageType(1) << L); }
    explicit operator double() const { return value / double(StorageType(1) << L); }

    friend FastFixed operator/(FastFixed a, int b) {
        return FastFixed::from_raw(a.value / b);
    }

    friend FastFixed operator*(FastFixed a, int b) {
        return FastFixed::from_raw(a.value * b);
    }

    friend FastFixed operator*(int a, FastFixed b) {
        return b * a;
    }

    template<size_t M2, size_t L2>
    explicit operator FastFixed<M2,L2>() const {
        if constexpr (L2 >= L) {
            return FastFixed<M2,L2>::from_raw(static_cast<typename FastFixed<M2,L2>::StorageType>(value) << (L2 - L));
        } else {
            constexpr size_t shift = L - L2;
            if constexpr (shift >= M2) {
                auto temp = value >> (shift - M2 + 1);
                return FastFixed<M2,L2>::from_raw(static_cast<typename FastFixed<M2,L2>::StorageType>(temp) >> 1);
            } else {
                return FastFixed<M2,L2>::from_raw(static_cast<typename FastFixed<M2,L2>::StorageType>(value) >> shift);
            }
        }
    }

    template<size_t M2, size_t L2>
    explicit operator Fixed<M2,L2>() const {
        if constexpr (L2 >= L) {
            return Fixed<M2,L2>::from_raw(static_cast<typename Fixed<M2,L2>::StorageType>(value) << (L2 - L));
        } else {
            constexpr size_t shift = L - L2;
            if constexpr (shift >= M2) {
                auto temp = value >> (shift - M2 + 1);
                return Fixed<M2,L2>::from_raw(static_cast<typename Fixed<M2,L2>::StorageType>(temp) >> 1);
            } else {
                return Fixed<M2,L2>::from_raw(static_cast<typename Fixed<M2,L2>::StorageType>(value) >> shift);
            }
        }
    }
};

template<size_t M, size_t L>
FastFixed<M,L> operator+(FastFixed<M,L> x, FastFixed<M,L> y) {
    return FastFixed<M,L>::from_raw(x.value + y.value);
}

template<size_t M, size_t L>
FastFixed<M,L> operator-(FastFixed<M,L> x, FastFixed<M,L> y) {
    return FastFixed<M,L>::from_raw(x.value - y.value);
}

template<size_t M, size_t L>
FastFixed<M,L> &operator+=(FastFixed<M,L> &x, FastFixed<M,L> y) { return x = x + y; }

template<size_t M, size_t L>
FastFixed<M,L> &operator-=(FastFixed<M,L> &x, FastFixed<M,L> y) { return x = x - y; }

template<size_t M, size_t L>
FastFixed<M,L> &operator*=(FastFixed<M,L> &x, FastFixed<M,L> y) { return x = x * y; }

template<size_t M, size_t L>
FastFixed<M,L> &operator/=(FastFixed<M,L> &x, FastFixed<M,L> y) { return x = x / y; }

template<size_t M, size_t L>
FastFixed<M,L> operator-(FastFixed<M,L> x) { return FastFixed<M,L>::from_raw(-x.value); }

template<size_t M, size_t L>
std::ostream &operator<<(std::ostream &out, FastFixed<M,L> x) {
    return out << static_cast<double>(x);
}

template<size_t M, size_t L>
FastFixed<M,L> abs(FastFixed<M,L> x) {
    FastFixed<M,L> result = x;
    if (result.value < 0) result.value = -result.value;
    return result;
}


template<size_t M, size_t L>
FastFixed<M,L> operator*(FastFixed<M,L> x, FastFixed<M,L> y) {
    using ST = typename FastFixed<M,L>::StorageType;
    if constexpr (M <= 32) {
        return FastFixed<M,L>::from_raw((static_cast<int_fast64_t>(x.value) * y.value) >> L);
    } else {
        ST high = (x.value >> L) * y.value;
        ST low = (x.value & ((ST(1) << L) - 1)) * y.value >> L;
        return FastFixed<M,L>::from_raw(high + low);
    }
}

template<size_t M, size_t L>
FastFixed<M,L> operator/(FastFixed<M,L> x, FastFixed<M,L> y) {
    using ST = typename FastFixed<M,L>::StorageType;
    if constexpr (M <= 32) {
        return FastFixed<M,L>::from_raw((static_cast<int_fast64_t>(x.value) << L) / y.value);
    } else {
        return FastFixed<M,L>::from_raw((x.value << L) / y.value);
    }
}
