#pragma once

#include <cmath>
#include <array>

// Some of the less up-to-date compilers don't have `std::numbers` as of writing
// this, so this file serves as a replacement.
template <typename T>
T constexpr pi{static_cast<T>(3.14159265358979323846264338328)};

template <typename T>
T constexpr pi_halves{static_cast<T>(1.57079632679489661923132169164)};

template <typename T>
T constexpr two_pi{static_cast<T>(6.28318530717958647692528676656)};

template <typename T>
T constexpr sqrt_two{static_cast<T>(1.41421356237309504880168872421)};

template <typename T>
T constexpr sqrt_one_half{static_cast<T>(.707106781186547524400844362105)};

template <typename T>
T constexpr three_hundred_and_sixty{static_cast<T>(360)};

template <typename T, bool inv = false>
std::array<T, 2> constexpr rotate(std::array<T, 2> const &v, T const &theta) {
  auto const &[x0, x1]{v};
  T const y0{std::cos(theta)}, y1{std::sin(theta)};
  if constexpr (not inv)
    return {x0 * y0 - x1 * y1, x0 * y1 + x1 * y0};
  else
    return {x0 * y0 + x1 * y1, x1 * y0 - x0 * y1};
}

// `std::cyl_bessel_j` does not exist on AppleClang as of writing this, so
// here's a polynomial approximation for the zeroth-order function.
template <typename T>
T bessel_i0(T x) {
   T abs_x = std::abs(x);
   if (abs_x < static_cast<T>(3.75)) {
      T const y = x * x / static_cast<T>(14.0625);
      return static_cast<T>(1.0) + y * (
        static_cast<T>(3.5156229) + y * (
          static_cast<T>(3.0899424) + y * (
            static_cast<T>(1.2067492) + y * (
              static_cast<T>(0.2659732) + y * (
                static_cast<T>(0.360768e-1) + y *
                  static_cast<T>(0.45813e-2))))));
   } else {
      T const y = static_cast<T>(3.75) / abs_x;
      return (exp(abs_x) / sqrt(abs_x)) * (
        static_cast<T>(0.39894228) + y * (
          static_cast<T>(0.1328592e-1) + y * (
            static_cast<T>(0.225319e-2) + y * (
              static_cast<T>(-0.157565e-2) + y * (
                static_cast<T>(0.916281e-2) + y * (
                  static_cast<T>(-0.2057706e-1) + y * (
                    static_cast<T>(0.2635537e-1) + y * (
                      static_cast<T>(-0.1647633e-1) + y *
                        static_cast<T>(0.392377e-2)))))))));
   }
}

