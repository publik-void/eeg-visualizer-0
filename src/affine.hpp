#pragma once

#include <cmath>
#include <array>

template <typename T, bool inv = false>
std::array<T, 2> constexpr rotate(std::array<T, 2> const &v, T const &theta) {
  auto const &[x0, x1]{v};
  T const y0{std::cos(theta)}, y1{std::sin(theta)};
  if constexpr (not inv)
    return {x0 * y0 - x1 * y1, x0 * y1 + x1 * y0};
  else
    return {x0 * y0 + x1 * y1, x1 * y0 - x0 * y1};
}


