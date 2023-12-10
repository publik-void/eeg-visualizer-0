// Stack-allocateable, inlineable basic 2D linear algebra

#pragma once

#include <cstddef>
#include <cmath>
#include <optional>
#include <array>

#include "SFML/System.hpp"

// `constexpr` (unrolled) for-loop
template <auto n, auto i = std::size_t{0}, auto inc = std::size_t{1}, class F>
void constexpr for_constexpr(F &&f) {
  if constexpr (i < n) {
    f(std::integral_constant<decltype(i), i>{});
    for_constexpr<n, i + inc, inc>(f);
  }
}

template <auto n, auto m, auto i = std::size_t{0}, auto i_inc = std::size_t{1},
    auto j = std::size_t{0}, auto j_inc = std::size_t{1}, class F>
void constexpr for2d_constexpr(F &&f) {
  for_constexpr<n, i, i_inc>([&f](auto const i_idx){
      for_constexpr<m, j, j_inc>([&f, &i_idx](auto const j_idx){
          f(i_idx, j_idx);
        });
    });
}

//template <class F, class T, class... Ts>
//void constexpr for_constexpr(F &&f, T &&t, Ts &&...ts) {
//  std::size_t constexpr n{std::tuple_size_v<std::decay_t<T>>};
//  for_constexpr<n>([&](auto const i){
//    f(std::get<i.value>(t), std::get<i.value>(ts)...);
//  });
//}

template<typename T>
std::array<T, 2> sf_vec_to_array(sf::Vector2<T> const &v) {
  return {v.x, v.y};
}

template<typename T>
sf::Vector2<T> array_to_sf_vec(std::array<T, 2> const &v) {
  return {std::get<0>(v), std::get<1>(v)};
}

template <typename T>
T abs2(T const &x) {
  // TODO: Does the compiler optimize out the `std::abs` for scalar types?
  return std::abs(x * x);
}

template <typename T, std::size_t n>
T abs2(std::array<T, n> const &a) {
  T y{static_cast<T>(0)};
  for_constexpr<n>([&](auto const i){
      T const &x{std::get<i.value>(a)};
      y += abs2(x);
    });
  return y;
}

// TODO: This happens to be row-major order. Is that a good/bad choice?
template <typename T, std::size_t n, std::size_t m>
using array2d = std::array<std::array<T, m>, n>;

template <std::size_t i, std::size_t j, class T, std::size_t n, std::size_t m>
T constexpr &get(array2d<T, n, m> &a) noexcept {
  return std::get<j>(std::get<i>(a));
}

template <std::size_t i, std::size_t j, class T, std::size_t n, std::size_t m>
T constexpr &&get(array2d<T, n, m> &&a) noexcept {
  return std::get<j>(std::get<i>(a));
}

template <std::size_t i, std::size_t j, class T, std::size_t n, std::size_t m>
T constexpr const &get(array2d<T, n, m> const &a) noexcept {
  return std::get<j>(std::get<i>(a));
}

template <std::size_t i, std::size_t j, class T, std::size_t n, std::size_t m>
T constexpr const &&get(array2d<T, n, m> const &&a) noexcept {
  return std::get<j>(std::get<i>(a));
}

template <typename T, std::size_t n>
std::array<T, n> constexpr copy(std::array<T, n> const a) {
  return a;
}

template <typename T, std::size_t n, std::size_t m>
array2d<T, n, m> constexpr fill(T const &c = static_cast<T>(0),
    std::optional<T> const &diag = {}, array2d<T, n, m> &&result = {}) {
  for2d_constexpr<n, m>([&](auto const i, auto const j){
      auto const _i{i.value}, _j{j.value};
      get<_i, _j>(result) = (_i == _j) ? diag.value_or(c) : c;
    });
  return result;
}

template <typename T, std::size_t n, std::size_t m>
array2d<T, m, n> constexpr transpose(array2d<T, n, m> const &a,
    array2d<T, m, n> &&result = {}) {
  for2d_constexpr<n, m>([&](auto const i, auto const j){
      auto const _i{i.value}, _j{j.value};
      get<_j, _i>(result) = get<_i, _j>(a);
    });
  return result;
}

template <typename T, std::size_t n>
array2d<T, n, 1> constexpr array_to_array2d(std::array<T, n> const &v) {
  return transpose(std::array<std::array<T, n>, 1>{v});
}

template <typename T, std::size_t n>
std::array<T, n> constexpr array2d_to_array(array2d<T, n, 1> const &a) {
  return std::get<0>(transpose(a));
}

template<typename T>
array2d<T, 2, 1> sf_vec_to_array2d(sf::Vector2<T> const &v) {
  return array_to_array2d(sf_vec_to_array(v));
}

template<typename T>
sf::Vector2<T> array2d_to_sf_vec(array2d<T, 2, 1> const &v) {
  return array_to_sf_vec(array2d_to_array(v));
}

template <typename T>
T constexpr neg(T const &a, T &&result = {}) {
  return result = -a;
}

template <typename T, std::size_t n>
std::array<T, n> constexpr neg(std::array<T, n> const &a,
    std::array<T, n> &&result = {}) {
  for_constexpr<n>([&](auto const i){
      auto const _i{i.value};
      get<_i>(result) = neg(get<_i>(a));
    });
  return result;
}

template <typename T>
T constexpr add(T const &a, T const &b, T &&result = {}) {
  return result = a + b;
}

template <typename T, std::size_t n>
std::array<T, n> constexpr add(
    std::array<T, n> const &a, std::array<T, n> const &b,
    std::array<T, n> &&result = {}) {
  for_constexpr<n>([&](auto const i){
      auto const _i{i.value};
      std::get<_i>(result) =
        add(std::get<_i>(a), std::get<_i>(b), std::move(std::get<_i>(result)));
    });
  return result;
}

template <typename T, std::size_t n, std::size_t m, std::size_t l>
array2d<T, n, l> constexpr muladd(
    array2d<T, n, m> const &a, array2d<T, m, l> const &b,
    array2d<T, n, l> &&c = fill<T, n, l>()) {
  // This will mainly be used for small matrices, so we'll keep it simple.
  for2d_constexpr<n, l>([&](auto const i, auto const k){
      auto const _i{i.value}, _k{k.value};
      for_constexpr<m>([&](auto const j){
        auto const _j{j.value};
          get<_i, _k>(c) += get<_i, _j>(a) * get<_j, _k>(b);
        });
    });
  return c;
}

template <typename T, std::size_t n, std::size_t m>
struct Affine {
  array2d<T, n, m> a;
  array2d<T, n, 1> b;

  array2d<T, n, 1> constexpr operator()(array2d<T, m, 1> const &x) const {
    return muladd(this->a, x, copy(this->b));
  }

  std::array<T, n> constexpr operator()(std::array<T, n> const &x) const {
    return array2d_to_array((*this)(array_to_array2d(x)));
  }

  sf::Vector2f operator()(sf::Vector2f const &x) const {
    return array_to_sf_vec((*this)(sf_vec_to_array(x)));
  }

  template <typename U>
  std::optional<U> constexpr operator()(std::optional<U> const &x_opt) const {
    if (x_opt.has_value()) return (*this)(x_opt.value());
    else return {};
  }

  template <std::size_t l>
  Affine<T, n, l> constexpr operator()(Affine<T, m, l> const &that) const {
    return {muladd(this->a, that.a), muladd(this->a, that.b, copy(this->b))};
  }

  static Affine<T, n, m> constexpr identity() {
    return {fill<T, n, m>(static_cast<T>(0), static_cast<T>(1)),
      fill<T, n, 1>()};
  }
};

//template <typename T, bool inv = false>
//std::array<T, 2> constexpr rotate(std::array<T, 2> const &v, T const &theta) {
//  auto const &[x0, x1]{v};
//  T const y0{std::cos(theta)}, y1{std::sin(theta)};
//  if constexpr (not inv)
//    return {x0 * y0 - x1 * y1, x0 * y1 + x1 * y0};
//  else
//    return {x0 * y0 + x1 * y1, x1 * y0 - x0 * y1};
//}

//template <bool inv = false>
//sf::Vector2f rotate(sf::Vector2f const &v, float const theta) {
//  return array_to_sf_vec(rotate<float, inv>(sf_vec_to_array(v), theta));
//}

