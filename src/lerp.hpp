#pragma once

#include <cstddef>
#include <array>

template <bool inverse, typename S, typename T>
S _lerp(S const &x, std::array<T, 2> const &ys, S &&result = {}) {
  auto const &[y0, y1]{ys};
  if constexpr (not inverse) return result = (y1 - y0) * x + y0;
  else return result = (x - y0) / (y1 - y0);
}

template <bool inverse, typename S, typename T>
S _lerp(S const &x, T const &xs, T const &ys,
    S &&result = {}) {
  return _lerp<inverse>(_lerp<not inverse>(x,
      inverse ? ys : xs, std::move(result)),
    inverse ? xs : ys, std::move(result));
}

template <bool inverse, typename S, std::size_t n, typename T>
std::array<S, n> _lerp(std::array<S, n> const &x, T const &ys,
    std::array<S, n> &&result = {}) {
  for (std::size_t i{0}; i < n; ++i)
    _lerp<inverse>(x[i], ys, std::move(result[i]));
  return result;
}

// The idea with the above templating was to make some of this code less verbose
// and redundant. Unfortunately, one consequence of this is that several
// overloads with`std::initializer_list` need to be defined in the below,
// because the type of the latter is not deducible by templates, for some
// reason. Well, at least the above templating efforts make it simpler to extend
// this lerping function.
namespace {
  template <typename T, std::size_t n = 2>
  std::array<T, n> constexpr il_to_ar(std::initializer_list<T> const &il,
      std::array<T, n> &&ar = {}) {
    assert(il.size() == n);
    auto itr_il{il.begin()}; auto itr_ar{ar.begin()};
    while (/*itr_il < il.end() and */itr_ar < ar.end()) {
      *itr_ar = *itr_il;
      ++itr_il; ++itr_ar;
    }
    return ar;
  }
}

template <typename... Params, typename S, typename T>
S lerp(S const &x, T const &ys, S &&result = {}) {
  return _lerp<false, Params...>(x, ys, std::move(result));
}

template <typename... Params, typename S, typename T>
S ilerp(S const &x, T const &ys, S &&result = {}) {
  return _lerp<true, Params...>(x, ys, std::move(result));
}

template <typename... Params, typename S, typename T>
S lerp(S const &x, std::initializer_list<T> const &&ys, S &&result = {}) {
  return _lerp<false, Params...>(x, il_to_ar(ys), std::move(result));
}

template <typename... Params, typename S, typename T>
S ilerp(S const &x, std::initializer_list<T> const &&ys, S &&result = {}) {
  return _lerp<true, Params...>(x, il_to_ar(ys), std::move(result));
}

template <typename... Params, typename S, typename T>
S lerp(S const &x, T const &xs, T const &ys, S &&result = {}) {
  return _lerp<false, Params...>(x, xs, ys, std::move(result));
}

template <typename... Params, typename S, typename T>
S ilerp(S const &x, T const &xs, T const &ys, S &&result = {}) {
  return _lerp<true, Params...>(x, xs, ys, std::move(result));
}

template <typename... Params, typename S, typename T, typename U>
S lerp(S const &x, std::initializer_list<T> const &&xs,
    U const &ys, S &&result = {}) {
  return _lerp<false, Params...>(x, il_to_ar(xs), ys, std::move(result));
}

template <typename... Params, typename S, typename T, typename U>
S ilerp(S const &x, std::initializer_list<T> const &&xs,
    U const &ys, S &&result = {}) {
  return _lerp<true, Params...>(x, il_to_ar(xs), ys, std::move(result));
}

template <typename... Params, typename S, typename T, typename U>
S lerp(S const &x, U const &xs,
    std::initializer_list<T> const &&ys, S &&result = {}) {
  return _lerp<false, Params...>(x, xs, il_to_ar(ys), std::move(result));
}

template <typename... Params, typename S, typename T, typename U>
S ilerp(S const &x, U const &xs,
    std::initializer_list<T> const &&ys, S &&result = {}) {
  return _lerp<true, Params...>(x, xs, il_to_ar(ys), std::move(result));
}

template <typename... Params, typename S, typename T>
S lerp(S const &x, std::initializer_list<T> const &&xs,
    std::initializer_list<T> const &&ys, S &&result = {}) {
  return _lerp<false, Params...>(x, il_to_ar(xs), il_to_ar(ys),
    std::move(result));
}

template <typename... Params, typename S, typename T>
S ilerp(S const &x, std::initializer_list<T> const &&xs,
    std::initializer_list<T> const &&ys, S &&result = {}) {
  return _lerp<true, Params...>(x, il_to_ar(xs), il_to_ar(ys),
    std::move(result));
}

template <typename T>
struct Lerp {
  T y0, yd;

  T y1() const { return yd + y0; }
  std::array<T, 2> ys() const { return {{y0, this->y1()}}; }

  template <typename S>
  S operator()(S const &x, S &&result = {}) const {
    return result = x * yd + y0;
  }

  template <typename S, std::size_t n>
  std::array<S, n> operator()(std::array<S, n> const &x,
      std::array<S, n> &&result = {}) const {
    for (std::size_t i{0}; i < n; ++i) result[i] = (*this)(x[i]);
    return result;
  }

  Lerp<T> operator*(Lerp<T> const &that) const {
    Lerp<T> result;
    result.y0 = that.y0 * this->yd + this->y0;
    result.yd = that.yd * this->yd;
    return result;
  }

  Lerp() {}
  Lerp(T const &y0, T const &y1) : y0{y0}, yd{y1 - y0} {}
  Lerp(std::array<T, 2> const &ys) : Lerp(std::get<0>(ys), std::get<1>(ys)) {}
};

template <typename T>
Lerp<T> iLerp(T const &y0, T const &y1, Lerp<T> &&result = {}) {
  T const neg_yd{y0 - y1};
  result.y0 = y0 / neg_yd;
  result.yd = static_cast<T>(-1) / neg_yd;
  return result;
}

template <typename T>
Lerp<T> iLerp(std::array<T, 2> const &ys) {
  return iLerp(std::get<0>(ys), std::get<1>(ys));
}

