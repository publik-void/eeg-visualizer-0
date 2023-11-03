#pragma once

#include <numbers>
#include <cmath>
#include <utility>
#include <array>
#include <algorithm>
#include <execution>
#include <iterator>

//#include <iostream>

namespace dsp {

template <std::size_t n, typename T>
std::array<T, n> fill(T const &value = static_cast<T>(0),
    std::array<T, n> a = {}) {
  std::fill(std::begin(a), std::end(a), value);
  return a;
}

enum struct FilterTopology {
  // …
  transposed_direct_form_2_b0a0_normalized
};

FilterTopology constexpr
  tdf2_b0a0{FilterTopology::transposed_direct_form_2_b0a0_normalized};

template <typename T>
struct FilterSpec {};

template <typename T, std::size_t ob = 0, std::size_t oa = 0,
  FilterTopology topology = tdf2_b0a0>
struct IIRSpec : public FilterSpec<T> {};

template <typename T, std::size_t ob, std::size_t oa>
struct IIRSpec<T, ob, oa, tdf2_b0a0> : public FilterSpec<T> {
  T static constexpr one{static_cast<T>(1)};
  std::size_t static constexpr nb{ob + 1}, na{oa + 1};

  std::array<T, ob> bs;
  std::array<T, oa> as;
  T scale;

  IIRSpec(std::array<T, nb> const &bs = fill<nb>(one),
      std::array<T, na> const &as = fill<na>(one), T const &scale = one) {
    T const b0{bs[0]}, a0{as[0]};
    T const b0_inv{one / b0}, a0_inv{one / a0};
    std::transform(std::next(std::cbegin(bs)), std::cend(bs),
      std::begin(this->bs), [&](T const x){ return x * b0_inv; });
    std::transform(std::next(std::cbegin(as)), std::cend(as),
      std::begin(this->as), [&](T const x){ return x * -a0_inv; });
    this->scale = scale * b0 * a0_inv;
  }

  template <typename TB, typename TA>
  IIRSpec(std::array<TB, nb> const &bs,
      std::array<TA, na> const &as = fill<na>(static_cast<TA>(1)),
      T const &scale = one) : IIRSpec(
        [&](){std::array<T, nb> _bs;
              std::copy(std::cbegin(bs), std::cend(bs), std::begin(_bs));
              return _bs; }(),
        [&](){std::array<T, na> _as;
              std::copy(std::cbegin(as), std::cend(as), std::begin(_as));
              return _as; }(),
          scale) {}
};

using BiquadDFT2SingleSpec = IIRSpec< float, 2, 2, tdf2_b0a0>;
using BiquadDFT2DoubleSpec = IIRSpec<double, 2, 2, tdf2_b0a0>;

template <typename T>
struct FilterState {};

template <typename T, std::size_t ob = 0, std::size_t oa = 0,
  FilterTopology topology = tdf2_b0a0>
struct IIRState : public FilterState<T> {};

template <typename T, std::size_t ob, std::size_t oa>
struct IIRState<T, ob, oa, tdf2_b0a0> : public FilterState<T> {
  std::size_t static constexpr nb{ob + 1}, na{oa + 1};
  std::size_t static constexpr ns{std::max(nb, na)};
  std::size_t static constexpr n_channels_default{1};

  std::vector<std::array<T, ns>> sss;

  IIRState(std::size_t const n_channels = n_channels_default) :
      sss{std::vector<std::array<T, ns>>(n_channels, fill<ns, T>())} {}
};

template <typename T, std::size_t ob, std::size_t oa>
inline std::vector<T *> firsts(IIRState<T, ob, oa, tdf2_b0a0> &state) {
  std::vector<T *> s_firsts(state.sss.size());
  std::transform(std::begin(state.sss), std::end(state.sss),
    std::begin(s_firsts), [&](auto &ss){ return ss.data(); });
  return s_firsts;
}

using BiquadDFT2SingleState = IIRState< float, 2, 2, tdf2_b0a0>;
using BiquadDFT2DoubleState = IIRState<double, 2, 2, tdf2_b0a0>;

template <typename T, std::size_t ob, std::size_t oa, FilterTopology topology>
IIRState<T, ob, oa, topology> make_state(
    IIRSpec<T, ob, oa, topology> const &, std::size_t const n_channels =
      IIRState<T, ob, oa, topology>::n_channels_default) {
  return IIRState<T, ob, oa, topology>(n_channels);
}

template <class XItr, class YItr, class CItr, class SItrItr,
  FilterTopology topology = tdf2_b0a0>
inline void tick(XItr const &x_first, XItr const &x_last, YItr const &y_first,
    CItr const &b_first, CItr const &b_last,
    CItr const &a_first, CItr const &a_last,
    SItrItr const &s_first_first, std::size_t const ns) {
  if constexpr (topology == tdf2_b0a0) {
    bool const ss_empty{b_first >= b_last and a_first >= a_last};

    std::transform(//std::execution::unseq,
      x_first, x_last, s_first_first, y_first,
      [&](auto const &x, auto const &s_first){
        auto const y{ss_empty ? x : x + (*s_first)};

        auto s_itr{s_first};
        for (size_t i{0}; i < ns - 1; ++i) {
          auto const s{s_itr};
          *s = *(++s_itr);
        }
        *s_itr = static_cast<decltype(y)>(0);

        std::transform(//std::execution::unseq,
          b_first, b_last, s_first, s_first,
          [&](auto const &b, auto const &s){ return s + b * x; });
        std::transform(//std::execution::unseq,
          a_first, a_last, s_first, s_first,
          [&](auto const &a, auto const &s){ return s + a * y; });

        // TODO: subnormal canceling
        return y;
      });
  }
}

template <typename T, class XItr, class YItr, std::size_t ob, std::size_t oa,
  FilterTopology topology>
inline void tick(XItr const &x_first, XItr const &x_last, YItr const &y_first,
    IIRSpec<T, ob, oa, topology> const &spec,
    IIRState<T, ob, oa, topology> &state) {
  auto const s_firsts(firsts(state));

  return tick(x_first, x_last, y_first,
    std::cbegin(spec.bs), std::cend(spec.bs),
    std::cbegin(spec.as), std::cend(spec.as),
    std::cbegin(s_firsts), state.ns);
}

enum struct BiquadShape {
  zero,
  identity,
  lowpass,
  highpass
  // …
};

// Robert Bristow-Johnson biquads
// `https://github.com/WebAudio/Audio-EQ-Cookbook.git`

template <typename T, BiquadShape shape>
inline auto biquad_rbj(
    T const &omega_0 = std::numbers::pi_v<T> / static_cast<T>(2),
    T const &q = static_cast<T>(1) / std::numbers::sqrt2_v<T>,
    [[maybe_unused]] /*TODO*/ T const &a = static_cast<T>(1)) {
  T const zero{static_cast<T>(0)};
  T const one_half{static_cast<T>(.5)};
  T const one{static_cast<T>(1)};
  T const two{static_cast<T>(2)};
  T const sin_omega_0{std::sin(omega_0)};
  T const cos_omega_0{std::cos(omega_0)};
  T const alpha{sin_omega_0 / (two * q)};

  std::array<T, 3> bs, as;
  if constexpr (shape == BiquadShape::zero) {
    bs[0] = zero; bs[1] = zero; bs[2] = zero;
    as[0] = one; as[1] = zero; as[2] = zero;
  } else if constexpr (shape == BiquadShape::identity) {
    bs[0] = one; bs[1] = zero; bs[2] = zero;
    as[0] = one; as[1] = zero; as[2] = zero;
  } else if constexpr (shape == BiquadShape::lowpass) {
    bs[1] = one - cos_omega_0;
    bs[0] = bs[1] * one_half;
    bs[2] = bs[0];
    as[0] = one + alpha;
    as[1] = -two * cos_omega_0;
    as[2] = one - alpha;
  } else if constexpr (shape == BiquadShape::highpass) {
    T const one_plus_cos_omega_0{one + cos_omega_0};
    bs[0] = one_plus_cos_omega_0 * one_half;
    bs[1] = -one_plus_cos_omega_0;
    bs[2] = bs[0];
    as[0] = one + alpha;
    as[1] = -two * cos_omega_0;
    as[2] = one - alpha;
  }
  // …
  return std::make_pair(bs, as);
}

enum struct WindowShape {
  identity,
  blackman,
  kaiser
  // …
};

template<typename T, WindowShape shape>
inline T window(std::size_t const i, std::size_t const n, T const parameter_0) {
  T const zero{static_cast<T>(0)};
  T const one_half{static_cast<T>(.5)};
  T const one{static_cast<T>(1)};
  T const two{static_cast<T>(2)};
  if constexpr (shape == WindowShape::identity) {
    return one;
  } else if constexpr (shape == WindowShape::blackman) {
    T const alpha{parameter_0};
    T const alpha_2{alpha * one_half};
    T const alpha_0{one_half - alpha_2};
    T const alpha_1{one_half};
    T const two_pi_i_over_n_minus_one{
      two * std::numbers::pi_v<T> * i / (n - 1)};
    return alpha_0 - alpha_1 * std::cos(two_pi_i_over_n_minus_one) +
      alpha_2 * std::cos(two * two_pi_i_over_n_minus_one);
  } else if constexpr (shape == WindowShape::kaiser) {
    //auto const square{[](T const &x){ return x * x; }};
    //auto const i0{[](T const &x){
    //  return std::cyl_bessel_j(zero, x); }};
    //T const alpha{parameter_0};
    //T const pi_alpha{std::numbers::pi_v<T> * alpha};
    //return i0(pi_alpha * std::sqrt(one - square((two * i) / (n - 1) - one))) /
    //  i0(pi_alpha);
    // NOTE: `std::cyl_bessel_j` does not exist on AppleClang 14?
  }
}

} // namespace dsp
