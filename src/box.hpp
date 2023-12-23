#pragma once

#include <array>

#include "lerp.hpp"
#include "affine.hpp"

template <typename T>
struct Box {
  T x0, y0, x1, y1;

  std::array<T, 2> p0() const { return {this->x0, this->y0}; }
  std::array<T, 2> p1() const { return {this->x1, this->y1}; }
  std::array<T, 2> xs() const { return {this->x0, this->x1}; }
  std::array<T, 2> ys() const { return {this->y0, this->y1}; }
  T xd() const { return this->x1 - this->x0; }
  T yd() const { return this->y1 - this->y0; }
  std::array<T, 2> center() const { return lerp(
    std::array<T, 2>{static_cast<T>(.5), static_cast<T>(.5)}, *this); }
  std::array<T, 2> p00() const { return {this->x0, this->y0}; }
  std::array<T, 2> p01() const { return {this->x1, this->y0}; }
  std::array<T, 2> p10() const { return {this->x0, this->y1}; }
  std::array<T, 2> p11() const { return {this->x1, this->y1}; }
  std::array<std::array<T, 2>, 4> vertices() const {
    return {{this->p00(), this->p01(), this->p10(), this->p11()}};
  }
  std::array<std::array<std::array<T, 2>, 2>, 4> edges(
      std::array<T, 2> const &r = {static_cast<T>(0), static_cast<T>(0)}
    ) const {
    return {{{{{{this->x0 - std::get<0>(r), this->y0}},
               {{this->x1 - std::get<0>(r), this->y0}}}},
             {{{{this->x0, this->y0 + std::get<1>(r)}},
               {{this->x0, this->y1 + std::get<1>(r)}}}},
             {{{{this->x0 + std::get<0>(r), this->y1}},
               {{this->x1 + std::get<0>(r), this->y1}}}},
             {{{{this->x1, this->y0 - std::get<1>(r)}},
               {{this->x1, this->y1 - std::get<1>(r)}}}}}};
  }

  constexpr Box() {}
  constexpr Box(T const &x0, T const &y0, T const &x1, T const &y1) :
      x0{x0}, y0{y0}, x1{x1}, y1{y1} {}
  constexpr Box(std::array<T, 2> const &p0, std::array<T, 2> const &p1) :
      x0{std::get<0>(p0)}, y0{std::get<1>(p0)},
      x1{std::get<0>(p1)}, y1{std::get<1>(p1)} {}
};

template <typename T>
Box<T> operator+(Box<T> const &self, std::array<T, 2> const &p) {
  return Box{self.x0 + std::get<0>(p), self.y0 + std::get<1>(p),
             self.x1 + std::get<0>(p), self.y1 + std::get<1>(p)};
}

template <typename T>
Box<T> operator-(Box<T> const &self, std::array<T, 2> const &p) {
  return Box{self.x0 - std::get<0>(p), self.y0 - std::get<1>(p),
             self.x1 - std::get<0>(p), self.y1 - std::get<1>(p)};
}

template <typename T>
Box<T> operator*(Box<T> const &self, std::array<T, 2> const &p) {
  return Box{self.x0, self.y0, lerp(std::get<0>(p), {self.x0, self.x1}),
                               lerp(std::get<1>(p), {self.y0, self.y1})};
}

template <typename T>
Box<T> &operator+=(Box<T> &self, auto const &x) { return self = self + x; }

template <typename T>
Box<T> &operator-=(Box<T> &self, auto const &x) { return self = self - x; }

template <typename T>
Box<T> &operator*=(Box<T> &self, auto const &x) { return self = self * x; }

template <typename T>
std::array<T, 2> offset(std::array<T, 2> const &p, Box<T> const &box) {
  auto const &[p0, p1]{p};
  return {{p0 + box.x0, p1 + box.y0}};
}

template <typename T>
std::array<T, 2> scale(std::array<T, 2> const &p, Box<T> const &box) {
  auto const &[p0, p1]{p};
  return {{p0 * box.xd(), p1 * box.yd()}};
}

template <bool inverse, typename T>
std::array<T, 2> _lerp(std::array<T, 2> const &p, Box<T> const &box,
    std::array<T, 2> &&result = {}) {
  auto const &[p0, p1]{p};
  return result = {{_lerp<inverse>(p0, box.xs()),
                    _lerp<inverse>(p1, box.ys())}};
}

template <bool inverse, typename T>
array2d<T, 2, 1> _lerp(array2d<T, 2, 1> const &p, Box<T> const &box,
    array2d<T, 2, 1> &&result = {}) {
  return result = array_to_array2d(_lerp<inverse>(array2d_to_array(p), box));
}

template <bool inverse, typename T>
Box<T> _lerp(Box<T> const &ps, Box<T> const &box, Box<T> &&result = {}) {
  return result = {_lerp<inverse>(ps.p0(), box), _lerp<inverse>(ps.p1(), box)};
}

template <typename T>
Box<T> box_to_box(Box<T> const &from, Box<T> const &to) {
  Lerp<T> const from_to_x{Lerp<T>(to.xs()) * iLerp<T>(from.xs())};
  Lerp<T> const from_to_y{Lerp<T>(to.ys()) * iLerp<T>(from.ys())};
  return {from_to_x.y0, from_to_y.y0, from_to_x.y1(), from_to_y.y1()};
}

