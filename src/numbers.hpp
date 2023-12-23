#pragma once

// Some of the less up-to-date compilers don't have `std::numbers` as of writing
// this, so this file serves as a replacement.
template <typename T>
T constexpr pi{static_cast<T>(3.14159265358979323846264338328)};

template <typename T>
T constexpr pi_halves{static_cast<T>(1.57079632679489661923132169164)};

template <typename T>
T constexpr three_pi_halves{static_cast<T>(4.71238898038468985769396507492)};

template <typename T>
T constexpr two_pi{static_cast<T>(6.28318530717958647692528676656)};

template <typename T>
T constexpr one_over_pi{static_cast<T>(0.318309886183790671537767526745)};

template <typename T>
T constexpr one_over_two_pi{static_cast<T>(0.159154943091895335768883763373)};

template <typename T>
T constexpr sqrt_two{static_cast<T>(1.41421356237309504880168872421)};

template <typename T>
T constexpr sqrt_one_half{static_cast<T>(.707106781186547524400844362105)};

template <typename T>
T constexpr three_hundred_and_sixty{static_cast<T>(360)};

template <typename T>
T constexpr two_pi_over_three_hundred_and_sixty{static_cast<T>(
  0.0174532925199432957692369076849)};

template <typename T>
T constexpr three_hundred_and_sixty_over_two_pi{static_cast<T>(
  57.2957795130823208767981548141)};

