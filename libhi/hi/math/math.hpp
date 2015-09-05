#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif _USE_MATH_DEFINES
#include <cmath>
#undef _USE_MATH_DEFINES

#include <hi/math/promote_traits.hpp>

#include <mmintrin.h>   // MML
#include <xmmintrin.h>  // SSE
#include <emmintrin.h>  // SSE2
//#include <pmmintrin.h> // SSE3
//#include <tmmintrin.h> // SSE4

///
/// 一般的なユーティリティ
///
namespace hi {
// x^2
template <typename T>
inline T square_of(T const &x) {
  return x * x;
}

// x^3
template <typename T>
inline T cube_of(T const &x) {
  return x * hi::square_of(x);
}

// x^4
template <typename T>
inline T fourth_power_of(T const &x) {
  return hi::square_of(hi::square_of(x));
}

// x^5
template <typename T>
inline T fifth_power_of(T const &x) {
  return x * fourth_power_of(x);
}

// x^n
template <typename T>
inline T pow(T const &x, std::size_t const &n) {
  T y = T(1);
  std::size_t i = (n + (8 - 1)) / 8;
  switch (n & (8 - 1))  // Duff's Device
  {
    do {
    case 0:
      y *= x;
    case 7:
      y *= x;
    case 6:
      y *= x;
    case 5:
      y *= x;
    case 4:
      y *= x;
    case 3:
      y *= x;
    case 2:
      y *= x;
    case 1:
      y *= x;
    } while (--i > 0);
  }
  return y;
}

/// 階乗を計算する．
template <typename T>
T factorial(T const &n) {
  T value = 1;
  for (T i = 0; i < n; ++i) {
    value *= (n - i);
  }
  return value;
}

/// 離散的累乗を計算する．
/// x^_n = x * (x-1) * (x-2) * ... * (x-(n-1))
template <typename T>
T discretePow(T const &x, T const &n) {
  T value = 1;
  for (T i = 0; i < n; ++i) {
    value *= (x - i);
  }
  return value;
}

/// 二項係数を計算する．
template <typename T>
T binomialCoefficient(T const &n, T const &r) {
  return hi::discretePow(n, r) / hi::factorial(r);
}

// min{x, y}
template <typename X, typename Y>
inline typename promote_traits<X, Y>::type min(X const &x, Y const &y) {
  return (x < y) ? x : y;
}

// min{x, y, z}
template <typename X, typename Y, typename Z>
inline typename promote_traits<typename promote_traits<X, Y>::type, Z>::type
min(X const &x, Y const &y, Z const &z) {
  return hi::min(hi::min(x, y), z);
}

// min{x, y, z, w}
template <typename X, typename Y, typename Z, typename W>
inline typename promote_traits<typename promote_traits<X, Y>::type,
                               typename promote_traits<Z, W>::type>::type
min(X const &x, Y const &y, Z const &z, W const &w) {
  return hi::min(hi::min(x, y), hi::min(z, w));
}

// max{x, y}
template <typename X, typename Y>
inline typename promote_traits<X, Y>::type max(X const &x, Y const &y) {
  return (x < y) ? y : x;
}

// max{x, y, z}
template <typename X, typename Y, typename Z>
inline typename promote_traits<typename promote_traits<X, Y>::type, Z>::type
max(X const &x, Y const &y, Z const &z) {
  return hi::max(hi::max(x, y), z);
}

// max{x, y, z, w}
template <typename X, typename Y, typename Z, typename W>
inline typename promote_traits<typename promote_traits<X, Y>::type,
                               typename promote_traits<Z, W>::type>::type
max(X const &x, Y const &y, Z const &z, W const &w) {
  return hi::max(hi::max(x, y), hi::max(z, w));
}

/// <summary>
/// 逆数(ReCiProcal)を計算する。
/// </summary>
template <typename T>
inline T rcp(T const &x) {
  return T(1) / x;
}

/// <summary>
/// 逆平方根(Reciprocal SQuare RooT)を計算する。
/// </summary>
template <typename T>
inline T rsqrt(T const &x) {
  return T(1) / std::sqrt(x);
}

template <typename T>
inline T lerp(T const &a, T const &b, T const &t) {
  return a + (b - a) * t;
}

template <typename T>
inline T to_degrees(T const &radians) {
  return T(180) * radians / T(M_PI);
}

template <typename T>
inline T to_radians(T const &degrees) {
  return T(M_PI) * degrees / T(180);
}

}  // end of namespace hi

namespace hi {
// x = min{x, y}
template <typename X, typename Y>
inline void min(X &x, Y const &y) {
  if (x > y) {
    x = static_cast<X>(y);
  }
}

// x = max{x, y}
template <typename X, typename Y>
inline void max(X &x, Y const &y) {
  if (x < y) {
    x = static_cast<X>(y);
  }
}
}  // end of namespace hi

#ifdef _USE_SSE_OPTIMIZATION
namespace hi {
template <>
inline float min(float const &a, float const &b) {
  float c;
  _mm_store_ss(&c, _mm_min_ss(_mm_load_ss(&a), _mm_load_ss(&b)));
  return c;
}

template <>
inline float max(float const &a, float const &b) {
  float c;
  _mm_store_ss(&c, _mm_max_ss(_mm_load_ss(&a), _mm_load_ss(&b)));
  return c;
}
}  // namespace hi
#endif _USE_SSE_OPTIMIZATION

#ifdef _USE_SSE2_OPTIMIZATION
namespace hi {
template <>
inline double min(double const &a, double const &b) {
  double c;
  _mm_store_sd(&c, _mm_min_sd(_mm_load_sd(&a), _mm_load_sd(&b)));
  return c;
}

template <>
inline double max(double const &a, double const &b) {
  double c;
  _mm_store_ss(&c, _mm_max_ss(_mm_load_sd(&a), _mm_load_sd(&b)));
  return c;
}
}  // namespace hi
#endif _USE_SSE2_OPTIMIZATION

/// <summary>
/// ベクトル演算用のユーティリティ
/// </summary>
namespace hi {
template <typename T>
inline typename T::value_type length_squared(T const &a) {
  return hi::dot(a, a);
}

template <typename T>
inline typename T::value_type length(T const &a) {
  return std::sqrt(hi::length_squared(a));
}

template <typename T>
inline typename T::value_type distance_squared(T const &a, T const &b) {
  return hi::length_squared(a - b);
}

template <typename T>
inline typename T::value_type distance(T const &a, T const &b) {
  return std::sqrt(hi::distance_squared(a, b));
}

template <typename T>
inline T normalize(T const &a) {
  return a * hi::rsqrt(hi::length_squared(a));
}

template <typename T>
inline typename T::value_type average(T const &a) {
  return hi::sum(a) * hi::rcp<T::value_type>(a.size());
}

}  // end of namespace hi

namespace hi {
template <typename T>
inline T divexplog(T const &a, T const &b) {
  return std::exp(std::log(a) - std::log(b));
}
/*
  // boundary test
  std::tcerr << _TEXT(" 0/ 0 = ") << hi::divexplog( 0.0,  0.0) << std::endl;
  std::tcerr << _TEXT(" 0/ 1 = ") << hi::divexplog( 0.0,  1.0) << std::endl;
  std::tcerr << _TEXT(" 1/ 0 = ") << hi::divexplog( 1.0,  0.0) << std::endl;
  std::tcerr << _TEXT(" 1/ 1 = ") << hi::divexplog( 1.0,  1.0) << std::endl;
  std::tcerr << _TEXT(" 0/-0 = ") << hi::divexplog( 0.0, -0.0) << std::endl;
  std::tcerr << _TEXT(" 0/-1 = ") << hi::divexplog( 0.0, -1.0) << std::endl;
  std::tcerr << _TEXT(" 1/-0 = ") << hi::divexplog( 1.0, -0.0) << std::endl;
  std::tcerr << _TEXT(" 1/-1 = ") << hi::divexplog( 1.0, -1.0) << std::endl;
  std::tcerr << _TEXT("-0/ 0 = ") << hi::divexplog(-0.0,  0.0) << std::endl;
  std::tcerr << _TEXT("-0/ 1 = ") << hi::divexplog(-0.0,  1.0) << std::endl;
  std::tcerr << _TEXT("-1/ 0 = ") << hi::divexplog(-1.0,  0.0) << std::endl;
  std::tcerr << _TEXT("-1/ 1 = ") << hi::divexplog(-1.0,  1.0) << std::endl;
  std::tcerr << _TEXT("-0/-0 = ") << hi::divexplog(-0.0, -0.0) << std::endl;
  std::tcerr << _TEXT("-0/-1 = ") << hi::divexplog(-0.0, -1.0) << std::endl;
  std::tcerr << _TEXT("-1/-0 = ") << hi::divexplog(-1.0, -0.0) << std::endl;
  std::tcerr << _TEXT("-1/-1 = ") << hi::divexplog(-1.0, -1.0) << std::endl;
*/

// 3 bit 間隔で値がある 30bit データから 10 bit データを取り出す
// in : XXx XXx XXx XXx XXx XXx XXx XXx XXx XXx (8)
// out: 000 000 000 000 000 000 00x xxx xxx xxx (8)
inline int pw_compact3(int x) {  // 1 111 111 111
  x = (x & 01111111111) | ((x >> 2) & 02222222222) |
      ((x >> 4) & 04444444444);                       // 7 007 007 007
  x = (x & 00007000007) | ((x >> 6) & 00070000070);   // 0 077 000 077
  x = (x & 00000000077) | ((x >> 12) & 00000007700);  // 0 000 007 777
  return x;
}

// 10 bit のデータを 3 bit 単位で 30 bit に拡張する
// in : 000 000 000 000 000 000 00x xxx xxx xxx (8)
// out: 00x 00x 00x 00x 00x 00x 00x 00x 00x 00x (8)
inline int pw_expand3(int x) {                        // 0 000 007 777
  x = (x & 00000000077) | ((x & 00000007700) << 12);  // 0 077 000 077
  x = (x & 00007000007) | ((x & 00070000070) << 6);   // 7 007 007 007
  x = (x & 01001001001) | ((x & 02222222222) << 2) |
      ((x & 04444444444) << 4);  // 1 111 111 111
  return x;
}

// 空間座標からハッシュ値 (モートンコード) を得る
inline int morton_code(int const x, int const y, int const z) {
  return pw_expand3(x) | (pw_expand3(y) << 1) | (pw_expand3(z) << 2);
}
}  // end of namespace hi
