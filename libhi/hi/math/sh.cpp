#include <hi/math/sh.hpp>

namespace {
/// ルジャンドルの倍関数
template <typename T>
T P(int const l, int const m, T const x) {
  T pmm = T(1);

  if (0 < m) {
    T const somx2 = std::sqrt((1 - x) * (1 + x));
    T fact = 1;

    for (int i = 1; i <= m; ++i) {
      pmm *= (-fact) * somx2;
      fact += 2;
    }
  }

  if (l == m) {
    return pmm;
  }

  T pmmp1 = x * (2 * m + 1) * pmm;

  if (l == (m + 1)) {
    return pmmp1;
  }

  T pll = 0;

  for (int ll = m + 2; ll <= l; ++ll) {
    pll = ((2 * ll - 1) * x * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
    pmm = pmmp1;
    pmmp1 = pll;
  }

  return pll;
}

/// 球面調和関数の正規化項
template <typename T>
T K(int const l, int m) {
  T lpm = 1;
  T lnm = 1;

  if (m < 0) {
    m = -m;
  }

  for (int i = l - m; 0 < i; --i) {
    lnm *= i;
  }
  for (int i = l + m; 0 < i; --i) {
    lpm *= i;
  }

  return std::sqrt(((2 * l + 1) * lnm) / (4 * T(M_PI) * lpm));
}
}

namespace hi {
#define MATERIALIZE(T)                                                         \
  template <>                                                                  \
  T SphericalHamonics(int const l, int const m, T const theta, T const phi) {  \
    return K<T>(l, m) *                                                        \
           ((0 < m)                                                            \
                ? T(M_SQRT2) * std::cos(m * phi) * P<T>(l, m, std::cos(theta)) \
                : (m < 0)                                                      \
                      ? T(M_SQRT2) * std::sin(-m * phi) *                      \
                            P<T>(l, -m, std::cos(theta))                       \
                      : P<T>(l, m, std::cos(theta)));                          \
  }
MATERIALIZE(float)
MATERIALIZE(double)
MATERIALIZE(long double)
#undef MATERIALIZE
}
// end of namespace hi
