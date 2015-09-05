// rgb2spd.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

#include "stdafx.h"

#define DIMENSION 31

namespace {
/// フーリエ級数の七次までの基底を初期化する．
void InitializeFourierSeries(__notnull hi::basic_spectrum *b) {
  for (int k = 0; k < DIMENSION; ++k) {
    for (std::size_t i = 0; i < hi::basic_spectrum::Size; ++i) {
      // x in [-pi, pi]
      hi::basic_spectrum::value_type const x =
          (2.0 * M_PI) * ((i + 0.5) / hi::basic_spectrum::Size - 0.5);
      b[k][i] = (k & 1) ? std::sin(x * (k / 2 + 1)) : std::cos(x * (k / 2));
    }

    // 正規化
    b[k] = hi::normalize(b[k]);
  }
}

inline int DiracDelta(int const x) { return x ? 0 : 1; }

inline double CosineSeriesNormalizationTerm(int const i) {
  return std::sqrt(hi::basic_spectrum::Size / (2.0 - DiracDelta(i)));
}

inline double triple_product_coefficient(int const i, int const j,
                                         int const k) {
  return (DiracDelta(i + j + k)     // i = j = k = 0
          + DiracDelta(i + j - k)   // i + j = k
          + DiracDelta(i - j + k)   // i + k = j
          + DiracDelta(i - j - k))  // i = j + k
         *
         hi::basic_spectrum::Size / 4.0 / CosineSeriesNormalizationTerm(i) /
         CosineSeriesNormalizationTerm(j) / CosineSeriesNormalizationTerm(k);
}

inline double triple_product_coefficient(hi::basic_spectrum const &i,
                                         hi::basic_spectrum const &j,
                                         hi::basic_spectrum const &k) {
  double sum = 0;
  for (int x = 0; x < hi::basic_spectrum::Size; ++x) {
    sum += i[x] * j[x] * k[x];
  }
  return sum;
}

void InitializeCosineSeries(__notnull hi::basic_spectrum *b) {
  for (int k = 0; k < DIMENSION; ++k) {
    for (std::size_t i = 0; i < hi::basic_spectrum::Size; ++i) {
      // x in [0, pi]
      hi::basic_spectrum::value_type const x =
          M_PI * (i + 0.5) / hi::basic_spectrum::Size;
      b[k][i] = std::cos(k * x);
    }
    b[k] *= 1.0 / CosineSeriesNormalizationTerm(k);  // 正規化
    // b[k] = hi::normalize(b[k]);
  }
}

///
/// ルジャンドル多項式の再帰的定義
///             P_{0}(x) = 1
///             P_{1}(x) = x
///   (k+1) * P_{k+1}(x) = (2 * k + 1) * x * P_{k}(x) - k * P_{k-1}(x)
///
double LegendrePolynomial(int const k, double const x) {
  if (0 == k) {
    return 1;
  }

  if (1 == k) {
    return x;
  }

  return ((2 * (k - 1) + 1) * x * LegendrePolynomial(k - 1, x) -
          (k - 1) * LegendrePolynomial(k - 2, x)) /
         k;
}

/// ルジャンドルの多項式の七次までの基底を初期化する．
void InitializeLegendrePolynomial(__notnull hi::basic_spectrum *b) {
  for (std::size_t i = 0; i < hi::basic_spectrum::Size; ++i) {
    // x in [-1, 1]
    hi::basic_spectrum::value_type const x =
        2 * (i + 0.5) / hi::basic_spectrum::Size - 1;

    for (int k = 0; k < DIMENSION; ++k) {
      b[k][i] = LegendrePolynomial(k, x);
    }
    /*
          hi::basic_spectrum::value_type const xSq = hi::square_of(x);
          b[0][i] = 1;
          b[1][i] = x;
          b[2][i] = (3 * xSq - 1) * 0.5;
          b[3][i] = (5 * xSq - 3) * 0.5 * x;
          b[4][i] = ((35 * xSq - 30) * xSq +  3) * 0.125;
          b[5][i] = ((63 * xSq - 70) * xSq + 15) * 0.125 * x;
          b[6][i] = (((231 * xSq - 315) * xSq + 105) * xSq -  5) * 0.0625;
          b[7][i] = (((429 * xSq - 639) * xSq + 315) * xSq - 35) * 0.0625 * x;
          b[8][i] = ((((6435 * xSq - 12012) * xSq + 6930) * xSq - 1260) * xSq +
       35) * 0.0078125;
    */
  }

  // 正規化
  for (std::size_t n = 0; n < DIMENSION; ++n) {
    b[n] = hi::normalize(b[n]);
  }
}

///
/// チェビシェフ多項式の再帰的な定義[-1, 1]
///     T_{0}(x) = 1
///     T_{1}(x) = x
///   T_{k+1}(x) = 2 * x * T_{k}(x) - T_{k-1}(x)
///
///
double ChebyshevPolynomials(int const k, double const x) {
  if (0 == k) {
    return 1;
  }

  if (1 == k) {
    return x;
  }

  return 2 * x * ChebyshevPolynomials(k - 1, x) -
         ChebyshevPolynomials(k - 2, x);
}

inline double ChebyshevPolynomialWeight(double const &x) {
  return hi::rcp(std::sqrt(1 - hi::square_of(x)));
}

/// チェビシェフの多項式の七次までの基底を初期化する．
/// チェビシェフの多項式は，重み w(x) = 1/sqrt(1-x^2) の基で直交関数系である．
void InitializeChebyshevPolynomial(__notnull hi::basic_spectrum *b) {
  for (std::size_t i = 0; i < hi::basic_spectrum::Size; ++i) {
    // x in (-1, 1)
    hi::basic_spectrum::value_type const x =
        2 * (i + 0.5) / hi::basic_spectrum::Size - 1;

    for (int k = 0; k < DIMENSION; ++k) {
      b[k][i] = ChebyshevPolynomials(k, x);
    }
    /*
          hi::basic_spectrum::value_type const xSq = hi::square_of(x);
          b0[i] = 1;
          b1[i] = x;
          b2[i] = (( 2 * xSq -  1)           );
          b3[i] = (( 4 * xSq -  3) * x       );
          b4[i] = (( 8 * xSq -  8) * xSq +  1);
          b5[i] = ((16 * xSq - 20) * xSq +  5) * x;
          b6[i] = ((32 * xSq - 48) * xSq + 18) * xSq - 1;
    */
  }

  // 正規化
  // hi::basic_spectrum * const b[] = { &b0, &b1, &b2, &b3, &b4, &b5, &b6 };
  for (std::size_t n = 0; n < DIMENSION; ++n) {
    hi::basic_spectrum::value_type value = 0;
    for (std::size_t i = 0; i < hi::basic_spectrum::Size; ++i) {
      hi::basic_spectrum::value_type const x =
          2 * (i + 0.5) / hi::basic_spectrum::Size - 1;
      value += hi::square_of(b[n][i]) * ChebyshevPolynomialWeight(x);
    }
    b[n] *= hi::rsqrt(value);
  }
}

/// 選点直交多項式の多項式の七次までの基底を初期化する．
void InitializeSelectedPoints(__notnull hi::basic_spectrum *b) {
  int const N = static_cast<int>(hi::basic_spectrum::Size) - 1;
  for (int k = 0; k < DIMENSION; ++k) {
    for (int x = 0; x <= N; ++x) {
      b[k][x] = 0;
      for (int a = 0; a <= k; ++a) {
        b[k][x] += ((a & 1) ? -1.0 : 1.0) *
                   hi::binomialCoefficient((double)k, (double)a) *
                   hi::binomialCoefficient((double)k + a, (double)a) *
                   hi::discretePow((double)x, (double)a) /
                   hi::discretePow((double)N, (double)a);
      }
    }
  }
#if 1
  // 正規化
  for (std::size_t k = 0; k < DIMENSION; ++k) {
    b[k] = hi::normalize(b[k]);
  }
#else
  // 正規化
  for (int k = 0; k < DIMENSION; ++k) {
    b[k] *= hi::rcp(std::sqrt(hi::discretePow(N + k + 1, k + 1) /
                              (hi::discretePow(N, k) * double(2 * k + 1))));
  }
#endif
}

void FactrizeD2(hi::basic_spectrum const &data, std::vector<double> &result) {
  std::vector<double> a(32, 0.0);
  std::vector<double> b(32, 0.0);

  std::copy(&data.data()[0], &data.data()[hi::basic_spectrum::Size], a.begin());
  a[hi::basic_spectrum::Size] = data[hi::basic_spectrum::Size - 1];

  for (int w = 32; w > 1; w /= 2) {
    int const o = w / 2;
    for (int i = 0; i < w; i += 2) {
      b[i / 2] = (a[i] + a[i + 1]) / 2;
      b[i / 2 + o] = (a[i] - a[i + 1]) / 2;
    }
    std::copy(b.begin(), b.begin() + w, a.begin());
  }
  a.swap(result);  // 結果を設定
}

// k個の係数で再構築する
void ReconstructD2(int const k, std::vector<double> const &data,
                   hi::basic_spectrum &result) {
  std::vector<double> a(data);
  std::vector<double> b(32, 0.0);

  // k個の係数だけ残す
  std::size_t flag = ~0U;
  for (int i = k; i < 32; ++i) {
    double minValue = std::numeric_limits<double>::infinity();
    int minIndex = 0;
    for (int j = 0; j < 32; ++j) {
      if ((flag & (1 << j)) && (minValue > std::abs(a[j]))) {
        minValue = std::abs(a[j]);
        minIndex = j;
      }
    }
    flag &= ~(1 << minIndex);
    a[minIndex] = 0;
  }

  for (int w = 2; w <= 32; w *= 2) {
    int const o = w / 2;
    for (int i = 0; i < w; i += 2) {
      b[i] = a[i / 2] + a[o + i / 2];
      b[i + 1] = a[i / 2] - a[o + i / 2];
    }
    std::copy(b.begin(), b.begin() + w, a.begin());
  }

  std::copy(a.begin(), a.begin() + hi::basic_spectrum::Size, result.data());
}

double const hd4[] = {
    (1 + std::sqrt(3.0)) / (4 * std::sqrt(2.0)),
    (3 + std::sqrt(3.0)) / (4 * std::sqrt(2.0)),
    (3 - std::sqrt(3.0)) / (4 * std::sqrt(2.0)),
    (1 - std::sqrt(3.0)) / (4 * std::sqrt(2.0)),
};
double const gd4[] = {hd4[3], -hd4[2], hd4[1], -hd4[0]};

// ドビッシーのマザーウェーブレット(D4)を使ったウェーブレット変換．
void FactrizeD4(hi::basic_spectrum const &data, std::vector<double> &result) {
  std::vector<double> a(32, 0.0);
  std::vector<double> b(32, 0.0);

  std::copy(&data.data()[0], &data.data()[hi::basic_spectrum::Size], a.begin());
  a[hi::basic_spectrum::Size] = data[hi::basic_spectrum::Size - 1];

  for (int w = 32; w > 1; w /= 2) {
    int const o = w / 2;
    for (int i = 0; i < w; i += 2) {
      b[i / 2] = hd4[0] * a[i + 0] + hd4[1] * a[i + 1] +
                 hd4[2] * a[(i + 2) % w] + hd4[3] * a[(i + 3) % w];
      b[i / 2 + o] = gd4[0] * a[i + 0] + gd4[1] * a[i + 1] +
                     gd4[2] * a[(i + 2) % w] + gd4[3] * a[(i + 3) % w];
    }
    std::copy(b.begin(), b.begin() + w, a.begin());
  }
  a.swap(result);  // 結果を設定
}

// k個の係数で再構築する
void ReconstructD4(int const k, std::vector<double> const &data,
                   hi::basic_spectrum &result) {
  std::vector<double> a(data);
  std::vector<double> b(32, 0.0);

  // k個の係数だけ残す
  std::size_t flag = ~0;
  for (int i = k; i < 32; ++i) {
    double minValue = std::numeric_limits<double>::infinity();
    int minIndex = 0;
    for (int j = 0; j < 32; ++j) {
      if ((flag & (1 << j)) && (minValue > std::abs(a[j]))) {
        minValue = std::abs(a[j]);
        minIndex = j;
      }
    }
    flag &= ~(1 << minIndex);
    a[minIndex] = 0;
  }

  for (int w = 2; w <= 32; w *= 2) {
    int const o = w / 2;
    for (int i = 0; i < w; i += 2) {
      int l = i ? (i / 2 - 1) : o - 1;
      b[i] = hd4[0] * a[i / 2] + hd4[2] * a[l] + gd4[0] * a[o + i / 2] +
             gd4[2] * a[o + l];
      b[i + 1] = hd4[1] * a[i / 2] + hd4[3] * a[l] + gd4[1] * a[o + i / 2] +
                 gd4[3] * a[o + l];
    }
    std::copy(b.begin(), b.begin() + w, a.begin());
  }
  std::copy(a.begin(), a.begin() + hi::basic_spectrum::Size, result.data());
}
}
// end of unnamed namespace

int _tmain(int argc, _TCHAR *argv[])
#if 0
// f(theta, phi) = 1 の関数をLight Probe Image領域で積分する．
{
  int const N = 100;

  double sum_of_f = 0;
  for (int y = 0; y < N; ++y)
  {
    for (int x = 0; x < N; ++x)
    {
      double const u = 2 * (y + 0.5) / N - 1;
      double const v = 2 * (x + 0.5) / N - 1;
      double const d = u * u + v * v;

      if (d > 1)
      {
        continue;
      }

      double const theta = std::atan2(v, u);
      double const phi = M_PI * std::sqrt(d);
      double const jacobian = std::sin(phi) * M_PI * hi::rsqrt(d);

      sum_of_f += jacobian;
    }
  }

  ::_ftprintf_s(stdout, _TEXT("%g\n"), 4 * M_PI);
  ::_ftprintf_s(stdout, _TEXT("%g\n"), sum_of_f * ((2.0*2.0)/(N*N)));

  return 0;
}

#elif 0
// COS級数の近似の正規直交行列を計算する
{
  std::vector<double> basis[8];

  for (int n = 0; n < 8; ++n) {
    basis[n].resize(8);

    double sum = 0;
    for (int i = 0; i < 8; ++i) {
      hi::basic_spectrum::value_type const x = M_PI * (i + 0.5) / 8;
      basis[n][i] = std::cos(n * x);
      sum += hi::square_of(basis[n][i]);
    }

    // 正規化
    sum = hi::rsqrt(sum);
    for (int i = 0; i < 8; ++i) {
      basis[n][i] *= sum;
    }
  }

  for (int n = 0; n < 8; ++n) {
    for (int i = 0; i < 8; ++i) {
      ::_ftprintf_s(stdout, _TEXT(" %g, "), basis[n][i]);
    }
    ::_ftprintf_s(stdout, _TEXT("\n"));
  }

  return 0;
}
#elif 0
//
// 等色関数を二次元に拡張したものの級数展開の係数を求める
//
{
  hi::basic_spectrum const XYZ_CMF[3] = {
      hi::basic_spectrum(
          0.014 / 9.80328, 0.044 / 9.80328, 0.134 / 9.80328, 0.284 / 9.80328,
          0.348 / 9.80328, 0.336 / 9.80328, 0.291 / 9.80328, 0.195 / 9.80328,
          0.096 / 9.80328, 0.032 / 9.80328, 0.005 / 9.80328, 0.009 / 9.80328,
          0.063 / 9.80328, 0.166 / 9.80328, 0.290 / 9.80328, 0.433 / 9.80328,
          0.595 / 9.80328, 0.762 / 9.80328, 0.916 / 9.80328, 1.026 / 9.80328,
          1.062 / 9.80328, 1.003 / 9.80328, 0.854 / 9.80328, 0.642 / 9.80328,
          0.448 / 9.80328, 0.284 / 9.80328, 0.165 / 9.80328, 0.087 / 9.80328,
          0.047 / 9.80328, 0.023 / 9.80328, 0.011 / 9.80328),
      hi::basic_spectrum(
          0.000 / 9.80328, 0.001 / 9.80328, 0.004 / 9.80328, 0.012 / 9.80328,
          0.023 / 9.80328, 0.038 / 9.80328, 0.060 / 9.80328, 0.091 / 9.80328,
          0.139 / 9.80328, 0.208 / 9.80328, 0.323 / 9.80328, 0.503 / 9.80328,
          0.710 / 9.80328, 0.862 / 9.80328, 0.954 / 9.80328, 0.995 / 9.80328,
          0.995 / 9.80328, 0.952 / 9.80328, 0.870 / 9.80328, 0.757 / 9.80328,
          0.631 / 9.80328, 0.503 / 9.80328, 0.381 / 9.80328, 0.265 / 9.80328,
          0.175 / 9.80328, 0.107 / 9.80328, 0.061 / 9.80328, 0.032 / 9.80328,
          0.017 / 9.80328, 0.008 / 9.80328, 0.004 / 9.80328),
      hi::basic_spectrum(
          0.068 / 9.80328, 0.207 / 9.80328, 0.646 / 9.80328, 1.386 / 9.80328,
          1.747 / 9.80328, 1.772 / 9.80328, 1.669 / 9.80328, 1.288 / 9.80328,
          0.813 / 9.80328, 0.465 / 9.80328, 0.272 / 9.80328, 0.158 / 9.80328,
          0.078 / 9.80328, 0.042 / 9.80328, 0.020 / 9.80328, 0.009 / 9.80328,
          0.004 / 9.80328, 0.002 / 9.80328, 0.002 / 9.80328, 0.001 / 9.80328,
          0.001 / 9.80328, 0.000 / 9.80328, 0.000 / 9.80328, 0.000 / 9.80328,
          0.000 / 9.80328, 0.000 / 9.80328, 0.000 / 9.80328, 0.000 / 9.80328,
          0.000 / 9.80328, 0.000 / 9.80328, 0.000 / 9.80328),
  };

  hi::basic_spectrum RGB_CMF[3];
  for (std::size_t i = 0; i < hi::basic_spectrum::Size; ++i) {
    RGB_CMF[0][i] = 1.9104 * XYZ_CMF[0][i] - 0.5338 * XYZ_CMF[1][i] -
                    0.2891 * XYZ_CMF[2][i];
    RGB_CMF[1][i] = -0.9844 * XYZ_CMF[0][i] + 1.9985 * XYZ_CMF[1][i] -
                    0.0279 * XYZ_CMF[2][i];
    RGB_CMF[2][i] = 0.0585 * XYZ_CMF[0][i] - 0.1187 * XYZ_CMF[1][i] +
                    0.9017 * XYZ_CMF[2][i];
  }

  hi::basic_spectrum basis_vectors[hi::basic_spectrum::Size];
  InitializeCosineSeries(hi::basic_spectrum::Size, basis_vectors);

  for (int k = 0; k < 3; ++k) {
    hi::basic_spectrum matrix0[hi::basic_spectrum::Size];
    for (int i = 0; i < hi::basic_spectrum::Size; ++i) {
      for (int j = 0; j < hi::basic_spectrum::Size; ++j) {
        matrix0[i][i] = 0;
      }
      matrix0[i][i] = RGB_CMF[k][i];
    }

    hi::basic_spectrum matrix1[hi::basic_spectrum::Size];
    hi::basic_spectrum matrix2[hi::basic_spectrum::Size];

    // 第一段階の係数を求めて転置する
    for (int i = 0; i < hi::basic_spectrum::Size; ++i) {
      for (int j = 0; j < hi::basic_spectrum::Size; ++j) {
        matrix1[j][i] = hi::dot(matrix0[i], basis_vectors[j]);
      }
    }

    // 第二段階の係数を求めて転置する
    for (int i = 0; i < hi::basic_spectrum::Size; ++i) {
      for (int j = 0; j < hi::basic_spectrum::Size; ++j) {
        matrix2[j][i] = hi::dot(matrix1[i], basis_vectors[j]);
      }
    }

    ::_ftprintf_s(stdout, _TEXT("Matrix[%d]\n"), k);
    for (int k = 0; k < 8; ++k) {
      for (int i = 0; i < 8; ++i) {
        ::_ftprintf_s(stdout, _TEXT(" %g, "), matrix2[k][i]);
      }
      ::_ftprintf_s(stdout, _TEXT("\n"));
    }
  }

  return 0;
}
#elif 0
//
// 等色関数の級数展開の係数を求める
//
{
  hi::basic_spectrum const XYZ_CMF[3] = {
      hi::basic_spectrum(
          0.014 / 9.80328, 0.044 / 9.80328, 0.134 / 9.80328, 0.284 / 9.80328,
          0.348 / 9.80328, 0.336 / 9.80328, 0.291 / 9.80328, 0.195 / 9.80328,
          0.096 / 9.80328, 0.032 / 9.80328, 0.005 / 9.80328, 0.009 / 9.80328,
          0.063 / 9.80328, 0.166 / 9.80328, 0.290 / 9.80328, 0.433 / 9.80328,
          0.595 / 9.80328, 0.762 / 9.80328, 0.916 / 9.80328, 1.026 / 9.80328,
          1.062 / 9.80328, 1.003 / 9.80328, 0.854 / 9.80328, 0.642 / 9.80328,
          0.448 / 9.80328, 0.284 / 9.80328, 0.165 / 9.80328, 0.087 / 9.80328,
          0.047 / 9.80328, 0.023 / 9.80328, 0.011 / 9.80328),
      hi::basic_spectrum(
          0.000 / 9.80328, 0.001 / 9.80328, 0.004 / 9.80328, 0.012 / 9.80328,
          0.023 / 9.80328, 0.038 / 9.80328, 0.060 / 9.80328, 0.091 / 9.80328,
          0.139 / 9.80328, 0.208 / 9.80328, 0.323 / 9.80328, 0.503 / 9.80328,
          0.710 / 9.80328, 0.862 / 9.80328, 0.954 / 9.80328, 0.995 / 9.80328,
          0.995 / 9.80328, 0.952 / 9.80328, 0.870 / 9.80328, 0.757 / 9.80328,
          0.631 / 9.80328, 0.503 / 9.80328, 0.381 / 9.80328, 0.265 / 9.80328,
          0.175 / 9.80328, 0.107 / 9.80328, 0.061 / 9.80328, 0.032 / 9.80328,
          0.017 / 9.80328, 0.008 / 9.80328, 0.004 / 9.80328),
      hi::basic_spectrum(
          0.068 / 9.80328, 0.207 / 9.80328, 0.646 / 9.80328, 1.386 / 9.80328,
          1.747 / 9.80328, 1.772 / 9.80328, 1.669 / 9.80328, 1.288 / 9.80328,
          0.813 / 9.80328, 0.465 / 9.80328, 0.272 / 9.80328, 0.158 / 9.80328,
          0.078 / 9.80328, 0.042 / 9.80328, 0.020 / 9.80328, 0.009 / 9.80328,
          0.004 / 9.80328, 0.002 / 9.80328, 0.002 / 9.80328, 0.001 / 9.80328,
          0.001 / 9.80328, 0.000 / 9.80328, 0.000 / 9.80328, 0.000 / 9.80328,
          0.000 / 9.80328, 0.000 / 9.80328, 0.000 / 9.80328, 0.000 / 9.80328,
          0.000 / 9.80328, 0.000 / 9.80328, 0.000 / 9.80328),
  };

  hi::basic_spectrum RGB_CMF[3];
  for (std::size_t i = 0; i < hi::basic_spectrum::Size; ++i) {
    RGB_CMF[0][i] = 1.9104 * XYZ_CMF[0][i] - 0.5338 * XYZ_CMF[1][i] -
                    0.2891 * XYZ_CMF[2][i];
    RGB_CMF[1][i] = -0.9844 * XYZ_CMF[0][i] + 1.9985 * XYZ_CMF[1][i] -
                    0.0279 * XYZ_CMF[2][i];
    RGB_CMF[2][i] = 0.0585 * XYZ_CMF[0][i] - 0.1187 * XYZ_CMF[1][i] +
                    0.9017 * XYZ_CMF[2][i];
  }

  hi::basic_spectrum basis_vectors[8];
  InitializeCosineSeries(8, basis_vectors);

  hi::spectrum_vector c[8];
  for (int i = 0; i < 8; ++i) {
    c[i][0] = hi::dot(RGB_CMF[0], basis_vectors[i]);
    c[i][1] = hi::dot(RGB_CMF[1], basis_vectors[i]);
    c[i][2] = hi::dot(RGB_CMF[2], basis_vectors[i]);
  }

  for (int i = 0; i < 8; ++i) {
    ::_ftprintf_s(stdout, _TEXT("%d, %g, %g, %g\n"), i, c[i][0], c[i][1],
                  c[i][2]);
  }

  return 0;
}
#elif 0
//
// RGB -> SPD 変換
//
{
  hi::spectrum_vector src(255.5 / 256);
  hi::basic_spectrum spd0;
  hi::rgb2spd(src, spd0);

  hi::basic_spectrum spd1(spd0);
  hi::clamp(spd1,
            std::numeric_limits<hi::basic_spectrum::value_type>::infinity());

  for (int i = 0; i < hi::basic_spectrum::Size; ++i) {
    ::_ftprintf_s(stdout, _TEXT("%d, %g, %g\n"), 400 + 10 * i, spd0[i],
                  spd1[i]);
  }

  hi::spectrum_vector xyz;
  hi::spd2xyz(spd1, xyz);

  hi::spectrum_vector dst;
  hi::CCIR601_1_XYZ2RGB(xyz.data(), dst.data());

  ::_ftprintf_s(stderr, _TEXT("src = % 9.8lf, % 9.8lf, % 9.8lf\n"), src[0],
                src[1], src[2]);
  ::_ftprintf_s(stderr, _TEXT("xyz = % 9.8lf, % 9.8lf, % 9.8lf\n"), xyz[0],
                xyz[1], xyz[2]);
  ::_ftprintf_s(stderr, _TEXT("dst = % 9.8lf, % 9.8lf, % 9.8lf\n"), dst[0],
                dst[1], dst[2]);
  ::_ftprintf_s(
      stderr, _TEXT("      % 9.8lf\n"),
      std::sqrt(hi::length_squared(dst - src) / hi::length_squared(src)));

  return 0;
}
#elif 0
//
// 三重積の計算
//
{
  // 級数展開のために汎用的に利用する基底ベクトル
  hi::basic_spectrum basis_vectors[DIMENSION];
  InitializeCosineSeries(basis_vectors);

  // 元データ
  hi::basic_spectrum e = hi::basic_spectrum::MccYellowGreen();
  hi::basic_spectrum f = hi::basic_spectrum::MccGreen();

  // 元データの級数展開
  std::vector<double> co_e(DIMENSION, 0.0);
  std::vector<double> co_f(DIMENSION, 0.0);
  for (int i = 0; i < DIMENSION; ++i) {
    co_e[i] = hi::dot(e, basis_vectors[i]);
    co_f[i] = hi::dot(f, basis_vectors[i]);
  }

  // 正解データ
  hi::basic_spectrum g = e * f;

  // 正解データの級数近似
  std::vector<double> co_a(DIMENSION, 0.0);
  for (int i = 0; i < DIMENSION; ++i) {
    co_a[i] = hi::dot(g, basis_vectors[i]);
  }

  // 元データの級数展開データから三重積により近似係数を計算
  std::vector<double> co_b(DIMENSION, 0.0);
  {
    // i = 0 の場合
    {
      double sum = 0;
      for (int l = 0; l < DIMENSION; ++l) {
        sum += triple_product_coefficient(0, l, l) * co_e[l] * co_f[l];
      }
      co_b[0] = sum;
    }

    // 0 < i < N-1 の場合
    for (int i = 1; i < DIMENSION - 1; ++i) {
      double sum = 0;

      // 二重線の部分
      for (int l = i; l < DIMENSION; ++l) {
        sum += triple_product_coefficient(i, l - i, l) * co_e[l - i] * co_f[l];
        sum += triple_product_coefficient(i, l, l - i) * co_e[l] * co_f[l - i];
      }

      // 対角線の分
      for (int l = 1; l < i; ++l) {
        sum += triple_product_coefficient(i, i - l, l) * co_e[i - l] * co_f[l];
      }

      co_b[i] = sum;
    }

    // i = N-1 の場合
    {
      double sum = 0;
      for (int l = 0; l < DIMENSION; ++l) {
        sum += triple_product_coefficient(DIMENSION - 1, DIMENSION - 1 - l, l) *
               co_e[DIMENSION - 1 - l] * co_f[l];
      }
      co_b[DIMENSION - 1] = sum;
    }
  }

  // 偏光の計算の形を利用して近似係数を求める
  //   S(lambda) = \int_{380}^{780} h(lambda, lambda') S(lambda') du(lambda')
  std::vector<double> co_c(DIMENSION, 0.0);
  {
    hi::basic_spectrum matrix0[hi::basic_spectrum::Size];
    for (int i = 0; i < hi::basic_spectrum::Size; ++i) {
      for (int j = 0; j < hi::basic_spectrum::Size; ++j) {
        matrix0[i][i] = 0;
      }
      matrix0[i][i] = f[i];
    }

    hi::basic_spectrum matrix1[hi::basic_spectrum::Size];
    hi::basic_spectrum matrix2[hi::basic_spectrum::Size];

    // 第一段階の係数を求めて転置する
    for (int i = 0; i < hi::basic_spectrum::Size; ++i) {
      for (int j = 0; j < hi::basic_spectrum::Size; ++j) {
        matrix1[j][i] = hi::dot(matrix0[i], basis_vectors[j]);
      }
    }

    // 第二段階の係数を求めて転置する
    for (int i = 0; i < hi::basic_spectrum::Size; ++i) {
      for (int j = 0; j < hi::basic_spectrum::Size; ++j) {
        matrix2[j][i] = hi::dot(matrix1[i], basis_vectors[j]);
      }
    }

    ::_ftprintf_s(stderr, _TEXT("Matrix\n"));
    for (int k = 0; k < DIMENSION; ++k) {
      co_c[k] = 0;
      for (int i = 0; i < DIMENSION; ++i) {
        ::_ftprintf_s(stderr, _TEXT(" % 8.5lf"), matrix2[k][i]);
        co_c[k] += matrix2[k][i] * co_e[i];
      }
      ::_ftprintf_s(stderr, _TEXT("\n"));
    }
  }

  // 計算した係数の表示
  ::_ftprintf_s(stdout, _TEXT("co_a, co_b, co_c\n"));
  for (std::size_t i = 0; i < DIMENSION; ++i) {
    ::_ftprintf_s(stdout, _TEXT("%g, %g, %g\n"), co_a[i], co_b[i], co_c[i]);
  }

  // 近似係数から波形を再構成
  hi::basic_spectrum approx_a(0.0);
  hi::basic_spectrum approx_b(0.0);
  hi::basic_spectrum approx_c(0.0);
  {
    hi::basic_spectrum vector;
    for (int i = 0; i < DIMENSION; ++i) {
      vector = basis_vectors[i];
      vector *= co_a[i];
      approx_a += vector;

      vector = basis_vectors[i];
      vector *= co_b[i];
      approx_b += vector;

      vector = basis_vectors[i];
      vector *= co_c[i];
      approx_c += vector;
    }
  }

  ::_ftprintf_s(
      stdout,
      _TEXT(
          "e, f, g=e*f, approx_a, approx_b, approx_c, abs_a, abs_b, abs_c\n"));
  for (int i = 0; i < hi::basic_spectrum::Size; ++i) {
    ::_ftprintf_s(stdout, _TEXT("%g, %g, %g, %g, %g, %g, %g, %g, %g\n"), e[i],
                  f[i], g[i], approx_a[i], approx_b[i], approx_c[i],
                  std::abs(approx_a[i] - g[i]), std::abs(approx_b[i] - g[i]),
                  std::abs(approx_c[i] - g[i]));
  }

  return 0;
}
#elif 1
//
// 誤差のテスト
//
{
  hi::basic_spectrum f;
  /*
    for (int i = 0; i < 24; ++i)
    {
      switch (i)
      {
      case  0: f = hi::basic_spectrum::MccDarkSkin();     ::_ftprintf_s(stdout,
  _TEXT("mmc-dark-skin\n"));     break;
      case  1: f = hi::basic_spectrum::MccLightSkin();    ::_ftprintf_s(stdout,
  _TEXT("mmc-light-skin\n"));    break;
      case  2: f = hi::basic_spectrum::MccBlueSky();      ::_ftprintf_s(stdout,
  _TEXT("mmc-blue-sky\n"));      break;
      case  3: f = hi::basic_spectrum::MccFoliage();      ::_ftprintf_s(stdout,
  _TEXT("mmc-foliage\n"));       break;
      case  4: f = hi::basic_spectrum::MccBlueFlower();   ::_ftprintf_s(stdout,
  _TEXT("mmc-blue-flower\n"));   break;
      case  5: f = hi::basic_spectrum::MccBluishGreen();  ::_ftprintf_s(stdout,
  _TEXT("mmc-bluish-green\n"));  break;
      case  6: f = hi::basic_spectrum::MccOrange();       ::_ftprintf_s(stdout,
  _TEXT("mmc-orange\n"));        break;
      case  7: f = hi::basic_spectrum::MccPurplishBlue(); ::_ftprintf_s(stdout,
  _TEXT("mmc-purplish-blue\n")); break;
      case  8: f = hi::basic_spectrum::MccModerateRed();  ::_ftprintf_s(stdout,
  _TEXT("mmc-moderate-red\n"));  break;
      case  9: f = hi::basic_spectrum::MccPurple();       ::_ftprintf_s(stdout,
  _TEXT("mmc-purple\n"));        break;
      case 10: f = hi::basic_spectrum::MccYellowGreen();  ::_ftprintf_s(stdout,
  _TEXT("mmc-yellow-green\n"));  break;
      case 11: f = hi::basic_spectrum::MccOrangeYellow(); ::_ftprintf_s(stdout,
  _TEXT("mmc-orange-yellow\n")); break;
      case 12: f = hi::basic_spectrum::MccBlue();         ::_ftprintf_s(stdout,
  _TEXT("mmc-blue\n"));          break;
      case 13: f = hi::basic_spectrum::MccGreen();        ::_ftprintf_s(stdout,
  _TEXT("mmc-green\n"));         break;
      case 14: f = hi::basic_spectrum::MccRed();          ::_ftprintf_s(stdout,
  _TEXT("mmc-red\n"));           break;
      case 15: f = hi::basic_spectrum::MccYellow();       ::_ftprintf_s(stdout,
  _TEXT("mmc-yellow\n"));        break;
      case 16: f = hi::basic_spectrum::MccMagenta();      ::_ftprintf_s(stdout,
  _TEXT("mmc-magenta\n"));       break;
      case 17: f = hi::basic_spectrum::MccCyan();         ::_ftprintf_s(stdout,
  _TEXT("mmc-cyan\n"));          break;
      case 18: f = hi::basic_spectrum::MccWhite();        ::_ftprintf_s(stdout,
  _TEXT("mmc-white\n"));         break;
      case 19: f = hi::basic_spectrum::MccNeutral80();    ::_ftprintf_s(stdout,
  _TEXT("mmc-neutral.8\n"));     break;
      case 20: f = hi::basic_spectrum::MccNeutral65();    ::_ftprintf_s(stdout,
  _TEXT("mmc-neutral.65\n"));    break;
      case 21: f = hi::basic_spectrum::MccNeutral50();    ::_ftprintf_s(stdout,
  _TEXT("mmc-neutral.5\n"));     break;
      case 22: f = hi::basic_spectrum::MccNeutral35();    ::_ftprintf_s(stdout,
  _TEXT("mmc-neutral.35\n"));    break;
      case 23: f = hi::basic_spectrum::MccBlack();        ::_ftprintf_s(stdout,
  _TEXT("mmc-black\n"));         break;
      }
  /*/
  for (int i = 0; i < 12; ++i) {
    switch (i) {
    case 0:
      f = hi::basic_spectrum::LightSourceA();
      ::_ftprintf_s(stdout, _TEXT("Source A\n"));
      break;
    case 1:
      f = hi::basic_spectrum::LightSourceB();
      ::_ftprintf_s(stdout, _TEXT("Source B\n"));
      break;
    case 2:
      f = hi::basic_spectrum::LightSourceC();
      ::_ftprintf_s(stdout, _TEXT("Source C\n"));
      break;
    case 3:
      f = hi::basic_spectrum::LightD65();
      ::_ftprintf_s(stdout, _TEXT("D65\n"));
      break;
    case 4:
      f = hi::basic_spectrum::LightD100();
      ::_ftprintf_s(stdout, _TEXT("D100r\n"));
      break;
    case 5:
      f = hi::basic_spectrum::LightStandardWarmWhite();
      ::_ftprintf_s(stdout, _TEXT("Standard Warm White\n"));
      break;
    case 6:
      f = hi::basic_spectrum::LightWhite();
      ::_ftprintf_s(stdout, _TEXT("White\n"));
      break;
    case 7:
      f = hi::basic_spectrum::LightStandardWhite();
      ::_ftprintf_s(stdout, _TEXT("Standard White\n"));
      break;
    case 8:
      f = hi::basic_spectrum::LightDaylight();
      ::_ftprintf_s(stdout, _TEXT("Daylight\n"));
      break;
    case 9:
      f = hi::basic_spectrum::LightWarmWhiteDeLuxe();
      ::_ftprintf_s(stdout, _TEXT("Warm White DeLuxe\n"));
      break;
    case 10:
      f = hi::basic_spectrum::LightSoftWhite();
      ::_ftprintf_s(stdout, _TEXT("Soft White\n"));
      break;
    case 11:
      f = hi::basic_spectrum::LightCoolWhiteDeLuxe();
      ::_ftprintf_s(stdout, _TEXT("Cool White DeLuxe\n"));
      break;
    case 12:
      f = hi::basic_spectrum::LightMercuryArcLamp();
      ::_ftprintf_s(stdout, _TEXT("Mercury ArcLamp\n"));
      break;
    }
    //*/

    //
    // 基底を初期化する
    //
    hi::basic_spectrum b_fs[DIMENSION];
    InitializeFourierSeries(b_fs);
    hi::basic_spectrum b_cs[DIMENSION];
    InitializeCosineSeries(b_cs);
    hi::basic_spectrum b_lp[DIMENSION];
    InitializeLegendrePolynomial(b_lp);
    hi::basic_spectrum b_sp[DIMENSION];
    InitializeSelectedPoints(b_sp);
    hi::basic_spectrum b_cp[DIMENSION];
    InitializeChebyshevPolynomial(b_cp);

    //
    // 係数
    //
    double c_fs[DIMENSION];
    double c_cs[DIMENSION];
    double c_lp[DIMENSION];
    double c_sp[DIMENSION];
    double c_cp[DIMENSION];
    std::vector<double> c_d2;
    std::vector<double> c_d4;

    // 係数を求める
    for (std::size_t k = 0; k < DIMENSION; ++k) {
      c_fs[k] = hi::dot(f, b_fs[k]);
      c_cs[k] = hi::dot(f, b_cs[k]);
      c_lp[k] = hi::dot(f, b_lp[k]);
      c_sp[k] = hi::dot(f, b_sp[k]);

      // チェビシェフ多項式の場合
      c_cp[k] = 0;
      for (std::size_t j = 0; j < hi::basic_spectrum::Size; ++j) {
        double const x = 2 * (j + 0.5) / hi::basic_spectrum::Size - 1;
        c_cp[k] += f[j] * b_cp[k][j] * ChebyshevPolynomialWeight(x);
      }
    }
    // ウェーブレット変換の場合
    FactrizeD2(f, c_d2);
    FactrizeD4(f, c_d4);

    //
    // 復元された関数
    //
    hi::basic_spectrum g_fs(0);
    hi::basic_spectrum g_cs(0);
    hi::basic_spectrum g_lp(0);
    hi::basic_spectrum g_sp(0);
    hi::basic_spectrum g_cp(0);
    hi::basic_spectrum g_d2(0);
    hi::basic_spectrum g_d4(0);

    typedef hi::basic_vector3<hi::basic_spectrum::value_type> vector_t;
    vector_t reference;
    hi::spd2xyz(f, reference);  // f -> XYZ
    // double const cScaleXYZ = 1 / hi::length_squared(reference);
    double const cScale = 1 / hi::length_squared(f);

    // 誤差: Sun (2000) の尺度
    vector_t e_fs(0);
    vector_t e_cs(0);
    vector_t e_lp(0);
    vector_t e_sp(0);
    vector_t e_cp(0);
    vector_t e_d2(0);
    vector_t e_d4(0);

#if 1
    ::_ftprintf_s(stdout, _TEXT(", FS, CS, LP, SP, CP, D2, D4\n"));
    for (int k = 0; k < DIMENSION; ++k) {
      // 再構築
      b_fs[k] *= c_fs[k];
      g_fs += b_fs[k];
      b_cs[k] *= c_cs[k];
      g_cs += b_cs[k];
      b_lp[k] *= c_lp[k];
      g_lp += b_lp[k];
      b_sp[k] *= c_sp[k];
      g_sp += b_sp[k];
      b_cp[k] *= c_cp[k];
      g_cp += b_cp[k];
      ReconstructD2(k + 1, c_d2, g_d2);  // WTの場合
      ReconstructD4(k + 1, c_d4, g_d4);  // D4の場合

      /*
            hi::spd2xyz(g_fs - f, e_fs);
            hi::spd2xyz(g_cs - f, e_cs);
            hi::spd2xyz(g_lp - f, e_lp);
            hi::spd2xyz(g_sp - f, e_sp);
            hi::spd2xyz(g_cp - f, e_cp);
            hi::spd2xyz(g_d2 - f, e_d2);
            hi::spd2xyz(g_d4 - f, e_d4);

          ::_ftprintf_s(stdout, _TEXT("%d, %g, %g, %g, %g, %g, %g, %g\n"), k+1
              , std::sqrt(hi::length_squared(e_fs) * cScaleXYZ)
              , std::sqrt(hi::length_squared(e_cs) * cScaleXYZ)
              , std::sqrt(hi::length_squared(e_lp) * cScaleXYZ)
              , std::sqrt(hi::length_squared(e_sp) * cScaleXYZ)
              , std::sqrt(hi::length_squared(e_cp) * cScaleXYZ)
              , std::sqrt(hi::length_squared(e_d2) * cScaleXYZ)
              , std::sqrt(hi::length_squared(e_d4) * cScaleXYZ)
             );
      /*/
      ::_ftprintf_s(stdout, _TEXT("%d, %g, %g, %g, %g, %g, %g, %g\n"), k + 1,
                    std::sqrt(hi::length_squared(g_fs - f) * cScale),
                    std::sqrt(hi::length_squared(g_cs - f) * cScale),
                    std::sqrt(hi::length_squared(g_lp - f) * cScale),
                    std::sqrt(hi::length_squared(g_sp - f) * cScale),
                    std::sqrt(hi::length_squared(g_cp - f) * cScale),
                    std::sqrt(hi::length_squared(g_d2 - f) * cScale),
                    std::sqrt(hi::length_squared(g_d4 - f) * cScale));
      //*/
    }
#else
    for (int k = 0; k < DIMENSION; ++k) {
      // 再構築
      b_fs[k] *= c_fs[k];
      g_fs += b_fs[k];
      b_cs[k] *= c_cs[k];
      g_cs += b_cs[k];
      b_lp[k] *= c_lp[k];
      g_lp += b_lp[k];
      b_sp[k] *= c_sp[k];
      g_sp += b_sp[k];
      b_cp[k] *= c_cp[k];
      g_cp += b_cp[k];
      ReconstructD2(k + 1, c_d2, g_d2);  // WTの場合
      ReconstructD4(k + 1, c_d4, g_d4);  // D4の場合

      ::_ftprintf_s(stdout, _TEXT("K=%d\n"), k + 1);
      ::_ftprintf_s(stdout, _TEXT("S, FS, CS, LP, SP, CP, D2, D4\n"));
      for (std::size_t i = 0; i < hi::basic_spectrum::Size; ++i) {
        ::_ftprintf_s(stdout, _TEXT("%g, %g, %g, %g, %g, %g, %g, %g\n"), f[i],
                      g_fs[i], g_cs[i], g_lp[i], g_sp[i], g_cp[i], g_d2[i],
                      g_d4[i]);
      }
    }
#endif
    ::_ftprintf_s(stdout, _TEXT("\n"));
  }

  return 0;
}
#else
{
  //
  // 直交関数近似
  //
  hi::basic_spectrum w(0.007161355, 0.021800007, 0.068032871, 0.145965261,
                       0.183983630, 0.186616481, 0.175769135, 0.135644485,
                       0.085620315, 0.048971029, 0.032948156, 0.051309358,
                       0.072424739, 0.087929754, 0.097314368, 0.101496642,
                       0.101496642, 0.097110355, 0.093930510, 0.105210375,
                       0.108901967, 0.102851857, 0.087572768, 0.065833392,
                       0.045939813, 0.029122560, 0.016919797, 0.008921348,
                       0.004819579, 0.002358517, 0.001127986);
  hi::basic_spectrum f;

  for (int i = 0; i < 25; ++i) {
    switch (i) {
    case 0:
      f = hi::basic_spectrum::MccDarkSkin();
      ::_ftprintf_s(stdout, _TEXT("mmc-dark-skin\n"));
      break;
    case 1:
      f = hi::basic_spectrum::MccLightSkin();
      ::_ftprintf_s(stdout, _TEXT("mmc-light-skin\n"));
      break;
    case 2:
      f = hi::basic_spectrum::MccBlueSky();
      ::_ftprintf_s(stdout, _TEXT("mmc-blue-sky\n"));
      break;
    case 3:
      f = hi::basic_spectrum::MccFoliage();
      ::_ftprintf_s(stdout, _TEXT("mmc-foliage\n"));
      break;
    case 4:
      f = hi::basic_spectrum::MccBlueFlower();
      ::_ftprintf_s(stdout, _TEXT("mmc-blue-flower\n"));
      break;
    case 5:
      f = hi::basic_spectrum::MccBluishGreen();
      ::_ftprintf_s(stdout, _TEXT("mmc-bluish-green\n"));
      break;
    case 6:
      f = hi::basic_spectrum::MccOrange();
      ::_ftprintf_s(stdout, _TEXT("mmc-orange\n"));
      break;
    case 7:
      f = hi::basic_spectrum::MccPurplishBlue();
      ::_ftprintf_s(stdout, _TEXT("mmc-purplish-blue\n"));
      break;
    case 8:
      f = hi::basic_spectrum::MccModerateRed();
      ::_ftprintf_s(stdout, _TEXT("mmc-moderate-red\n"));
      break;
    case 9:
      f = hi::basic_spectrum::MccPurple();
      ::_ftprintf_s(stdout, _TEXT("mmc-purple\n"));
      break;
    case 10:
      f = hi::basic_spectrum::MccYellowGreen();
      ::_ftprintf_s(stdout, _TEXT("mmc-yellow-green\n"));
      break;
    case 11:
      f = hi::basic_spectrum::MccOrangeYellow();
      ::_ftprintf_s(stdout, _TEXT("mmc-orange-yellow\n"));
      break;
    case 12:
      f = hi::basic_spectrum::MccBlue();
      ::_ftprintf_s(stdout, _TEXT("mmc-blue\n"));
      break;
    case 13:
      f = hi::basic_spectrum::MccGreen();
      ::_ftprintf_s(stdout, _TEXT("mmc-green\n"));
      break;
    case 14:
      f = hi::basic_spectrum::MccRed();
      ::_ftprintf_s(stdout, _TEXT("mmc-red\n"));
      break;
    case 15:
      f = hi::basic_spectrum::MccYellow();
      ::_ftprintf_s(stdout, _TEXT("mmc-yellow\n"));
      break;
    case 16:
      f = hi::basic_spectrum::MccMagenta();
      ::_ftprintf_s(stdout, _TEXT("mmc-magenta\n"));
      break;
    case 17:
      f = hi::basic_spectrum::MccCyan();
      ::_ftprintf_s(stdout, _TEXT("mmc-cyan\n"));
      break;
    case 18:
      f = hi::basic_spectrum::MccWhite();
      ::_ftprintf_s(stdout, _TEXT("mmc-white\n"));
      break;
    case 19:
      f = hi::basic_spectrum::MccNeutral80();
      ::_ftprintf_s(stdout, _TEXT("mmc-neutral.8\n"));
      break;
    case 20:
      f = hi::basic_spectrum::MccNeutral65();
      ::_ftprintf_s(stdout, _TEXT("mmc-neutral.65\n"));
      break;
    case 21:
      f = hi::basic_spectrum::MccNeutral50();
      ::_ftprintf_s(stdout, _TEXT("mmc-neutral.5\n"));
      break;
    case 22:
      f = hi::basic_spectrum::MccNeutral35();
      ::_ftprintf_s(stdout, _TEXT("mmc-neutral.35\n"));
      break;
    case 23:
      f = hi::basic_spectrum::MccBlack();
      ::_ftprintf_s(stdout, _TEXT("mmc-black\n"));
      break;
    case 24:
      f = 0;
      f[5] = 0.8;
      f[17] = 0.7;
      f[23] = 0.6;
      ::_ftprintf_s(stdout, _TEXT("mmc-spikes\n"));
      break;
    }

    hi::basic_spectrum b_fs[DIMENSION];
    InitializeFourierSeries(b_fs);
    hi::basic_spectrum b_cs[DIMENSION];
    InitializeCosineSeries(b_cs);
    hi::basic_spectrum b_lp[DIMENSION];
    InitializeLegendrePolynomial(b_lp);
    hi::basic_spectrum b_cp[DIMENSION];
    InitializeChebyshevPolynomial(b_cp[0], b_cp[1], b_cp[2], b_cp[3], b_cp[4],
                                  b_cp[5], b_cp[6]);
    hi::basic_spectrum b_sp[DIMENSION];
    InitializeSelectedPoints(b_sp);

    double c_fs[DIMENSION];
    double c_cs[DIMENSION];
    // double c_lp[DIMENSION];
    // double c_cp[DIMENSION];
    // double c_sp[DIMENSION];
    std::vector<double> c_d2;
    std::vector<double> c_d4;

    for (std::size_t k = 0; k < DIMENSION; ++k) {
      c_fs[k] = hi::dot(f, b_fs[k]);
      c_cs[k] = hi::dot(f, b_cs[k]);
      /*
      c_lp[k] = hi::dot(f, b_lp[k]);
      c_cp[k] = 0; // チェビシェフ多項式の場合
      for (std::size_t i = 0; i < hi::basic_spectrum::Size; ++i)
      {
        double const x = 2 * (i + 0.5) / hi::basic_spectrum::Size - 1;
        c_cp[k] += f[i] * b_cp[k][i] / std::sqrt(1-x*x);
      }
      c_sp[k] = hi::dot(f, b_sp[k]);
      */
    }
    // WTの場合
    FactrizeD2(f, c_d2);
    FactrizeD4(f, c_d4);
    /*
        ::_ftprintf_s(stdout, _TEXT("FS, %g, %g, %g, %g, %g, %g, %g\n"),
       c_fs[0], c_fs[1], c_fs[2], c_fs[3], c_fs[4], c_fs[5], c_fs[6]);
        ::_ftprintf_s(stdout, _TEXT("CS, %g, %g, %g, %g, %g, %g, %g\n"),
       c_cs[0], c_cs[1], c_cs[2], c_cs[3], c_cs[4], c_cs[5], c_cs[6]);
        ::_ftprintf_s(stdout, _TEXT("LP, %g, %g, %g, %g, %g, %g, %g\n"),
       c_lp[0], c_lp[1], c_lp[2], c_lp[3], c_lp[4], c_lp[5], c_lp[6]);
        ::_ftprintf_s(stdout, _TEXT("CP, %g, %g, %g, %g, %g, %g, %g\n"),
       c_cp[0], c_cp[1], c_cp[2], c_cp[3], c_cp[4], c_cp[5], c_cp[6]);
        ::_ftprintf_s(stdout, _TEXT("SP, %g, %g, %g, %g, %g, %g, %g\n"),
       c_sp[0], c_sp[1], c_sp[2], c_sp[3], c_sp[4], c_sp[5], c_sp[6]);
        ::_ftprintf_s(stdout, _TEXT("WT, %g, %g, %g, %g, %g, %g, %g\n"),
       c_d2[0], c_d2[1], c_d2[2], c_d2[3], c_d2[4], c_d2[5], c_d2[6]);
        ::_ftprintf_s(stdout, _TEXT("\n"));
    */
    hi::basic_spectrum g_fs(0);
    hi::basic_spectrum g_cs(0);
    hi::basic_spectrum g_d2(0);
    hi::basic_spectrum g_d4(0);
    // hi::basic_spectrum g_lp(0);
    // hi::basic_spectrum g_cp(0);
    // hi::basic_spectrum g_sp(0);

    /*
        typedef hi::basic_vector3<hi::basic_spectrum::value_type> vector_t;
        vector_t reference;
        hi::spd2xyz(f, reference); // f -> XYZ
        double const cScale = 1 / hi::length_squared(reference);

        // Sun (2000) の誤差尺度
        vector_t e_fs(0);
        vector_t e_cs(0);
        vector_t e_d2(0);
        vector_t e_d4(0);
      //vector_t e_lp(0);
      //vector_t e_cp(0);
      //vector_t e_sp(0);

        ::_ftprintf_s(stdout, _TEXT(", FS, CS, D2, D4\n"));
        for (int k = 0; k < DIMENSION; ++k)
        {
          // 再構築
          b_fs[k] *= c_fs[k]; g_fs += b_fs[k];
          b_cs[k] *= c_cs[k]; g_cs += b_cs[k];
        //b_lp[k] *= c_lp[k]; g_lp += b_lp[k];
        //b_cp[k] *= c_cp[k]; g_cp += b_cp[k];
        //b_sp[k] *= c_sp[k]; g_sp += b_sp[k];
          ReconstructD2(k+1, c_d2, g_d2); // WTの場合
          ReconstructD4(k+1, c_d4, g_d4); // D4の場合

          hi::spd2xyz(g_fs - f, e_fs);
          hi::spd2xyz(g_cs - f, e_cs);
          hi::spd2xyz(g_d2 - f, e_d2);
          hi::spd2xyz(g_d4 - f, e_d4);
        //hi::spd2xyz(g_lp - f, e_lp);
        //hi::spd2xyz(g_cp - f, e_cp);
        //hi::spd2xyz(g_sp - f, e_sp);

        ::_ftprintf_s(stdout, _TEXT("%d, %g, %g, %g, %g\n"), k+1
            , std::sqrt(hi::length_squared(e_fs) * cScale)
            , std::sqrt(hi::length_squared(e_cs) * cScale)
            , std::sqrt(hi::length_squared(e_d2) * cScale)
            , std::sqrt(hi::length_squared(e_d4) * cScale)
          //, std::sqrt(hi::length_squared(e_lp) * cScale)
          //, std::sqrt(hi::length_squared(e_cp) * cScale)
          //, std::sqrt(hi::length_squared(e_sp) * cScale)
           );
        }
        ::_ftprintf_s(stdout, _TEXT("\n"));
    /*/
    double e_fs_abs_l1[DIMENSION], e_fs_abs_l2[DIMENSION],
        e_fs_abs_lm[DIMENSION], e_fs_rel_l1[DIMENSION], e_fs_rel_l2[DIMENSION],
        e_fs_rel_lm[DIMENSION];
    double e_cs_abs_l1[DIMENSION], e_cs_abs_l2[DIMENSION],
        e_cs_abs_lm[DIMENSION], e_cs_rel_l1[DIMENSION], e_cs_rel_l2[DIMENSION],
        e_cs_rel_lm[DIMENSION];
    double e_d2_abs_l1[DIMENSION], e_d2_abs_l2[DIMENSION],
        e_d2_abs_lm[DIMENSION], e_d2_rel_l1[DIMENSION], e_d2_rel_l2[DIMENSION],
        e_d2_rel_lm[DIMENSION];
    double e_d4_abs_l1[DIMENSION], e_d4_abs_l2[DIMENSION],
        e_d4_abs_lm[DIMENSION], e_d4_rel_l1[DIMENSION], e_d4_rel_l2[DIMENSION],
        e_d4_rel_lm[DIMENSION];
    // double e_lp_abs_l1[DIMENSION], e_lp_abs_l2[DIMENSION],
    // e_lp_abs_lm[DIMENSION], e_lp_rel_l1[DIMENSION], e_lp_rel_l2[DIMENSION],
    // e_lp_rel_lm[DIMENSION];
    // double e_cp_abs_l1[DIMENSION], e_cp_abs_l2[DIMENSION],
    // e_cp_abs_lm[DIMENSION], e_cp_rel_l1[DIMENSION], e_cp_rel_l2[DIMENSION],
    // e_cp_rel_lm[DIMENSION];
    // double e_sp_abs_l1[DIMENSION], e_sp_abs_l2[DIMENSION],
    // e_sp_abs_lm[DIMENSION], e_sp_rel_l1[DIMENSION], e_sp_rel_l2[DIMENSION],
    // e_sp_rel_lm[DIMENSION];

    for (int k = 0; k < DIMENSION; ++k) {
      b_fs[k] *= c_fs[k];
      g_fs += b_fs[k];
      b_cs[k] *= c_cs[k];
      g_cs += b_cs[k];
      ReconstructD2(k + 1, c_d2, g_d2);
      ReconstructD4(k + 1, c_d4, g_d4);
      // b_lp[k] *= c_lp[k]; g_lp += b_lp[k];
      // b_cp[k] *= c_cp[k]; g_cp += b_cp[k];
      // b_sp[k] *= c_sp[k]; g_sp += b_sp[k];

      e_fs_abs_l1[k] = e_fs_abs_l2[k] = e_fs_abs_lm[k] = e_fs_rel_l1[k] =
          e_fs_rel_l2[k] = e_fs_rel_lm[k] = 0;
      e_cs_abs_l1[k] = e_cs_abs_l2[k] = e_cs_abs_lm[k] = e_cs_rel_l1[k] =
          e_cs_rel_l2[k] = e_cs_rel_lm[k] = 0;
      e_d2_abs_l1[k] = e_d2_abs_l2[k] = e_d2_abs_lm[k] = e_d2_rel_l1[k] =
          e_d2_rel_l2[k] = e_d2_rel_lm[k] = 0;
      e_d4_abs_l1[k] = e_d4_abs_l2[k] = e_d4_abs_lm[k] = e_d4_rel_l1[k] =
          e_d4_rel_l2[k] = e_d4_rel_lm[k] = 0;
      // e_lp_abs_l1[k] = e_lp_abs_l2[k] = e_lp_abs_lm[k] = e_lp_rel_l1[k] =
      // e_lp_rel_l2[k] = e_lp_rel_lm[k] = 0;
      // e_cp_abs_l1[k] = e_cp_abs_l2[k] = e_cp_abs_lm[k] = e_cp_rel_l1[k] =
      // e_cp_rel_l2[k] = e_cp_rel_lm[k] = 0;
      // e_sp_abs_l1[k] = e_sp_abs_l2[k] = e_sp_abs_lm[k] = e_sp_rel_l1[k] =
      // e_sp_rel_l2[k] = e_sp_rel_lm[k] = 0;

      ::_ftprintf_s(
          stdout,
          _TEXT("approx ~ %d,  f")
              _TEXT(",    g_fs   ,    g_cs   ,    g_d2   ,    g_d4")
                  _TEXT(", |f-g_fs|  , |f-g_cs|  , |f-g_d2|  , |f-g_d4|") _TEXT(
                      ", |f-g_fs|/f, |f-g_cs|/f, |f-g_d2|/f, |f-g_d4|/f\n"),
          k);
      for (std::size_t i = 0; i < hi::basic_spectrum::Size; ++i) {
        double const e_fs_abs = std::abs(f[i] - g_fs[i]);
        double const e_cs_abs = std::abs(f[i] - g_cs[i]);
        double const e_d2_abs = std::abs(f[i] - g_d2[i]);
        double const e_d4_abs = std::abs(f[i] - g_d4[i]);
        // double const e_lp_abs = std::abs(f[i] - g_lp[i]);
        // double const e_cp_abs = std::abs(f[i] - g_cp[i]);
        // double const e_sp_abs = std::abs(f[i] - g_sp[i]);

        double const e_fs_rel = (f[i] > 0) ? e_fs_abs / f[i] : 0;
        double const e_cs_rel = (f[i] > 0) ? e_cs_abs / f[i] : 0;
        double const e_d2_rel = (f[i] > 0) ? e_d2_abs / f[i] : 0;
        double const e_d4_rel = (f[i] > 0) ? e_d4_abs / f[i] : 0;
        // double const e_lp_rel = (f[i]>0) ? e_lp_abs / f[i] : 0;
        // double const e_cp_rel = (f[i]>0) ? e_cp_abs / f[i] : 0;
        // double const e_sp_rel = (f[i]>0) ? e_sp_abs / f[i] : 0;

        e_fs_abs_l1[k] += e_fs_abs;
        e_fs_abs_l2[k] += hi::square_of(e_fs_abs);
        e_fs_abs_lm[k] = std::max(e_fs_abs_lm[k], e_fs_abs);
        e_fs_rel_l1[k] += e_fs_rel;
        e_fs_rel_l2[k] += hi::square_of(e_fs_rel);
        e_fs_rel_lm[k] = std::max(e_fs_rel_lm[k], e_fs_rel);
        e_cs_abs_l1[k] += e_cs_abs;
        e_cs_abs_l2[k] += hi::square_of(e_cs_abs);
        e_cs_abs_lm[k] = std::max(e_cs_abs_lm[k], e_cs_abs);
        e_cs_rel_l1[k] += e_cs_rel;
        e_cs_rel_l2[k] += hi::square_of(e_cs_rel);
        e_cs_rel_lm[k] = std::max(e_cs_rel_lm[k], e_cs_rel);
        e_d2_abs_l1[k] += e_d2_abs;
        e_d2_abs_l2[k] += hi::square_of(e_d2_abs);
        e_d2_abs_lm[k] = std::max(e_d2_abs_lm[k], e_d2_abs);
        e_d2_rel_l1[k] += e_d2_rel;
        e_d2_rel_l2[k] += hi::square_of(e_d2_rel);
        e_d2_rel_lm[k] = std::max(e_d2_rel_lm[k], e_d2_rel);
        e_d4_abs_l1[k] += e_d4_abs;
        e_d4_abs_l2[k] += hi::square_of(e_d4_abs);
        e_d4_abs_lm[k] = std::max(e_d4_abs_lm[k], e_d4_abs);
        e_d4_rel_l1[k] += e_d4_rel;
        e_d4_rel_l2[k] += hi::square_of(e_d4_rel);
        e_d4_rel_lm[k] = std::max(e_d4_rel_lm[k], e_d4_rel);
        // e_lp_abs_l1[k] += e_lp_abs; e_lp_abs_l2[k] +=
        // hi::square_of(e_lp_abs); e_lp_abs_lm[k] = std::max(e_lp_abs_lm[k],
        // e_lp_abs);
        // e_lp_rel_l1[k] += e_lp_rel; e_lp_rel_l2[k] +=
        // hi::square_of(e_lp_rel); e_lp_rel_lm[k] = std::max(e_lp_rel_lm[k],
        // e_lp_rel);
        // e_cp_abs_l1[k] += e_cp_abs; e_cp_abs_l2[k] +=
        // hi::square_of(e_cp_abs); e_cp_abs_lm[k] = std::max(e_cp_abs_lm[k],
        // e_cp_abs);
        // e_cp_rel_l1[k] += e_cp_rel; e_cp_rel_l2[k] +=
        // hi::square_of(e_cp_rel); e_cp_rel_lm[k] = std::max(e_cp_rel_lm[k],
        // e_cp_rel);
        // e_sp_abs_l1[k] += e_sp_abs; e_sp_abs_l2[k] +=
        // hi::square_of(e_sp_abs); e_sp_abs_lm[k] = std::max(e_sp_abs_lm[k],
        // e_sp_abs);
        // e_sp_rel_l1[k] += e_sp_rel; e_sp_rel_l2[k] +=
        // hi::square_of(e_sp_rel); e_sp_rel_lm[k] = std::max(e_sp_rel_lm[k],
        // e_sp_rel);

        ::_ftprintf_s(stdout,
                      _TEXT("%d, %g") _TEXT(", %g, %g, %g, %g")
                          _TEXT(", %g, %g, %g, %g") _TEXT(", %g, %g, %g, %g\n"),
                      400 + i * 10, f[i], g_fs[i], g_cs[i], g_d2[i], g_d4[i],
                      e_fs_abs, e_cs_abs, e_d2_abs, e_d4_abs, e_fs_rel,
                      e_cs_rel, e_d2_rel, e_d4_rel);
      }

      e_fs_abs_l2[k] = std::sqrt(e_fs_abs_l2[k]);
      e_fs_rel_l2[k] = std::sqrt(e_fs_rel_l2[k]);
      e_cs_abs_l2[k] = std::sqrt(e_cs_abs_l2[k]);
      e_cs_rel_l2[k] = std::sqrt(e_cs_rel_l2[k]);
      e_d2_abs_l2[k] = std::sqrt(e_d2_abs_l2[k]);
      e_d2_rel_l2[k] = std::sqrt(e_d2_rel_l2[k]);
      e_d4_abs_l2[k] = std::sqrt(e_d4_abs_l2[k]);
      e_d4_rel_l2[k] = std::sqrt(e_d4_rel_l2[k]);
      // e_lp_abs_l2[k] = std::sqrt(e_lp_abs_l2[k]); e_lp_rel_l2[k] =
      // std::sqrt(e_lp_rel_l2[k]);
      // e_cp_abs_l2[k] = std::sqrt(e_cp_abs_l2[k]); e_cp_rel_l2[k] =
      // std::sqrt(e_cp_rel_l2[k]);
      // e_sp_abs_l2[k] = std::sqrt(e_sp_abs_l2[k]); e_sp_rel_l2[k] =
      // std::sqrt(e_sp_rel_l2[k]);
    }
    ::_ftprintf_s(stdout, _TEXT("\n"));

    ::_ftprintf_s(
        stdout,
        _TEXT(
            ", e_fs_abs_l1, e_cs_abs_l1, e_lp_abs_l1, e_cp_abs_l1, "
            "e_sp_abs_l1, e_d2_abs_l1")
            _TEXT(
                ", e_fs_abs_l2, e_cs_abs_l2, e_lp_abs_l2, e_cp_abs_l2, "
                "e_sp_abs_l2, e_d2_abs_l2")
                _TEXT(
                    ", e_fs_abs_lm, e_cs_abs_lm, e_lp_abs_lm, e_cp_abs_lm, "
                    "e_sp_abs_lm, e_d2_abs_lm")
                    _TEXT(
                        ", e_fs_rel_l1, e_cs_rel_l1, e_lp_rel_l1, "
                        "e_cp_rel_l1, e_sp_rel_l1, e_d2_rel_l1")
                        _TEXT(
                            ", e_fs_rel_l2, e_cs_rel_l2, e_lp_rel_l2, "
                            "e_cp_rel_l2, e_sp_rel_l2, e_d2_rel_l2")
                            _TEXT(
                                ", e_fs_rel_lm, e_cs_rel_lm, e_lp_rel_lm, "
                                "e_cp_rel_lm, e_sp_rel_lm, e_d2_rel_lm\n"));
    for (std::size_t k = 0; k < DIMENSION; ++k) {
      ::_ftprintf_s(
          stdout,
          _TEXT("%d") _TEXT(", %g, %g, %g, %g") _TEXT(", %g, %g, %g, %g")
              _TEXT(", %g, %g, %g, %g") _TEXT(", %g, %g, %g, %g")
                  _TEXT(", %g, %g, %g, %g") _TEXT(", %g, %g, %g, %g\n"),
          k + 1, e_fs_abs_l1[k], e_cs_abs_l1[k], e_d2_abs_l1[k], e_d4_abs_l1[k],
          e_fs_abs_l2[k], e_cs_abs_l2[k], e_d2_abs_l2[k], e_d4_abs_l2[k],
          e_fs_abs_lm[k], e_cs_abs_lm[k], e_d2_abs_lm[k], e_d4_abs_lm[k],
          e_fs_rel_l1[k], e_cs_rel_l1[k], e_d2_rel_l1[k], e_d4_rel_l1[k],
          e_fs_rel_l2[k], e_cs_rel_l2[k], e_d2_rel_l2[k], e_d4_rel_l2[k],
          e_fs_rel_lm[k], e_cs_rel_lm[k], e_d2_rel_lm[k], e_d4_rel_lm[k]);
    }
    ::_ftprintf_s(stdout, _TEXT("\n"));
    //*/
  }
  return 0;
}
#endif
