// pvpm_test.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

#include "stdafx.h"

// utility function
double sgn(double x) { return (x < 0) ? -1 : (x > 0) ? 1 : 0; }

// test function 1 (x in [0,1] and Y = 0.75)
double step_function(double x) { return 0.75 + 0.25 * sgn(x - 0.5); }

// test function 2 (x in [0,1] and Y = 1/3)
double quadratic_function(double x) { return x * x; }

// test function 3 (x in [0,1] and Y = pi)
double pi_function(double x) { return 4 * std::sqrt(1 - x * x); }

// test function 4 (x in [0,1] and Y = 1/2)
double linear_function(double x) { return x; }

typedef double float_t;
typedef double (*test_function)(double x);
typedef hi::basic_vector3<float_t> float3_t;

// 1D integration tests
int _tmain(int argc, _TCHAR *argv[]) {
  _TCHAR *name[] = {_T("step function"), _T("quadratic function"),
                    _T("pi function"), _T("linear function")};
  test_function fn[] = {step_function, quadratic_function, pi_function,
                        linear_function};
  double const U[] = {1, 1, 4, 1};                  // 上限
  double const Y[] = {0.75, 1.0 / 3.0, M_PI, 0.5};  // 積分値

  // サンプリング点の数
  static std::size_t const M = 1;

  // サンプリング点
  float_t x_samp[M];

  // サンプリング点で共通のパラメータ
  float_t pow;  // 積分値
  float_t rad;  // 有効半径
  float_t cnt;  // カーネル関数の寄与

  float_t powU;
  float_t radU;
  float_t cntU;

  // フォトンの数
  static std::size_t const N = 10;
  static std::size_t const K = static_cast<std::size_t>(std::sqrt(double(N)));

  // フォトンのパラメータ
  float3_t photons[N];
  float_t distance[N];

  // 収束パラメータ
  static float_t const Alpha = 0.99;

  std::mt19937_64 generator;
  for (std::size_t k = 0; k < 4; ++k) {
    ::_ftprintf_s(stdout, _T("%s\n"), name[k]);

    // 初期化
    pow = 0;
    rad = 1;
    cnt = 0;

    powU = 0;
    radU = 1;
    cntU = 0;

    for (std::size_t i = 1; i <= 10000; ++i) {
// サンプリング点の初期化
#if 1
      for (std::size_t m = 0; m < M; ++m) {
        x_samp[m] = (m + generator.next<float_t>()) / M;
      }
#else
      for (std::size_t m = 0; m < M; ++m) {
        x_samp[m] = m / double(M - 1);
      }
#endif

      // フォトンの射出
      for (std::size_t n = 0; n < N; ++n) {
        // 棄却サンプリングによるフォトンのサンプリング
        double x;
        double y;
        do {
          x = generator.next<float_t>();
          y = fn[k](x);
        } while (generator.next<float_t>() * U[k] > y);

        photons[n][0] = x;  // x
        photons[n][1] = y;  // f(x)
      }

      // フォトンの半径を決定
      for (std::size_t n = 0; n < N; ++n) {
        for (std::size_t m = 0; m < N; ++m) {
          distance[m] = std::abs(photons[m][0] - photons[n][0]);
        }
        std::sort(distance, distance + N);
        photons[n][2] = distance[K] * 10;
      }

// フォトンの収集
#if 1
      for (std::size_t m = 0; m < M; ++m) {  // サンプリング点のループ
        for (std::size_t n = 0; n < N; ++n) {  // フォトンのループ

          double const r =
              std::abs(x_samp[m] - photons[n][0]);  // フォトンまでの距離

          // 1D Tent Kernel: (1/h)*(1-r/h)
          {
            double const h = photons[n][2] * rad;  // フォトンの影響半径
            if (r >= h) {
              goto next_sampler;  // no contribution
            }

            float_t const wgt = 1 - r / h;
            if (wgt <= 0) {
              goto next_sampler;  // no contribution
            }
            float_t const acm = wgt * Alpha;
            float_t const d = (cnt + acm) / (cnt + wgt);
            if (d > 1) {
              goto next_sampler;  // no contribution
            }

            float_t const s = 1 - std::sqrt(1 - d);
            if (s <= 0) {
              goto next_sampler;
            }

            // 更新
            pow += (wgt / h);
            pow *= d;
            rad *= s;
            cnt += acm;
          }
        next_sampler : {
          double const hU = photons[n][2] * radU;  // フォトンの影響半径
          if (r >= hU) {
            continue;  // no contribution
          }

          float_t const wgtU = 0.5;
          float_t const acmU = wgtU * Alpha;
          float_t const dU = (cntU + acmU) / (cntU + wgtU);
          float_t const sU = dU;

          powU += (wgtU / hU);
          powU *= dU;
          radU *= sU;
          cntU += acmU;
        }
        }
      }
      {
        double const I = pow / (i * N * M);
        double const IU = powU / (i * N * M);
        //::_ftprintf_s(stdout, _T("%g\t%g\n"), std::abs(I-1), std::abs(IU-1));
        //::_ftprintf_s(stdout, _T("%g\n"), std::abs(I-1));
        ::_ftprintf_s(stdout, _T("%g\n"), std::abs(IU - 1));
      }
#else
      for (std::size_t m = 0; m < M; ++m) {  // サンプリング点のループ
        double y1 = 0;
        double y2 = 0;
        for (std::size_t n = 0; n < N; ++n) {    // フォトンのループ
          double const h = photons[n][2] * rad;  // フォトンの影響半径
          double const r =
              std::abs(x_samp[m] - photons[n][0]);  // フォトンまでの距離
          if (r >= h) {
            continue;  // no contribution
          }

          // 更新
          y1 += ((0.5) / h) * Y[k];        // ユニフォームカーネル
          y2 += ((1 - r / h) / h) * Y[k];  // テントカーネル
        }
        y1 /= i * N;
        y2 /= i * N;

        double const y0 = fn[k](x_samp[m]);
        ::_ftprintf_s(stdout, _T("%d\t%g\t%g\t%g\n"), m, y0, y1, y2);
      }
#endif
    }
  }

  // relative error

  return 0;
}
