
#include "Evaluator.hpp"
#include "SamplingState.hpp"
#include "core/shell.hpp"
#include "core/Scene.hpp"
#include "core/Camera.hpp"
#include "core/Film.hpp"

namespace tgir {
void RELT::run() {
  std::tcerr << _TEXT("** Replica-Exchange Light Transport **") << std::endl;

  typedef hi::basic_vector4<tgir::Real> RealVectorK;
  typedef hi::basic_vector3<tgir::Real> RealVectorKm1;
  typedef hi::basic_vector4<hi::uint> UintVectorK;
  typedef hi::basic_vector3<hi::uint> UintVectorKm1;
  static std::size_t const K = RealVectorK::N;
  static std::size_t const N =
      TGIR_CONFIG_kSampleCount / (TGIR_CONFIG_kSystemCount * K);
#if 0
    static RealVectorKm1 const
      vMutationWidth(
        hi::rcp<tgir::Real>(32),
        hi::rcp<tgir::Real>(64),
        hi::rcp<tgir::Real>(96) ); // 変異幅
#else
  static RealVectorKm1 const vMutationWidth(
      hi::rcp<tgir::Real>(64), hi::rcp<tgir::Real>(96),
      hi::rcp<tgir::Real>(128));  // 変異幅
#endif

  tgir::Scene const &scene = tgir::Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  tgir::Film m(W, H);

  std::size_t const M = m.GetCount();

  tgir::Film I[K];
  for (std::size_t k = 0; k < K; ++k) {
    I[k].Resize(W, H);
  }

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  for (std::size_t k = 0; k < K; ++k) {
    ::_stprintf_s(filename, MAX_PATH, _TEXT("%s") _TEXT("%u/"),
                  tgir::GetOutdirPath().c_str(), k);
    hi::mkdir(filename);
  }

  struct Replica {
    typedef tgir::SamplingState<RealVectorK> SamplingState;

    std::size_t current;
    SamplingState state[2];

    Replica() : current(0) {}

    // a set of the pdfs
    inline static void PDFs(tgir::Path const &path, RealVectorK &g) {
#if 1
      g[0] = path.PDFs[2];
      g[1] = path.PDFs[1];
      g[2] = path.PDFs[0];
      g[3] = 1;
#else
      tgir::Real const brightness = 1 - std::exp(-path.PDFs[0]);
      g[0] = path.PDFs[0];  // liner space value
      g[1] = brightness;    // tonemapped value
      g[2] = std::sqrt(
          brightness);  // tonemapped value & gamma 2.0 corrected value
      g[3] = 1;         // uniform distribution
#endif
    }

    // return (a/b) which limited in [0,1]
    inline static tgir::Real Divide(tgir::Real const x, tgir::Real const y) {
      return (y < x) ? (y / x) : ((x > 0) ? 1 : 0);
    }
  };

  struct System {
    Replica replicas[K];
    std::size_t indices[K];

    System() {
      for (std::size_t k = 0; k < K; ++k) {
        indices[k] = k;
      }
    }

  } copies[TGIR_CONFIG_kSystemCount];

  struct ThreadData {
    ThreadData()
        : path(TGIR_CONFIG_kMaxPathLength),
          b(0),
          acceptanceRatio(0),
          exchangeRatio(0),
          exchangeCount(0) {}

    hi::mt19937 random;
    tgir::Path path;
    RealVectorK b;
    RealVectorK acceptanceRatio;
    RealVectorKm1 exchangeRatio;
    UintVectorKm1 exchangeCount;
  };

  std::vector<ThreadData> threadData(omp_get_max_threads());

  // 積分値の逆数
  RealVectorK b_hat(0);
#if !defined(CONFIG_MULTIPLE_INTEGRAION)
  RealVectorK b(0);
#endif

  // 計測開始
  int const renderingTime = tgir::GetRenderingTime();
  int const intervalTime = tgir::GetIntervalTime();
  int saveTime = intervalTime;
  std::clock_t beginTime = std::clock();

  //
  // Stage1: Initialization.
  //
  std::tcerr << _TEXT("Stage1: Initialization.") << std::endl;
#if defined(CONFIG_MULTIPLE_INTEGRAION) && 1
  //
  // Population Annealing
  //
  {
    std::size_t const SEEDS_COUNT = CONFIG_SEEDS * 2 / K;
    std::size_t current = 0;
    std::vector<Replica::SamplingState> seeds[2] = {
        std::vector<Replica::SamplingState>(SEEDS_COUNT),
        std::vector<Replica::SamplingState>(SEEDS_COUNT),
    };

    // cumulative distribution function
    std::vector<tgir::Real> cdf(SEEDS_COUNT + 1);
    cdf[0] = 0;

#define CUMULATIVE_DISTRIBUTION_FUNCTION(k)  \
  {                                          \
    for (int i = 1; i <= SEEDS_COUNT; ++i) { \
      if (cdf[i] > 0) {                      \
        ++n[k];                              \
      }                                      \
      cdf[i] += cdf[i - 1];                  \
    }                                        \
  }
#if 0
#define POPULATION_RESAMPLING()                                             \
  {                                                                         \
    for (std::size_t i = 0; i < SEEDS_COUNT; ++i) {                         \
      ThreadData &data = threadData[omp_get_thread_num()];                  \
      tgir::Real const xi =                                                 \
          cdf.back() * (i + data.random.next<tgir::Real>()) / SEEDS_COUNT;  \
      std::vector<tgir::Real>::const_iterator const it =                    \
          std::upper_bound(cdf.Begin(), cdf.End(), xi);                     \
      std::size_t const j = static_cast<std::size_t>(it - cdf.Begin() - 1); \
      seeds[current ^ 1][i] = seeds[current][j];                            \
    }                                                                       \
    current ^= 1;                                                           \
  }
#else
#define POPULATION_RESAMPLING()                                            \
  {                                                                        \
    std::size_t j = 0;                                                     \
    for (std::size_t i = 0; i < SEEDS_COUNT; ++i) {                        \
      ThreadData &data = threadData[omp_get_thread_num()];                 \
      tgir::Real const xi =                                                \
          cdf.back() * (i + data.random.next<tgir::Real>()) / SEEDS_COUNT; \
      while ((cdf[j] > xi) || (xi >= cdf[j + 1])) ++j;                     \
      seeds[current ^ 1][i] = seeds[current][j];                           \
    }                                                                      \
    current ^= 1;                                                          \
  }
#endif
#define CHOOSE_SAMPLES(k)                                                     \
  for (std::size_t c = 0; c < TGIR_CONFIG_kSystemCount; ++c) {                \
    Replica &replica = copies[c].replicas[copies[c].indices[(k)]];            \
    replica.state[replica.current] =                                          \
        seeds[current][seeds[current].size() * c / TGIR_CONFIG_kSystemCount]; \
  }
#define METROPOLIS_SAMPLING(k)                                                \
  for (int i = 0; i < SEEDS_COUNT; ++i) {                                     \
    ThreadData &data = threadData[omp_get_thread_num()];                      \
    Replica::SamplingState &oldSample = seeds[current][i];                    \
    Replica::SamplingState &newSample = seeds[current ^ 1][i];                \
    newSample.sample.Mutate(oldSample.sample, vMutationWidth[k],              \
                            data.random); /* mutate a path */                 \
    scene.Evalute(newSample.sample, data.path,                                \
                  newSample.values); /* evalute the path */                   \
    Replica::PDFs(data.path, newSample.q);                                    \
    tgir::Real const a = Replica::Divide(                                     \
        oldSample.q[k], newSample.q[k]); /* compute acceptance probability */ \
    data.b[k - 1] += (oldSample.q[k] > 0)                                     \
                         ? (1 - a) * oldSample.q[k - 1] / oldSample.q[k]      \
                         : 0;                                                 \
    data.b[k - 1] +=                                                          \
        (newSample.q[k] > 0) ? (a)*newSample.q[k - 1] / newSample.q[k] : 0;   \
    if (data.random.next<tgir::Real>() >= a) {                                \
      newSample.swap(oldSample);                                              \
    } /* rejected new sample */                                               \
    cdf[i + 1] =                                                              \
        (newSample.q[k] > 0) ? newSample.q[k - 1] / newSample.q[k] : 0;       \
  }                                                                           \
  current ^= 1;

    UintVectorK n(0);
    n[K - 1] = SEEDS_COUNT;

    // 初期サンプルの生成
    {
      int k = K - 1;
#pragma omp parallel for
      for (int i = 0; i < SEEDS_COUNT; ++i) {
        ThreadData &data = threadData[omp_get_thread_num()];
        Replica::SamplingState &newSample = seeds[current][i];
        newSample.sample.Init(data.random);
        scene.Evalute(newSample.sample, data.path, newSample.values);
        Replica::PDFs(data.path, newSample.q);
        data.b[k] += newSample.q[k];
        data.b[k - 1] += newSample.q[k - 1] / newSample.q[k];
        cdf[i + 1] = newSample.q[k - 1] / newSample.q[k];
      }
      CUMULATIVE_DISTRIBUTION_FUNCTION(k - 1);
      if (cdf.back() <= 0) {
        std::tcerr << _TEXT("初期サンプルの生成に失敗しました．") << std::endl;
        tgir::FinishRendering();
        return;  // レンダリングする必要なし
      }
      POPULATION_RESAMPLING();
      CHOOSE_SAMPLES(k - 1);
    }

    for (int k = K - 2; k > 0; --k) {
#pragma omp parallel for
      METROPOLIS_SAMPLING(k);
      CUMULATIVE_DISTRIBUTION_FUNCTION(k - 1);
      if (cdf.back() <= 0) {
        std::tcerr << _TEXT("初期サンプルの生成に失敗しました．") << std::endl;
        tgir::FinishRendering();
        return;  // レンダリングする必要なし
      }
      POPULATION_RESAMPLING();
      CHOOSE_SAMPLES(k - 1);
    }

    /// Metropolis sampling

    // 仮の正規化係数の計算
      for
        each(ThreadData const &data in threadData) { b_hat += data.b; }
      for (int k = 0; k < K; ++k) {
        b_hat[k] = std::log(b_hat[k]);
      }
      b_hat -= b_hat[K - 1];
      for (int k = K - 2; k >= 0; --k) {
        b_hat[k] += b_hat[k + 1];
      }
      for (int k = 0; k < K; ++k) {
        b_hat[k] = std::exp(b_hat[k]);
      }

      // 正規化定数の表示
      ::_ftprintf_s(
          stdout, _TEXT("  b_hat={%g, %g, %g, %g} n={%6u, %6u, %6u, %6u}\n"),
          b_hat[0], b_hat[1], b_hat[2], b_hat[3], n[0], n[1], n[2], n[3]);

      if ((b_hat[0] <= 0) || (b_hat[1] <= 0) || (b_hat[2] <= 0) ||
          (b_hat[3] <= 0)) {
        tgir::FinishRendering();
        return;  // レンダリングする必要なし
      }

      b_hat = hi::rcp(b_hat);  // 逆数を求める
#undef METROPOLIS_SAMPLING
#undef POPULATION_RESAMPLING
#undef CUMULATIVE_DISTRIBUTION_FUNCTION
  }
#else
  {
    std::vector<Replica::SamplingState> seeds(CONFIG_SEEDS);
    std::vector<tgir::Real> cdf[K];
    for (std::size_t k = 0; k < K; ++k) {
      cdf[k].Resize(CONFIG_SEEDS + 1);
    }
    {
// generate path seeds
#pragma omp parallel for
      for (int i = 0; i < CONFIG_SEEDS; ++i) {
        ThreadData &data = threadData[omp_get_thread_num()];
        seeds[i].sample.Init(data.random);
        Replica::PDFs(
            scene.Evalute(seeds[i].sample, data.path, seeds[i].values),
            data.path, seeds[i].q);
        for (std::size_t k = 0; k < K; ++k) {
          cdf[k][i + 1] = seeds[i].q[k];
        }
      }

      UintVectorK n(0);  // 生成できたサンプル数
      for (std::size_t k = 0; k < K; ++k) {
        cdf[k][0] = 0;
        for (int i = 1; i <= CONFIG_SEEDS; ++i) {
          if (cdf[k][i] > 0) {
            ++n[k];
          }
          cdf[k][i] += cdf[k][i - 1];
        }
#if !defined(CONFIG_MULTIPLE_INTEGRAION)
        b[k] = cdf[k].back();
#endif
        b_hat[k] = cdf[k].back() / CONFIG_SEEDS;
      }

      // 正規化定数の表示
      ::_ftprintf_s(
          stderr, _TEXT("  b_hat={%g, %g, %g, %g} n={%6u, %6u, %6u, %6u}\n"),
          b_hat[0], b_hat[1], b_hat[2], b_hat[3], n[0], n[1], n[2], n[3]);

      if ((b_hat[0] <= 0) || (b_hat[1] <= 0) || (b_hat[2] <= 0) ||
          (b_hat[3] <= 0)) {
        tgir::FinishRendering();
        return;  // レンダリングする必要なし
      }

      b_hat = hi::rcp(b_hat);  // 逆数を求める
    }

    // choose SamplingState
    for (std::size_t c = 0; c < TGIR_CONFIG_kSystemCount; ++c) {
      ThreadData &data = threadData[omp_get_thread_num()];

      Replica *const replica = copies[c].replicas;
      std::size_t *const indices = copies[c].current;

      for (std::size_t k = 0; k < K; ++k) {
        tgir::Real const xi = cdf[k].back() * data.random.next<tgir::Real>();
        std::vector<tgir::Real>::const_iterator const it =
            std::upper_bound(cdf[k].Begin(), cdf[k].End(), xi);
        Replica &replica = replica[indices[k]];  // k-th replica
        replica.state[replica.current] =
            seeds[static_cast<std::size_t>(it - cdf[k].Begin() - 1)];
      }
    }
  }
#endif

  {
    std::clock_t const currentTime = std::clock() - beginTime;
    std::clock_t const time =
        static_cast<std::clock_t>(currentTime / double(CLOCKS_PER_SEC) + 0.5);
    ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒\n"), time / (60 * 60),
                  (time / 60) % 60, time % 60);
    beginTime = std::clock();
  }

  //
  // Stage2: Computation.
  //
  std::tcerr << _TEXT("Stage2: Computation.") << std::endl;

  for (std::size_t i = 1;; ++i) {
#pragma omp parallel for
    for (int c = 0; c < TGIR_CONFIG_kSystemCount; ++c) {
      // キャッシュ
      ThreadData &data = threadData[omp_get_thread_num()];
      Replica *const replicas = copies[c].replicas;
      std::size_t *const indices = copies[c].indices;

      for (std::size_t n = 0; n < N; ++n) {
        //
        // Metropolis Process
        //
        for (std::size_t k = 0; k < (K - 1); ++k) {
          Replica &replica = replicas[indices[k]];
          Replica::SamplingState &oldSample = replica.state[replica.current];
          Replica::SamplingState &newSample =
              replica.state[replica.current ^ 1];

          // 変異
          for (std::size_t nTry = 0; nTry < 8; ++nTry) {
            newSample.sample.Mutate(oldSample.sample, vMutationWidth[k],
                                    data.random);  // mutate a path
            scene.Evalute(newSample.sample, data.path,
                          newSample.values);  // evalute the path
            Replica::PDFs(data.path, newSample.q);
            if (newSample.q[k] > 0) {
              break;
            }
          }

          tgir::Real const a = Replica::Divide(
              oldSample.q[k], newSample.q[k]);  // acceptance probability
          data.acceptanceRatio[k] += a;
#if defined(CONFIG_APPROX_INITIAL_B)
          I[k].Deposite(
              oldSample.values,
              (1 - a) *
                  hi::rcp(hi::dot(oldSample.q, b_hat)));  // deposite old sample
          I[k].Deposite(
              newSample.values,
              (a)*hi::rcp(hi::dot(newSample.q, b_hat)));  // deposite new sample
#else
          I[k].Deposite(
              oldSample.values,
              (1 - a) * hi::rcp(hi::sum(oldSample.q)));  // deposite old sample
          I[k].Deposite(
              newSample.values,
              (a)*hi::rcp(hi::sum(newSample.q)));  // deposite new sample
#endif
#if defined(CONFIG_MULTIPLE_INTEGRAION)
          if (k) {
            data.b[k - 1] += (oldSample.q[k] > 0)
                                 ? (1 - a) * oldSample.q[k - 1] / oldSample.q[k]
                                 : 0;
            data.b[k - 1] += (newSample.q[k] > 0)
                                 ? (a)*newSample.q[k - 1] / newSample.q[k]
                                 : 0;
          }
#endif
          if (data.random.next<tgir::Real>() < a) {
            replica.current ^= 1;
          }  // accept
        }
        {
          std::size_t const k = K - 1;
          Replica &replica = replicas[indices[k]];
          Replica::SamplingState &newSample =
              replica.state[replica.current ^ 1];

          // 変異
          newSample.sample.Init(data.random);  // mutate a path
          scene.Evalute(newSample.sample, data.path,
                        newSample.values);  // evalute the path
          Replica::PDFs(data.path, newSample.q);

#if defined(CONFIG_MULTIPLE_INTEGRAION)
          {
            data.b[k - 1] +=
                (newSample.q[k] > 0) ? newSample.q[k - 1] / newSample.q[k] : 0;
            data.b[k] += newSample.q[k];
          }
#else
          data.b += newSample.q;
#endif
          ++data.acceptanceRatio[k];
#if defined(CONFIG_APPROX_INITIAL_B)
          I[k].Deposite(newSample.values, /*-------*/ hi::rcp(hi::dot(
                            newSample.q, b_hat)));  // deposite new sample
#else
          I[k].Deposite(newSample.values, /*-------*/ hi::rcp(
                            hi::sum(newSample.q)));  // deposite new sample
#endif
          /*----------------------------------*/ {
            replica.current ^= 1;
          }  // accept
        }
        //
        // Exchagne Process
        //
        for (std::size_t k = 0; k < (K - 1); ++k) {
          // std::size_t const k = n % (K-1);
          Replica const &replica0 = replicas[indices[k + 0]];
          Replica const &replica1 = replicas[indices[k + 1]];
          Replica::SamplingState const &state0 =
              replica0.state[replica0.current];  // (k+0)^th replica
          Replica::SamplingState const &state1 =
              replica1.state[replica1.current];  // (k+1)^th replica
          tgir::Real const e = Replica::Divide(state0.q[k] * state1.q[k + 1],
                                               state0.q[k + 1] * state1.q[k]);
          data.exchangeRatio[k] += e;
          ++data.exchangeCount[k];
          if (data.random.next<tgir::Real>() < e) {
            std::swap(indices[k], indices[k + 1]);
          }  // exchagne
        }
      }
    }

    // 定期的な保存
    {
      std::clock_t const currentTime = std::clock() - beginTime;
      std::clock_t const time =
          static_cast<std::clock_t>(currentTime / double(CLOCKS_PER_SEC) + 0.5);
      hi::ulong const mutation = hi::ulong(i) * N * TGIR_CONFIG_kSystemCount;
      tgir::Real const fLogMutation = std::log(tgir::Real(mutation));
      tgir::Real const fMutationPerPixel =
          std::exp(fLogMutation - std::log(tgir::Real(M)));

#if defined(CONFIG_MULTIPLE_INTEGRAION)
      RealVectorK pdf(0);
#else
      RealVectorK pdf(b);
#endif
      RealVectorK averageAcceptanceRatio(0);
      RealVectorKm1 averageExchangeRatio(0);
      UintVectorKm1 averageExchangeCount(0);
        for
          each(ThreadData const &data in threadData) {
            pdf += data.b;
            averageAcceptanceRatio += data.acceptanceRatio;
            averageExchangeRatio += data.exchangeRatio;
            averageExchangeCount += data.exchangeCount;
          }

        for (int k = 0; k < K; ++k) {
          pdf[k] = std::log(pdf[k]);
        }
        pdf -= pdf[K - 1];
#if defined(CONFIG_MULTIPLE_INTEGRAION)
        for (int k = K - 2; k >= 0; --k) {
          pdf[k] += pdf[k + 1];
        }
#endif
        for (int k = 0; k < K; ++k) {
          pdf[k] = std::exp(pdf[k]);
        }
        for (int k = 0; k < K; ++k) {
          averageAcceptanceRatio[k] =
              std::exp(std::log(averageAcceptanceRatio[k]) - fLogMutation);
        }
        for (int k = 0; k < K - 1; ++k) {
          averageExchangeRatio[k] =
              (averageExchangeCount[k] > 0)
                  ? std::exp(std::log((averageExchangeRatio[k])) -
                             std::log(tgir::Real(averageExchangeCount[k])))
                  : tgir::Real(0);
        }

        ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒")
                                  _TEXT(" %g[mpp]") _TEXT(" b={%g, %g, %g, %g}")
                                      _TEXT(" a={%g, %g, %g, %g}")
                                          _TEXT(" e={%g, %g, %g}"),
                      time / (60 * 60), (time / 60) % 60, time % 60,
                      fMutationPerPixel, pdf[0], pdf[1], pdf[2], pdf[3],
                      averageAcceptanceRatio[0], averageAcceptanceRatio[1],
                      averageAcceptanceRatio[2], averageAcceptanceRatio[3],
                      averageExchangeRatio[0], averageExchangeRatio[1],
                      averageExchangeRatio[2]);

        // save
        bool const is_finished = tgir::IsFinished();
        if (is_finished || (time >= saveTime)) {
          saveTime += intervalTime;

          ::_stprintf_s(filename, MAX_PATH,
                        _TEXT("%s") _TEXT("%02dh%02dm%02ds") _TEXT(" %gmpp")
                            _TEXT(" b={%g, %g, %g}") _TEXT(" a={%g, %g, %g}")
                                _TEXT(" e={%g, %g, %g}") _TEXT(".pfm"),
                        tgir::GetOutdirPath().c_str(), time / (60 * 60),
                        (time / 60) % 60, time % 60, fMutationPerPixel, pdf[0],
                        pdf[1], pdf[2], averageAcceptanceRatio[0],
                        averageAcceptanceRatio[1], averageAcceptanceRatio[2],
                        averageExchangeRatio[0], averageExchangeRatio[1],
                        averageExchangeRatio[2]);

#if defined(CONFIG_APPROX_INITIAL_B)
          pdf *= b_hat;
#endif
          for (std::size_t p = 0; p < M; ++p) {
            m[p] = I[0][p] * pdf[0] + I[1][p] * pdf[1] + I[2][p] * pdf[2] +
                   I[3][p] /* pdf[3]*/;
          }

          //::_ftprintf_s( stderr, filename );
          m.Save(filename, fMutationPerPixel);

          for (std::size_t k = 0; k < K; ++k) {
            ::_stprintf_s(filename, MAX_PATH,
                          _TEXT("%s") _TEXT("%u/") _TEXT("%02dh%02dm%02ds")
                              _TEXT(" %gmpp") _TEXT(".pfm"),
                          tgir::GetOutdirPath().c_str(), k, time / (60 * 60),
                          (time / 60) % 60, time % 60, fMutationPerPixel);
            I[k].Save(filename, fMutationPerPixel / pdf[k]);
          }

          if (is_finished || (time >= renderingTime)) {
            tgir::FinishRendering();
#if defined(SHOW_RESULT)
            ::ShellExecute(hi::nullptr, _TEXT("open"), filename, hi::nullptr,
                           hi::nullptr, SW_SHOWDEFAULT);
#endif
            return;
          }

          beginTime = std::clock() - currentTime;  // 時計を戻す
        }
    }
  }
}

}  // end of namespace tgir
