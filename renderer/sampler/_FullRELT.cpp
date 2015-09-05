#include "Evaluator.hpp"
#include "SamplingState.hpp"
#include "core/shell.hpp"
#include "core/Scene.hpp"
#include "core/Camera.hpp"
#include "core/Film.hpp"

namespace {
#if defined(TGIR_CONFIG_TRIPLE_REPLICS)
typedef hi::basic_vector3<tgir::Real> RealVectorK;
typedef hi::basic_vector2<tgir::Real> RealVectorKm1;
typedef hi::basic_vector3<hi::uint> UintVectorK;
typedef hi::basic_vector2<hi::uint> UintVectorKm1;
#else
typedef hi::basic_vector4<tgir::Real> RealVectorK;
typedef hi::basic_vector3<tgir::Real> RealVectorKm1;
typedef hi::basic_vector4<hi::uint> UintVectorK;
typedef hi::basic_vector3<hi::uint> UintVectorKm1;
#endif

std::size_t const K = RealVectorK::N;
std::size_t const N = TGIR_CONFIG_kSampleCount / (TGIR_CONFIG_kSystemCount * K);
std::size_t const SEEDS_COUNT = CONFIG_SEEDS;

#if defined(TGIR_CONFIG_TRIPLE_REPLICS)
RealVectorKm1 const vMutationWidth(hi::rcp<tgir::Real>(96),
                                   hi::rcp<tgir::Real>(64));  // 変異幅
#else
RealVectorKm1 const vMutationWidth(hi::rcp<tgir::Real>(64),
                                   hi::rcp<tgir::Real>(96),
                                   hi::rcp<tgir::Real>(64));  // 変異幅
#endif

struct Replica {
  typedef tgir::SamplingState<RealVectorK> SamplingState;

  std::size_t current;
  SamplingState state[2];

  Replica() : current(0) {}
};

struct System {
  Replica replicas[K];
  std::size_t indices[K];

  System() {
    for (std::size_t k = 0; k < K; ++k) {
      indices[k] = k;
    }
  }
};

struct ThreadData {
  hi::mt19937 random;
  tgir::Path path;
  RealVectorK b;
  RealVectorK acceptanceRatio;
  RealVectorKm1 exchangeRatio;
  UintVectorKm1 exchangeCount;

  ThreadData()
      : path(TGIR_CONFIG_kMaxPathLength),
        b(0),
        acceptanceRatio(0),
        exchangeRatio(0),
        exchangeCount(0) {}
};

void CalculateCumulativeDistributionFunction(
    std::size_t const k, std::vector<tgir::Real> *const pCDF,
    UintVectorK *const n) {
  for (std::size_t i = 1; i <= SEEDS_COUNT; ++i) {
    if ((*pCDF)[i] > 0) {
      ++(*n)[k];
    }
    (*pCDF)[i] += (*pCDF)[i - 1];
  }
}

void DoPopulationResampling(std::vector<tgir::Real> const &rCDF,
                            std::vector<ThreadData> *const pThreadData,
                            std::vector<Replica::SamplingState> *const aSeeds,
                            std::size_t *const pCurrent) {
  std::size_t j = 0;
  for (std::size_t i = 0; i < SEEDS_COUNT; ++i) {
    ThreadData &data = (*pThreadData)[omp_get_thread_num()];
    tgir::Real const xi =
        rCDF.back() * (i + data.random.next<tgir::Real>()) / SEEDS_COUNT;
    while ((rCDF[j] > xi) || (xi >= rCDF[j + 1])) ++j;
    aSeeds[(*pCurrent) ^ 1][i] = aSeeds[(*pCurrent)][j];
  }
  (*pCurrent) ^= 1;
}

void ChooseSamples(std::size_t const k,
                   std::vector<Replica::SamplingState> const &seed,
                   System *const aCopies) {
  for (std::size_t c = 0; c < TGIR_CONFIG_kSystemCount; ++c) {
    Replica *const replica = &aCopies[c].replicas[aCopies[c].indices[k]];
    replica->state[replica->current] =
        seed[seed.size() * c / TGIR_CONFIG_kSystemCount];
  }
}

void DoMetropolisSampling(tgir::Scene const &rScene,
                          std::vector<ThreadData> *const pThreadData,
                          std::size_t const k,
                          std::vector<tgir::Real> *const pCDF,
                          std::vector<Replica::SamplingState> *const aSeeds,
                          std::size_t *const pCurrent) {
#pragma omp parallel for
  for (int i = 0; i < SEEDS_COUNT; ++i) {
    ThreadData &data = (*pThreadData)[omp_get_thread_num()];

    Replica::SamplingState &oldSample = aSeeds[(*pCurrent)][i];
    Replica::SamplingState &newSample = aSeeds[(*pCurrent) ^ 1][i];

// mutate a path
#if defined(TGIR_CONFIG_TRIPLE_REPLICS)
    if (k < 1)
#else
    if (k < 2)
#endif
    {
      newSample.sample.MutateFocusOnEyeSubPath(oldSample.sample,
                                               vMutationWidth[k], data.random);
    } else {
      newSample.sample.MutateFocusOnLightSubPath(
          oldSample.sample, vMutationWidth[k], data.random);
    }

    // evalute the path
    rScene.Evalute(newSample.sample, data.path, newSample.values);
    data.path.PDFsFullRELT(&newSample.q);

    // compute acceptance probability
    tgir::Real const a = tgir::SafeDivide(oldSample.q[k], newSample.q[k]);

    data.b[k - 1] += (oldSample.q[k] > 0)
                         ? (1 - a) * oldSample.q[k - 1] / oldSample.q[k]
                         : 0;
    data.b[k - 1] +=
        (newSample.q[k] > 0) ? (a)*newSample.q[k - 1] / newSample.q[k] : 0;

    if (data.random.next<tgir::Real>() >= a) {
      // rejected new sample
      newSample.swap(oldSample);
    }
    (*pCDF)[i + 1] =
        (newSample.q[k] > 0) ? newSample.q[k - 1] / newSample.q[k] : 0;
  }
  (*pCurrent) ^= 1;
}
}
// end of unnamed namespace

namespace tgir {
void FullRELT::run() {
  std::tcerr << _TEXT("** Full Replica-Exchange Light Transport **")
             << std::endl;

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

  System copies[TGIR_CONFIG_kSystemCount];

  std::vector<ThreadData> threadData(omp_get_max_threads());

  // 積分値の逆数
  RealVectorK b_hat(0);

  // 計測開始
  int const renderingTime = tgir::GetRenderingTime();
  int const intervalTime = tgir::GetIntervalTime();
  int saveTime = intervalTime;
  std::clock_t beginTime = std::clock();

  //
  // Stage1: Initialization.
  //
  std::tcerr << _TEXT("Stage1: Initialization.") << std::endl;

  //
  // Population Annealing
  //
  {
    std::size_t current = 0;
    std::vector<Replica::SamplingState> aSeeds[2] = {
        std::vector<Replica::SamplingState>(SEEDS_COUNT),
        std::vector<Replica::SamplingState>(SEEDS_COUNT),
    };

    // cumulative distribution function
    std::vector<tgir::Real> cdf(SEEDS_COUNT + 1);
    cdf[0] = 0;

    UintVectorK n(0);
    n[K - 1] = SEEDS_COUNT;

    // 初期サンプルの生成
    {
      int k = K - 1;
#pragma omp parallel for
      for (int i = 0; i < SEEDS_COUNT; ++i) {
        ThreadData &data = threadData[omp_get_thread_num()];
        Replica::SamplingState &newSample = aSeeds[current][i];
        newSample.sample.Init(data.random);
        scene.Evalute(newSample.sample, data.path, newSample.values);
        data.path.PDFsFullRELT(&newSample.q);
        data.b[k] += newSample.q[k];
        data.b[k - 1] += newSample.q[k - 1] / newSample.q[k];
        cdf[i + 1] = newSample.q[k - 1] / newSample.q[k];
      }
      CalculateCumulativeDistributionFunction(k - 1, &cdf, &n);
      if (cdf.back() <= 0) {
        std::tcerr << _TEXT("初期サンプルの生成に失敗しました．") << std::endl;
        tgir::FinishRendering();
        return;  // レンダリングする必要なし
      }
      DoPopulationResampling(cdf, &threadData, aSeeds, &current);
      ChooseSamples(k - 1, aSeeds[current], copies);
    }

    for (int k = K - 2; k > 0; --k) {
      DoMetropolisSampling(scene, &threadData, k, &cdf, aSeeds, &current);
      CalculateCumulativeDistributionFunction(k - 1, &cdf, &n);
      if (cdf.back() <= 0) {
        std::tcerr << _TEXT("初期サンプルの生成に失敗しました．") << std::endl;
        tgir::FinishRendering();
        return;  // レンダリングする必要なし
      }
      DoPopulationResampling(cdf, &threadData, aSeeds, &current);
      ChooseSamples(k - 1, aSeeds[current], copies);
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
#if defined(TGIR_CONFIG_TRIPLE_REPLICS)
    ::_ftprintf_s(stdout, _TEXT("  b_hat={%g, %g, %g} n={%6u, %6u, %6u}\n"),
                  b_hat[0], b_hat[1], b_hat[2], n[0], n[1], n[2]);
    if ((b_hat[0] <= 0) || (b_hat[1] <= 0) || (b_hat[2] <= 0)) {
      tgir::FinishRendering();
      return;  // レンダリングする必要なし
    }
#else
    ::_ftprintf_s(
        stdout, _TEXT("  b_hat={%g, %g, %g, %g} n={%6u, %6u, %6u, %6u}\n"),
        b_hat[0], b_hat[1], b_hat[2], b_hat[3], n[0], n[1], n[2], n[3]);
    if ((b_hat[0] <= 0) || (b_hat[1] <= 0) || (b_hat[2] <= 0) ||
        (b_hat[3] <= 0)) {
      tgir::FinishRendering();
      return;  // レンダリングする必要なし
    }
#endif

    b_hat = hi::rcp(b_hat);  // 逆数を求める
  }

  {
    std::clock_t const currentTime = std::clock() - beginTime;
    std::clock_t const time =
        static_cast<std::clock_t>(currentTime / double(CLOCKS_PER_SEC) + 0.5);
    ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒\n"), time / (60 * 60),
                  (time / 60) % 60, time % 60);
    beginTime = std::clock();
  }

#if 0
    RealVectorK const c_hat = b_hat * RealVectorK( 0.25, 0.5, 0.5, 1 );
#else
  RealVectorK const c_hat = b_hat;
#endif

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
#if defined(TGIR_CONFIG_TRIPLE_REPLICS)
          switch (k) {
          case 0:
            switch (n & 1) {
            case 0:
              newSample.sample.MutateFocusOnEyeSubPath(
                  oldSample.sample, vMutationWidth[k], data.random);
              break;
            case 1:
              newSample.sample.Mutate(oldSample.sample, vMutationWidth[k],
                                      data.random);
              break;
            }
          case 2:
            switch (n & 3) {
            case 0:
              newSample.sample.MutateFocusOnEyeSubPath(
                  oldSample.sample, vMutationWidth[k], data.random);
              break;
            case 1:
              newSample.sample.Mutate(oldSample.sample, vMutationWidth[k],
                                      data.random);
              break;
            case 2:
              newSample.sample.MutateFocusOnLightSubPath(
                  oldSample.sample, vMutationWidth[k], data.random);
              break;
            }
            break;
          }

          scene.Evalute(newSample.sample, data.path,
                        newSample.values);  // evalute the path
#else

          EMutationType type = MUTATION_FULL;

          switch (k) {
          case 0:
          case 1:
            if (n & 1)
              newSample.sample.MutateFocusOnEyeSubPath(
                  oldSample.sample, vMutationWidth[k], data.random);
            else
              newSample.sample.MutateFocusOnLightSubPath(
                  oldSample.sample, vMutationWidth[k], data.random);
            break;
          case 2:
            switch (n % 6) {
            case 0:
            case 1:
              newSample.sample.Mutate(oldSample.sample, vMutationWidth[k],
                                      data.random);
              break;
            case 2:
              newSample.sample.MutateFocusOnEyeSubPath(
                  oldSample.sample, vMutationWidth[k], data.random);
              break;
            case 3:
              newSample.sample.MutateFocusOnLightSubPath(
                  oldSample.sample, vMutationWidth[k], data.random);
              break;
            case 4:
              newSample.sample.MutateLightSubPath(
                  oldSample.sample, vMutationWidth[k], data.random);
              type = MUTATION_LIGHT;
              break;
            case 5:
              newSample.sample.MutateEyeSubPath(oldSample.sample,
                                                vMutationWidth[k], data.random);
              type = MUTATION_EYE;
              break;
            }
            break;
          }

          scene.Evalute(newSample.sample, data.path,
                        newSample.values);  // evalute the path
#endif

          data.path.PDFsFullRELT(&newSample.q);

          tgir::Real const a = tgir::SafeDivide(
              oldSample.q[k], newSample.q[k]);  // acceptance probability
          data.acceptanceRatio[k] += a;
#if defined(CONFIG_APPROX_INITIAL_B)
          I[k].Deposite(
              oldSample.values,
              (1 - a) *
                  hi::rcp(hi::dot(oldSample.q, c_hat)));  // deposite old sample
          I[k].Deposite(
              newSample.values,
              (a)*hi::rcp(hi::dot(newSample.q, c_hat)));  // deposite new sample
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
          data.path.PDFsFullRELT(&newSample.q);

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
                            newSample.q, c_hat)));  // deposite new sample
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
          tgir::Real const e = tgir::SafeDivide(state0.q[k] * state1.q[k + 1],
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

#if defined(TGIR_CONFIG_TRIPLE_REPLICS)
        ::_ftprintf_s(stderr,
                      _TEXT("\r  %02u時間%02u分%02u秒") _TEXT(" %8.4lf[mpp]")
                          _TEXT(" b={%8.7lf, %8.7lf}")
                              _TEXT(" a={%8.7lf, %8.7lf}")
                                  _TEXT(" e={%8.7lf, %8.7lf}"),
                      time / (60 * 60), (time / 60) % 60, time % 60,
                      fMutationPerPixel, pdf[0], pdf[1],
                      averageAcceptanceRatio[0], averageAcceptanceRatio[1],
                      averageExchangeRatio[0], averageExchangeRatio[1]);
#else
        ::_ftprintf_s(stderr,
                      _TEXT("\r  %02u時間%02u分%02u秒") _TEXT(" %8.4lf[mpp]")
                          _TEXT(" b={%8.7lf, %8.7lf, %8.7lf}")
                              _TEXT(" a={%8.7lf, %8.7lf, %8.7lf}")
                                  _TEXT(" e={%8.7lf, %8.7lf, %8.7lf}"),
                      time / (60 * 60), (time / 60) % 60, time % 60,
                      fMutationPerPixel, pdf[0], pdf[1], pdf[2],
                      averageAcceptanceRatio[0], averageAcceptanceRatio[1],
                      averageAcceptanceRatio[2], averageExchangeRatio[0],
                      averageExchangeRatio[1], averageExchangeRatio[2]);
#endif

        // save
        bool const is_finished = tgir::IsFinished();
        if (is_finished || (time >= saveTime)) {
          saveTime += intervalTime;

#if defined(TGIR_CONFIG_TRIPLE_REPLICS)
          ::_stprintf_s(filename, MAX_PATH,
                        _TEXT("%s") _TEXT("%02dh%02dm%02ds") _TEXT(" %gmpp")
                            _TEXT(" b={%g, %g}") _TEXT(" a={%g, %g}")
                                _TEXT(" e={%g, %g}") _TEXT(".pfm"),
                        tgir::GetOutdirPath().c_str(), time / (60 * 60),
                        (time / 60) % 60, time % 60, fMutationPerPixel, pdf[0],
                        pdf[1], averageAcceptanceRatio[0],
                        averageAcceptanceRatio[1], averageExchangeRatio[0],
                        averageExchangeRatio[1]);
#else
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
#endif

#if defined(CONFIG_APPROX_INITIAL_B)
          pdf *= b_hat;
#endif
          for (std::size_t p = 0; p < M; ++p) {
#if defined(TGIR_CONFIG_TRIPLE_REPLICS)
            m[p] = I[0][p] * pdf[0] + I[1][p] * pdf[1] + I[2][p] /* pdf[2] */;
#else
            m[p] = I[0][p] * pdf[0] + I[1][p] * pdf[1] + I[2][p] * pdf[2] +
                   I[3][p] /* pdf[3]*/;
#endif
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
}
// end of namespace tgir
