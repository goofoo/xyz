#include "Evaluator.hpp"
#include "SamplingState.hpp"
#include "core/shell.hpp"
#include "core/Scene.hpp"
#include "core/Camera.hpp"
#include "core/Film.hpp"

namespace tgir {
void FullMLT::run() {
  std::tcerr << _TEXT("** Full Metropolis Light Transport **") << std::endl;

  static std::size_t const N =
      TGIR_CONFIG_kSampleCount /
      TGIR_CONFIG_kSystemCount;  // １サイクルあたりのサンプル数
  static tgir::Real const vMutationWidth = hi::rcp<tgir::Real>(64);  // 変異幅

  tgir::Scene const &scene = tgir::Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  tgir::Film m(W, H);

  std::size_t const M = m.GetCount();

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  struct Replica {
    typedef tgir::SamplingState<tgir::Real> SamplingState;

    std::size_t current;
    SamplingState state[2];

    Replica() : current(0) {}
  } copies[TGIR_CONFIG_kSystemCount];

  struct ThreadData {
    ThreadData() : b(0), n(0), acceptanceRatio(0) {}

    ThreadData(ThreadData const &data)
        : b(data.b), n(data.n), acceptanceRatio(data.acceptanceRatio) {}

    std::mt19937_64 random;
    tgir::Real b;
    tgir::Real n;
    tgir::Real acceptanceRatio;
  };

  std::vector<ThreadData> threadData(omp_get_max_threads());

  // 積分値の逆数
  tgir::Real b = 0;

  // 計測開始
  int const renderingTime = tgir::GetRenderingTime();
  int const intervalTime = tgir::GetIntervalTime();
  int saveTime = intervalTime;
  std::clock_t beginTime = std::clock();

  //
  // Stage1: Initialization.
  //
  std::tcerr << _TEXT("Stage1: Initialization.") << std::endl;
  {
    std::vector<Replica::SamplingState> seeds(CONFIG_SEEDS);

    // cumulative distribution function
    std::vector<tgir::Real> cdf(CONFIG_SEEDS + 1);

    {
// generate path seeds
#pragma omp parallel for
      for (int i = 0; i < CONFIG_SEEDS; ++i) {
        ThreadData &data = threadData[omp_get_thread_num()];
        seeds[i].sample.Init(data.random);
        scene.Evalute(seeds[i].sample, &seeds[i].paths);
        seeds[i].paths.Values(&seeds[i].values);
        seeds[i].q = seeds[i].paths.PDFFullMLT();
        cdf[i + 1] = seeds[i].q;
      }

      std::size_t n = 0;  // 生成できたサンプル数
      cdf[0] = 0;
      for (int i = 1; i <= CONFIG_SEEDS; ++i) {
        if (cdf[i] > 0) {
          ++n;
        }
        cdf[i] += cdf[i - 1];
      }

      // 正規化定数の表示
      b = cdf.back() / CONFIG_SEEDS;
      ::_ftprintf_s(stderr, _TEXT("  b={%lf}  n={%6u}\n"), b, n);

      if (b <= 0) {
        tgir::FinishRendering();
        return;  // レンダリングする必要なし
      }

      b = hi::rcp(b);  // 逆数を計算しておく
    }

    // choose SamplingState
    for (int c = 0; c < TGIR_CONFIG_kSystemCount; ++c) {
      ThreadData &data = threadData[omp_get_thread_num()];
      tgir::Real const xi =
          cdf.back() *
          (c + data.std::uniform_real_distribution<tgir::Real>()(random)) /
          TGIR_CONFIG_kSystemCount;
      std::vector<tgir::Real>::const_iterator it =
          std::upper_bound(cdf.begin(), cdf.end(), xi);
      assert(it != cdf.begin());
      assert(it != cdf.end());

      std::size_t const i = it - cdf.begin() - 1;
      Replica &replica = copies[c];
      replica.state[replica.current] = seeds[i];
    }
  }

  //
  // Stage2: Computation.
  //
  std::tcerr << _TEXT("Stage2: Computation.") << std::endl;
  {
    std::clock_t const currentTime = std::clock() - beginTime;
    std::clock_t const time =
        static_cast<std::clock_t>(currentTime / double(CLOCKS_PER_SEC) + 0.5);
    ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒\n"), time / (60 * 60),
                  (time / 60) % 60, time % 60);
    beginTime = std::clock();
  }

  for (std::size_t i = 1;; ++i) {
#pragma omp parallel for
    for (int c = 0; c < TGIR_CONFIG_kSystemCount; ++c) {
      ThreadData &data = threadData[omp_get_thread_num()];
      Replica &replica = copies[c];

      for (std::size_t n = 0; n < N; ++n) {
        Replica::SamplingState &oldSample = replica.state[replica.current];
        Replica::SamplingState &newSample = replica.state[replica.current ^ 1];

        // mutation a primary sample point
        EMutationType type = MUTATION_FULL;
        int c = 0;
        switch (static_cast<unsigned>(
            data.std::uniform_real_distribution<tgir::Real>()(random) * 8)) {
        case 0:
        case 1:
          c = 1;
          newSample.sample.Init(data.random);
          break;
        case 2:
        case 3:
          newSample.sample.Mutate(oldSample.sample, vMutationWidth,
                                  data.random);
          break;
        case 4:
          newSample.sample.MutateFocusOnLightSubPath(
              oldSample.sample, vMutationWidth, data.random);
          break;
        case 5:
          newSample.sample.MutateFocusOnEyeSubPath(oldSample.sample,
                                                   vMutationWidth, data.random);
          break;
        case 6:
          newSample.sample.MutateLightSubPath(oldSample.sample, vMutationWidth,
                                              data.random);
          newSample.paths = oldSample.paths;
          type = MUTATION_LIGHT;
          break;
        case 7:
          newSample.sample.MutateEyeSubPath(oldSample.sample, vMutationWidth,
                                            data.random);
          newSample.paths = oldSample.paths;
          type = MUTATION_EYE;
          break;
        default:
          std::abort();
        }

        // evalute the paths
        scene.Evalute(newSample.sample, &newSample.paths, type);
        newSample.paths.Values(&newSample.values);
        newSample.q = newSample.paths.PDFFullMLT();
        if (c) {
          data.b += newSample.q;
          data.n++;
        }

        // acceptance probability
        tgir::Real const a = tgir::SafeDivide(oldSample.q, newSample.q);
        data.acceptanceRatio += a;

        // deposite
        if (oldSample.q > 0) {
          m.Deposite(oldSample.values, (1 - a) * hi::rcp(oldSample.q * b));
        }
        if (newSample.q > 0) {
          m.Deposite(newSample.values, (a)*hi::rcp(newSample.q * b));
        }

        // accept/reject
        if (data.std::uniform_real_distribution<tgir::Real>()(random) < a) {
          replica.current ^= 1;  // accept
        }
      }
    }

    // 定期的な保存
    {
      std::clock_t const currentTime = std::clock() - beginTime;
      std::clock_t const time =
          static_cast<std::clock_t>(currentTime / double(CLOCKS_PER_SEC) + 0.5);
      hi::ulong const mutation = hi::ulong(i) * TGIR_CONFIG_kSystemCount * N;
      tgir::Real const fLogMutation = std::log(tgir::Real(mutation));
      tgir::Real const fMutationPerPixel =
          std::exp(fLogMutation - std::log(tgir::Real(M)));

      tgir::Real actual_b = 0;
      tgir::Real actual_n = 0;
      tgir::Real averageAcceptanceRatio = 0;
        for
          each(ThreadData const &data in threadData) {
            actual_b += data.b;
            actual_n += data.n;
            averageAcceptanceRatio += data.acceptanceRatio;
          }
        averageAcceptanceRatio =
            std::exp(std::log(averageAcceptanceRatio) - fLogMutation);

        ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒")
                                  _TEXT(" % 8.2lf[mpp]") _TEXT(" a={%5.4lf}"),
                      time / (60 * 60), (time / 60) % 60, time % 60,
                      fMutationPerPixel, averageAcceptanceRatio);

        // save
        bool const is_finished = tgir::IsFinished();
        if (is_finished || (time >= saveTime)) {
          saveTime += intervalTime;

          ::_stprintf_s(
              filename, MAX_PATH,
              _TEXT("%s") _TEXT("%02dh%02dm%02ds") _TEXT(" %08.2lfmpp")
                  _TEXT(" b={%6.4lf}") _TEXT(" a={%5.4lf}") _TEXT(".pfm"),
              tgir::GetOutdirPath().c_str(), time / (60 * 60), (time / 60) % 60,
              time % 60, fMutationPerPixel, hi::rcp(b), averageAcceptanceRatio);

          tgir::Real const scale = b * actual_b / actual_n;
          m.Save(filename, fMutationPerPixel / scale);

          if (is_finished || (time >= renderingTime)) {
            tgir::FinishRendering();
#ifdef SHOW_RESULT
            ::ShellExecute(nullptr, _TEXT("open"), filename, nullptr, nullptr,
                           SW_SHOWDEFAULT);
#endif SHOW_RESULT
            return;
          }

          beginTime = std::clock() - currentTime;  // 時計を戻す
        }
    }
  }
}

}  // end of namespace tgir
