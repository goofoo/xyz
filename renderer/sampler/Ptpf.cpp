#include "Strategy.hpp"
#include "SamplingState.hpp"
#include "core/shell.hpp"
#include "core/Scene.hpp"
#include "core/Camera.hpp"
#include "core/Film.hpp"

namespace tgir {
void Ptpf1::run() {
  std::tcerr << _TEXT("** Path Tracing with Particle-Filter **") << std::endl;

  tgir::Scene const &scene = tgir::Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  if ((W % CONFIG_BLOCK_SIZE) || (H % CONFIG_BLOCK_SIZE)) {
    std::tcerr << _TEXT("不正な解像度です．") << std::endl
               << _TEXT("幅と高さをそれぞれ") << CONFIG_BLOCK_SIZE
               << _TEXT("の倍数にしてください．") << std::endl;
    return;
  }

  tgir::Film film(W, H);

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  struct ThreadData {
    ThreadData()
        : aColor(CONFIG_FILTER_SIZE),
          info(CONFIG_FILTER_SIZE),
          ray(CONFIG_FILTER_SIZE),
          tmp(CONFIG_FILTER_SIZE),
          cdf(CONFIG_FILTER_SIZE + 1) {}

    std::mt19937_64 random;
    std::vector<tgir::SpectrumVector> aColor;
    std::vector<tgir::ParticleInfo> info;
    std::vector<tgir::Particle> ray;
    std::vector<tgir::Particle> tmp;
    std::vector<tgir::Real> cdf;  // cumulative distribution function
  };

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  std::vector<ThreadData> threadData(omp_get_max_threads());

  // 計測開始
  int const renderingTime = tgir::GetRenderingTime();
  int const intervalTime = tgir::GetIntervalTime();
  int saveTime = intervalTime;
  std::clock_t beginTime = std::clock();

  //
  // Stage1: Computation.
  //
  std::tcerr << _TEXT("Stage1: Computation.") << std::endl;
  scene.n_FindIntersection = 0;
  for (std::size_t spp = 1;; ++spp) {
#pragma omp parallel for
    for (int y = 0; y < H; y += CONFIG_BLOCK_SIZE) {
      ThreadData &data = threadData[omp_get_thread_num()];
      for (int x = 0; x < W; x += CONFIG_BLOCK_SIZE) {
        scene.EvaluteParticleFilter(x, y, data.aColor, data.info, data.ray,
                                    data.tmp, data.cdf, data.random);
        for (std::size_t i = 0, dy = 0; dy < CONFIG_BLOCK_SIZE; ++dy)
          for (std::size_t dx = 0; dx < CONFIG_BLOCK_SIZE; ++dx, ++i) {
            film[(y + dy) * W + (x + dx)] += data.aColor[i];
          }
      }
    }

    // 定期的な保存
    {
      std::clock_t const currentTime = std::clock() - beginTime;
      std::clock_t const time =
          static_cast<std::clock_t>(currentTime / double(CLOCKS_PER_SEC) + 0.5);

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %+6u[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60, spp);

      // save
      bool const is_finished = tgir::IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += intervalTime;

        ::_stprintf_s(
            filename, MAX_PATH, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                    _TEXT(" %06uspp %012lu") _TEXT(".pfm"),
            tgir::GetOutdirPath().c_str(), time / (60 * 60), (time / 60) % 60,
            time % 60, spp, scene.n_FindIntersection);

        film.Save(filename, spp);

        if (is_finished || (time >= renderingTime)) {
          tgir::FinishRendering();
          break;
        }

        beginTime = std::clock() - currentTime;  // 時計を戻す
      }
    }
  }
#ifdef SHOW_RESULT
  ::ShellExecute(nullptr, _TEXT("open"), filename, nullptr, nullptr,
                 SW_SHOWDEFAULT);
#endif SHOW_RESULT
}
}  // end of namespace tgir

namespace tgir {
void Ptpf2::run() {
  std::tcerr << _TEXT("** Path Tracing with Particle-Filter **") << std::endl;

  tgir::Scene const &scene = tgir::Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  if ((W % CONFIG_BLOCK_SIZE) || (H % CONFIG_BLOCK_SIZE)) {
    std::tcerr << _TEXT("不正な解像度です．") << std::endl
               << _TEXT("幅と高さをそれぞれ") << CONFIG_BLOCK_SIZE
               << _TEXT("の倍数にしてください．") << std::endl;
    return;
  }

  tgir::Film film(W, H);

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  struct ThreadData {
    ThreadData()
        : info(CONFIG_FILTER_SIZE),
          ray(CONFIG_FILTER_SIZE),
          tmp(CONFIG_FILTER_SIZE),
          cdf(CONFIG_FILTER_SIZE + 1) {}

    std::mt19937_64 random;
    std::vector<tgir::ParticleInfo> info;
    std::vector<tgir::Particle> ray;
    std::vector<tgir::Particle> tmp;
    std::vector<tgir::Real> cdf;  // cumulative distribution function
  };

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  std::vector<ThreadData> threadData(omp_get_max_threads());

  // 計測開始
  int const renderingTime = tgir::GetRenderingTime();
  int const intervalTime = tgir::GetIntervalTime();
  int saveTime = intervalTime;
  std::clock_t beginTime = std::clock();

  //
  // Stage1: Computation.
  //
  tgir::SpectrumVector sColor;

  std::tcerr << _TEXT("Stage1: Computation.") << std::endl;
  scene.n_FindIntersection = 0;
  for (std::size_t spp = 1;; spp += CONFIG_FILTER_SIZE) {
#pragma omp parallel for private(sColor)
    for (int y = 0; y < H; ++y) {
      ThreadData &data = threadData[omp_get_thread_num()];
      for (int x = 0; x < W; ++x) {
        scene.EvaluteParticleFilter(x, y, sColor, data.info, data.ray, data.tmp,
                                    data.cdf, data.random);
        film[y * W + x] += sColor;
      }
    }

    // 定期的な保存
    {
      std::clock_t const currentTime = std::clock() - beginTime;
      std::clock_t const time =
          static_cast<std::clock_t>(currentTime / double(CLOCKS_PER_SEC) + 0.5);

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %+6u[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60, spp);

      // save
      bool const is_finished = tgir::IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += intervalTime;

        ::_stprintf_s(
            filename, MAX_PATH, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                    _TEXT(" %06uspp %012lu") _TEXT(".pfm"),
            tgir::GetOutdirPath().c_str(), time / (60 * 60), (time / 60) % 60,
            time % 60, spp, scene.n_FindIntersection);

        film.Save(filename, spp);

        if (is_finished || (time >= renderingTime)) {
          tgir::FinishRendering();
          break;
        }

        beginTime = std::clock() - currentTime;  // 時計を戻す
      }
    }
  }
#ifdef SHOW_RESULT
  ::ShellExecute(nullptr, _TEXT("open"), filename, nullptr, nullptr,
                 SW_SHOWDEFAULT);
#endif SHOW_RESULT
}
}  // end of namespace tgir
