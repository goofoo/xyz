#include "Evaluator.hpp"
#include "SamplingState.hpp"
#include "core/shell.hpp"
#include "core/Scene.hpp"
#include "core/Camera.hpp"
#include "core/Film.hpp"

namespace tgir {
void BPT::run() {
  std::tcerr << _TEXT("** Bidirectional Path Tracing **") << std::endl;

  tgir::Scene const &scene = tgir::Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  tgir::Film film(W, H);

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  struct ThreadData {
    ThreadData()
        : sample(TGIR_CONFIG_kMaxPathLength),
          path(TGIR_CONFIG_kMaxPathLength),
          values(TGIR_CONFIG_kMaxPathLength) {}

    std::mt19937_64 random;
    tgir::PrimarySample sample;  // サンプル
    tgir::Path path;             // サンプルで生成される経路
    std::vector<tgir::PixelDescriptor> values;  // サンプルの評価値
  };

  std::vector<ThreadData> threadData(omp_get_max_threads());

  // 計測開始
  int const renderingTime = tgir::GetRenderingTime();
  int const intervalTime = tgir::GetIntervalTime();

  std::tcerr << _TEXT("  rendering time: ") << renderingTime << std::endl;
  std::tcerr << _TEXT("  interval  time: ") << intervalTime << std::endl;

  int saveTime = intervalTime;
  std::clock_t beginTime = std::clock();

  //
  // Stage1: Computation.
  //
  std::tcerr << _TEXT("Stage1: Computation.") << std::endl;
  hi::basic_vector2<int> const pattern[4] = {
      hi::basic_vector2<int>(0, 0), hi::basic_vector2<int>(1, 1),
      hi::basic_vector2<int>(1, 0), hi::basic_vector2<int>(0, 1),
  };

  for (std::size_t i = 1;; ++i) {
    for (int n = 0; n < 4; ++n) {
#pragma omp parallel for
      for (int y = 0; y < H; y += 2) {
        ThreadData &data = threadData[omp_get_thread_num()];
        for (int x = 0; x < W; x += 2) {
          // mutation
          data.sample.Init(data.random);
          data.sample.Stratify(x + pattern[n][0], W, y + pattern[n][1], H);

          // evalution
          scene.Evalute(data.sample, &data.path);
          data.path.Values(&data.values);

          // deposite
          film.Deposite(data.values);
        }
      }

      // 定期的な保存
      {
        std::clock_t const currentTime = std::clock() - beginTime;
        std::clock_t const time = static_cast<std::clock_t>(
            currentTime / double(CLOCKS_PER_SEC) + 0.5);
        tgir::Real const spp = i + (n + 1) * tgir::Real(0.25);

        ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 % 8.2lf[spp]"),
                      time / (60 * 60), (time / 60) % 60, time % 60, spp);

        // save
        bool const is_finished = tgir::IsFinished();
        if (is_finished || (time >= saveTime)) {
          saveTime += intervalTime;

          ::_stprintf_s(filename, MAX_PATH,
                        _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                            _TEXT(" %08.2lfspp") _TEXT(".pfm"),
                        tgir::GetOutdirPath().c_str(), time / (60 * 60),
                        (time / 60) % 60, time % 60, spp);

          film.Save(filename, spp);

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
}
}  // end of namespace tgir
