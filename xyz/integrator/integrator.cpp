#include "integrator.hpp"
#include "../core/shell.hpp"
#include "../core/scene.hpp"
#include "../core/thinlenscamera.hpp"
#include "../core/imagefilm.hpp"

#include "../bsdf/bsdf.hpp"

using namespace xyz;

namespace {
std::mt19937_64 &get_instance_of_random_number_generator_for_single_thread() {
  static std::mt19937_64 instance_of_random_number_generator;
  return instance_of_random_number_generator;
}

xyz::float_t generate_random_number_for_single_thread() {
  return get_instance_of_random_number_generator_for_single_thread()
      .next<xyz::float_t>();
}

class RSample : public IPrimarySample {
 public:
  virtual float_t next() { return generate_random_number_for_single_thread(); }
};

class MSample : public IPrimarySample {
 public:
  inline MSample() : uIterator_(0U), bLargeStep_(true), newPoint_(2) {}

  inline MSample &operator=(__in MSample const &rhs) {
    if (this != &rhs) {
      oldPoint_ = rhs.oldPoint_;
      newPoint_ = rhs.newPoint_;
      uIterator_ = rhs.uIterator_;
      bLargeStep_ = rhs.bLargeStep_;
    }
    return *this;
  }

  virtual float_t next() {
    if (uIterator_ >= newPoint_.size()) {
      newPoint_.push_back(float_t());
    }

    if (bLargeStep_) {
      newPoint_[uIterator_] = generate_random_number_for_single_thread();
    } else {
      if (uIterator_ >= oldPoint_.size()) {
        oldPoint_.push_back(generate_random_number_for_single_thread());
      }

#if 1
      static float_t const s1 = 1.0 / 1024;
      static float_t const s2 = 1.0 / 64;
#else
      static float_t const s1 = 1.0 / 1024;
      static float_t const s2 = 1.0 / 32;
#endif
      static float_t const s3 = -std::log(s2 / s1);

      float_t const dv =
          s2 * std::exp(s3 * generate_random_number_for_single_thread());
      newPoint_[uIterator_] = oldPoint_[uIterator_];
      newPoint_[uIterator_] +=
          (generate_random_number_for_single_thread() < 0.5) ? dv : -dv;
      newPoint_[uIterator_] -= std::floor(newPoint_[uIterator_]);
    }

    return newPoint_[uIterator_++];
  }

  inline void init(__in bool const b) {
    uIterator_ = 0U;
    bLargeStep_ = b;
  }

  inline void accept() { oldPoint_.swap(newPoint_); }

 private:
  std::vector<float_t> oldPoint_;
  std::vector<float_t> newPoint_;
  std::size_t uIterator_;
  bool bLargeStep_;
};
}

void PathTracer::run() {
  ::_ftprintf_s(stderr, _TEXT("** Path Tracing **\n"));

  Scene const &scene = Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  ImageFilm film(W, H);

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  // std::vector<std::mt19937_64> randoms(omp_get_max_threads());
  RSample sample;

  // 計測開始
  int const nRenderingTime = GetRenderingTime();
  int const nIntervalTime = GetIntervalTime();
  int saveTime = nIntervalTime;
  std::clock_t begin_time = std::clock();

  //
  // Stage1: Computation.
  //
  CIE_XYZ_Color color;

  ::_ftprintf_s(stderr, _TEXT(">>>> Computation <<<<\n"));
  scene.n_FindIntersection = 0;
  for (std::size_t spp = 1;; ++spp) {
    //#pragma omp parallel for private(color)
    for (int y = 0; y < H; ++y) {
      // std::mt19937_64 & random = randoms[omp_get_thread_num()];
      for (int x = 0; x < W; ++x) {
        scene.EvalutePathTracing((x + sample.next()) / W,
                                 (y + sample.next()) / H, &color, sample);
        film[y * W + x] += color;
      }
    }

    // 定期的な保存
    {
      std::clock_t const current_time = std::clock() - begin_time;
      std::clock_t const time = static_cast<std::clock_t>(
          current_time / double(CLOCKS_PER_SEC) + 0.5);

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %+6u[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60, spp);

      // save
      bool const is_finished = IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += nIntervalTime;

        ::_stprintf_s(
            filename, MAX_PATH, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                    _TEXT(" %06uspp %012lu") _TEXT(".pfm"),
            GetOutdirPath().c_str(), time / (60 * 60), (time / 60) % 60,
            time % 60, spp, scene.n_FindIntersection);

        film.Save(filename, spp);

        if (is_finished || (time >= nRenderingTime)) {
          FinishRendering();
          break;
        }

        begin_time = std::clock() - current_time;  // 時計を戻す
      }
    }
  }
}

void ImportanceDrivenPathTracer::run() {
  ::_ftprintf_s(stderr, _TEXT("** Importance Driven Path Tracing **\n"));

  Scene const &scene = Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  ImageFilm film(W, H);

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  // std::vector<std::mt19937_64> randoms(omp_get_max_threads());
  RSample sample;

  // 計測開始
  int const nRenderingTime = GetRenderingTime();
  int const nIntervalTime = GetIntervalTime();
  int saveTime = nIntervalTime;
  std::clock_t begin_time = std::clock();

  ::_ftprintf_s(stderr, _TEXT(">>>> Build Importance Map <<<<\n"));
  {
    Scene::GetInstance().CreateLightImportanceMap(
        ::get_instance_of_random_number_generator_for_single_thread());
  }

  CIE_XYZ_Color color;

  ::_ftprintf_s(stderr, _TEXT(">>>> Computation <<<<\n"));
  scene.n_FindIntersection = 0;
  for (std::size_t spp = 1;; ++spp) {
    //#pragma omp parallel for private(color)
    for (int y = 0; y < H; ++y) {
      // std::mt19937_64 & random = randoms[omp_get_thread_num()];
      for (int x = 0; x < W; ++x) {
        scene.EvaluteImportanceDrivenPathTracing(
            (x + sample.next()) / W, (y + sample.next()) / H, &color, sample);
        film[y * W + x] += color;
      }
    }

    // 定期的な保存
    {
      std::clock_t const current_time = std::clock() - begin_time;
      std::clock_t const time = static_cast<std::clock_t>(
          current_time / double(CLOCKS_PER_SEC) + 0.5);

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %+6u[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60, spp);

      // save
      bool const is_finished = IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += nIntervalTime;

        ::_stprintf_s(
            filename, MAX_PATH, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                    _TEXT(" %06uspp %012lu") _TEXT(".pfm"),
            GetOutdirPath().c_str(), time / (60 * 60), (time / 60) % 60,
            time % 60, spp, scene.n_FindIntersection);

        film.Save(filename, spp);

        if (is_finished || (time >= nRenderingTime)) {
          FinishRendering();
          break;
        }

        begin_time = std::clock() - current_time;  // 時計を戻す
      }
    }
  }
}

void PathTracer_MetropolisLightTransport::run() {
  ::_ftprintf_s(stderr,
                _TEXT("** Path Tracing & Metropolis Light Transport **\n"));

  Scene const &scene = Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  ImageFilm film(W, H);

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  // std::vector<std::mt19937_64> randoms(omp_get_max_threads());
  pixel_descriptor_t oldColor;
  pixel_descriptor_t newColor;
  MSample sample;

  // 計測開始
  int const nRenderingTime = GetRenderingTime();
  int const nIntervalTime = GetIntervalTime();
  int saveTime = nIntervalTime;
  std::clock_t begin_time = std::clock();

  // 初期化
  ::_ftprintf_s(stderr, _TEXT(">>>> Initialization <<<<\n"));
  float_t image_contribution_luminance = float_t(0);
  {
    static std::size_t const R = 10000;

    std::vector<MSample> samples(R);
    std::vector<pixel_descriptor_t> colors(R);
    std::vector<float_t> cumulative_distribution_function(R + 1);
    cumulative_distribution_function[0] = float_t(0);

    for (std::size_t i = 0; i < R; ++i) {
      samples[i].init(true);

      float_t const x = sample.next();
      float_t const y = sample.next();
      colors[i].first =
          static_cast<std::size_t>(y * H) * W + static_cast<std::size_t>(x * W);

      scene.EvalutePathTracing(x, y, &colors[i].second, sample);

      cumulative_distribution_function[i + 1] =
          cumulative_distribution_function[i] + colors[i].second[1];
    }

    image_contribution_luminance = cumulative_distribution_function.back() / R;
    if (image_contribution_luminance <= float_t(0)) {
      FinishRendering();
      return;
    }

    {
      // float_t const xi = (c + ::generate_random_number_for_single_thread()) /
      // XYZ_CONFIG_kSystemCount;
      float_t const xi = ::generate_random_number_for_single_thread();
      float_t const value = xi * cumulative_distribution_function.back();
      std::vector<float_t>::const_iterator const it =
          std::upper_bound(cumulative_distribution_function.begin(),
                           cumulative_distribution_function.end(), value);

      assert(it != cumulative_distribution_function.begin());
      assert(it != cumulative_distribution_function.end());

      std::size_t const uIndex =
          it - cumulative_distribution_function.begin() - 1;
      sample = samples[uIndex];
      sample.accept();  // new を old に設定
      oldColor = colors[uIndex];
    }
  }
  float_t const c = hi::rcp(image_contribution_luminance);

  static float_t const kLargeStepProb =
      float_t(0.5);  // ラージステップを行う確率

  ::_ftprintf_s(stderr, _TEXT(">>>> Computation <<<<\n"));
  scene.n_FindIntersection = 0;
  for (std::size_t spp = 1;; ++spp) {
    //#pragma omp parallel for private(color)
    for (int n = 0, Q = W * H; n < Q; ++n) {
      bool const b = (::generate_random_number_for_single_thread() <
                      kLargeStepProb);  // 一様分布からサンプリングするか？
      sample.init(b);

      float_t const x = sample.next();
      float_t const y = sample.next();
      newColor.first =
          static_cast<std::size_t>(y * H) * W + static_cast<std::size_t>(x * W);

      scene.EvalutePathTracing(x, y, &newColor.second, sample);

      float_t const a =
          std::min(CIE_XYZ_Color::value_type(1),
                   newColor.second[1] / oldColor.second[1]);  // 受理確率
      film[oldColor.first] +=
          oldColor.second *
          ((1 - a) / (oldColor.second[1] * c + kLargeStepProb));
      film[newColor.first] +=
          newColor.second *
          ((b + a) / (newColor.second[1] * c + kLargeStepProb));

      if (::generate_random_number_for_single_thread() < a) {
        // 受理
        sample.accept();
        oldColor = newColor;
      }
    }

    // 定期的な保存
    {
      std::clock_t const current_time = std::clock() - begin_time;
      std::clock_t const time = static_cast<std::clock_t>(
          current_time / double(CLOCKS_PER_SEC) + 0.5);

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %+6u[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60, spp);

      // save
      bool const is_finished = IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += nIntervalTime;

        ::_stprintf_s(filename, MAX_PATH,
                      _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                          _TEXT(" %06uspp %012lu avg=%g") _TEXT(".pfm"),
                      GetOutdirPath().c_str(), time / (60 * 60),
                      (time / 60) % 60, time % 60, spp,
                      scene.n_FindIntersection, image_contribution_luminance);

        film.Save(filename, spp);

        if (is_finished || (time >= nRenderingTime)) {
          FinishRendering();
          return;
        }

        begin_time = std::clock() - current_time;  // 時計を戻す
      }
    }
  }
}

void ImportanceDrivenPathTracer_MetropolisLightTransport::run() {
  ::_ftprintf_s(stderr, _TEXT(
                            "** Importance Driven Path Tracing & Metropolis "
                            "Light Transport **\n"));

  Scene const &scene = Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  ImageFilm film(W, H);

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  // std::vector<std::mt19937_64> randoms(omp_get_max_threads());
  pixel_descriptor_t oldColor;
  pixel_descriptor_t newColor;
  MSample sample;

  // 計測開始
  int const nRenderingTime = GetRenderingTime();
  int const nIntervalTime = GetIntervalTime();
  int saveTime = nIntervalTime;
  std::clock_t begin_time = std::clock();

  ::_ftprintf_s(stderr, _TEXT(">>>> Build Importance Map <<<<\n"));
  {
    Scene::GetInstance().CreateLightImportanceMap(
        ::get_instance_of_random_number_generator_for_single_thread());
  }

  // 初期化
  ::_ftprintf_s(stderr, _TEXT(">>>> Initialization <<<<\n"));
  float_t image_contribution_luminance = float_t(0);
  {
    static std::size_t const R = 10000;

    std::vector<MSample> samples(R);
    std::vector<pixel_descriptor_t> colors(R);
    std::vector<float_t> cumulative_distribution_function(R + 1);
    cumulative_distribution_function[0] = float_t(0);

    for (std::size_t i = 0; i < R; ++i) {
      samples[i].init(true);

      float_t const x = sample.next();
      float_t const y = sample.next();
      colors[i].first =
          static_cast<std::size_t>(y * H) * W + static_cast<std::size_t>(x * W);

      scene.EvaluteImportanceDrivenPathTracing(x, y, &colors[i].second, sample);

      cumulative_distribution_function[i + 1] =
          cumulative_distribution_function[i] + colors[i].second[1];
    }

    image_contribution_luminance = cumulative_distribution_function.back() / R;
    if (image_contribution_luminance <= float_t(0)) {
      FinishRendering();
      return;
    }

    {
      // float_t const xi = (c + ::generate_random_number_for_single_thread()) /
      // XYZ_CONFIG_kSystemCount;
      float_t const xi = ::generate_random_number_for_single_thread();
      float_t const value = xi * cumulative_distribution_function.back();
      std::vector<float_t>::const_iterator const it =
          std::upper_bound(cumulative_distribution_function.begin(),
                           cumulative_distribution_function.end(), value);

      assert(it != cumulative_distribution_function.begin());
      assert(it != cumulative_distribution_function.end());

      std::size_t const uIndex =
          it - cumulative_distribution_function.begin() - 1;
      sample = samples[uIndex];
      sample.accept();  // new を old に設定
      oldColor = colors[uIndex];
    }
  }
  float_t const c = hi::rcp(image_contribution_luminance);

  static float_t const kLargeStepProb =
      float_t(0.5);  // ラージステップを行う確率

  ::_ftprintf_s(stderr, _TEXT(">>>> Computation <<<<\n"));
  scene.n_FindIntersection = 0;
  for (std::size_t spp = 1;; ++spp) {
    //#pragma omp parallel for private(color)
    for (int n = 0, Q = W * H; n < Q; ++n) {
      bool const b = (::generate_random_number_for_single_thread() <
                      kLargeStepProb);  // 一様分布からサンプリングするか？
      sample.init(b);

      float_t const x = sample.next();
      float_t const y = sample.next();
      newColor.first =
          static_cast<std::size_t>(y * H) * W + static_cast<std::size_t>(x * W);

      scene.EvaluteImportanceDrivenPathTracing(x, y, &newColor.second, sample);

      float_t const a =
          std::min(CIE_XYZ_Color::value_type(1),
                   newColor.second[1] / oldColor.second[1]);  // 受理確率
      film[oldColor.first] +=
          oldColor.second *
          ((1 - a) / (oldColor.second[1] * c + kLargeStepProb));
      film[newColor.first] +=
          newColor.second *
          ((b + a) / (newColor.second[1] * c + kLargeStepProb));

      if (::generate_random_number_for_single_thread() < a) {
        // 受理
        sample.accept();
        oldColor = newColor;
      }
    }

    // 定期的な保存
    {
      std::clock_t const current_time = std::clock() - begin_time;
      std::clock_t const time = static_cast<std::clock_t>(
          current_time / double(CLOCKS_PER_SEC) + 0.5);

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %+6u[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60, spp);

      // save
      bool const is_finished = IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += nIntervalTime;

        ::_stprintf_s(filename, MAX_PATH,
                      _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                          _TEXT(" %06uspp %012lu avg=%g") _TEXT(".pfm"),
                      GetOutdirPath().c_str(), time / (60 * 60),
                      (time / 60) % 60, time % 60, spp,
                      scene.n_FindIntersection, image_contribution_luminance);

        film.Save(filename, spp);

        if (is_finished || (time >= nRenderingTime)) {
          FinishRendering();
          return;
        }

        begin_time = std::clock() - current_time;  // 時計を戻す
      }
    }
  }
}

void PathTracerWithGoWithTheWinnersStrategy::run() {
  ::_ftprintf_s(
      stderr, _TEXT("** Path Tracing with Go with the Winners strategy **\n"));

  Scene const &scene = Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  ImageFilm film(W, H);

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  // std::vector<std::mt19937_64> randoms(omp_get_max_threads());
  RSample sample;

  // 計測開始
  int const nRenderingTime = GetRenderingTime();
  int const nIntervalTime = GetIntervalTime();
  int saveTime = nIntervalTime;
  std::clock_t begin_time = std::clock();

  //
  // Stage1: Computation.
  //
  CIE_XYZ_Color color;

  ::_ftprintf_s(stderr, _TEXT(">>>> Computation <<<<\n"));
  scene.n_FindIntersection = 0;
  for (std::size_t spp = 1;; ++spp) {
    //#pragma omp parallel for private(color)
    for (int y = 0; y < H; ++y) {
      // std::mt19937_64 & random = randoms[omp_get_thread_num()];
      for (int x = 0; x < W; ++x) {
        scene.EvalutePathTracingWithGoWithTheWinnersStrategy(
            (x + sample.next()) / W, (y + sample.next()) / H, &color, sample);
        film[y * W + x] += color;
      }
    }

    // 定期的な保存
    {
      std::clock_t const current_time = std::clock() - begin_time;
      std::clock_t const time = static_cast<std::clock_t>(
          current_time / double(CLOCKS_PER_SEC) + 0.5);

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %+6u[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60, spp);

      // save
      bool const is_finished = IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += nIntervalTime;

        ::_stprintf_s(
            filename, MAX_PATH, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                    _TEXT(" %06uspp %012lu") _TEXT(".pfm"),
            GetOutdirPath().c_str(), time / (60 * 60), (time / 60) % 60,
            time % 60, spp, scene.n_FindIntersection);

        film.Save(filename, spp);

        if (is_finished || (time >= nRenderingTime)) {
          FinishRendering();
          break;
        }

        begin_time = std::clock() - current_time;  // 時計を戻す
      }
    }
  }
}

void ImportanceDrivenPathTracerWithGoWithTheWinnersStrategy::run() {
  ::_ftprintf_s(stderr, _TEXT(
                            "** Importance Driven Path Tracing with Go with "
                            "the Winners strategy **\n"));

  Scene const &scene = Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  ImageFilm film(W, H);

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  // std::vector<std::mt19937_64> randoms(omp_get_max_threads());
  RSample sample;

  // 計測開始
  int const nRenderingTime = GetRenderingTime();
  int const nIntervalTime = GetIntervalTime();
  int saveTime = nIntervalTime;
  std::clock_t begin_time = std::clock();

  ::_ftprintf_s(stderr, _TEXT(">>>> Build Importance Map <<<<\n"));
  {
    Scene::GetInstance().CreateLightImportanceMap(
        ::get_instance_of_random_number_generator_for_single_thread());
  }

  //
  // Stage1: Computation.
  //
  CIE_XYZ_Color color;

  ::_ftprintf_s(stderr, _TEXT(">>>> Computation <<<<\n"));
  scene.n_FindIntersection = 0;
  for (std::size_t spp = 1;; ++spp) {
    //#pragma omp parallel for private(color)
    for (int y = 0; y < H; ++y) {
      // std::mt19937_64 & random = randoms[omp_get_thread_num()];
      for (int x = 0; x < W; ++x) {
        scene.EvaluteImportanceDrivenPathTracingWithGoWithTheWinnersStrategy(
            (x + sample.next()) / W, (y + sample.next()) / H, &color, sample);
        film[y * W + x] += color;
      }
    }

    // 定期的な保存
    {
      std::clock_t const current_time = std::clock() - begin_time;
      std::clock_t const time = static_cast<std::clock_t>(
          current_time / double(CLOCKS_PER_SEC) + 0.5);

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %+6u[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60, spp);

      // save
      bool const is_finished = IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += nIntervalTime;

        ::_stprintf_s(
            filename, MAX_PATH, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                    _TEXT(" %06uspp %012lu") _TEXT(".pfm"),
            GetOutdirPath().c_str(), time / (60 * 60), (time / 60) % 60,
            time % 60, spp, scene.n_FindIntersection);

        film.Save(filename, spp);

        if (is_finished || (time >= nRenderingTime)) {
          FinishRendering();
          break;
        }

        begin_time = std::clock() - current_time;  // 時計を戻す
      }
    }
  }
}

void PathTracerWithGoWithTheWinnersStrategy_MetropolisLightTransport::run() {
  ::_ftprintf_s(stderr, _TEXT(
                            "** Path Tracing with Go with the Winners "
                            "strategy & Metropolis Light Transport **\n"));

  Scene const &scene = Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  ImageFilm film(W, H);

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  // std::vector<std::mt19937_64> randoms(omp_get_max_threads());
  MSample sample;
  std::size_t oldIndex, newIndex;
  CIE_XYZ_Color oldColor, newColor;

  // 計測開始
  int const nRenderingTime = GetRenderingTime();
  int const nIntervalTime = GetIntervalTime();
  int saveTime = nIntervalTime;
  std::clock_t begin_time = std::clock();

  // 初期化
  ::_ftprintf_s(stderr, _TEXT(">>>> Initialization <<<<\n"));
  float_t image_contribution_luminance = float_t(0);
  {
    static std::size_t const R = 10000;

    std::vector<MSample> samples(R);
    std::vector<CIE_XYZ_Color> colors(R);
    std::vector<std::size_t> indices(R);
    std::vector<float_t> cumulative_distribution_function(R + 1);
    cumulative_distribution_function[0] = float_t(0);

    for (std::size_t i = 0; i < R; ++i) {
      samples[i].init(true);

      float_t const x = sample.next();
      float_t const y = sample.next();
      indices[i] =
          static_cast<std::size_t>(y * H) * W + static_cast<std::size_t>(x * W);

      scene.EvalutePathTracingWithGoWithTheWinnersStrategy(x, y, &colors[i],
                                                           sample);

      cumulative_distribution_function[i + 1] =
          cumulative_distribution_function[i] + colors[i][1];
    }

    image_contribution_luminance = cumulative_distribution_function.back() / R;
    if (image_contribution_luminance <= float_t(0)) {
      FinishRendering();
      return;
    }

    {
      // float_t const xi = (c + ::generate_random_number_for_single_thread()) /
      // XYZ_CONFIG_kSystemCount;
      float_t const xi = ::generate_random_number_for_single_thread();
      float_t const value = xi * cumulative_distribution_function.back();
      std::vector<float_t>::const_iterator const it =
          std::upper_bound(cumulative_distribution_function.begin(),
                           cumulative_distribution_function.end(), value);

      assert(it != cumulative_distribution_function.begin());
      assert(it != cumulative_distribution_function.end());

      std::size_t const uIndex =
          it - cumulative_distribution_function.begin() - 1;
      sample = samples[uIndex];
      sample.accept();  // new を old に設定
      oldIndex = indices[uIndex];
      oldColor = colors[uIndex];
    }
  }
  float_t const c = hi::rcp(image_contribution_luminance);

  static float_t const kLargeStepProb =
      float_t(0.5);  // ラージステップを行う確率

  ::_ftprintf_s(stderr, _TEXT(">>>> Computation <<<<\n"));
  scene.n_FindIntersection = 0;
  for (std::size_t spp = 1;; ++spp) {
    //#pragma omp parallel for private(color)
    for (int n = 0, Q = W * H; n < Q; ++n) {
      bool const b = (::generate_random_number_for_single_thread() <
                      kLargeStepProb);  // 一様分布からサンプリングするか？
      sample.init(b);

      float_t const x = sample.next();
      float_t const y = sample.next();
      newIndex =
          static_cast<std::size_t>(y * H) * W + static_cast<std::size_t>(x * W);

      scene.EvalutePathTracingWithGoWithTheWinnersStrategy(x, y, &newColor,
                                                           sample);

      float_t const a = std::min(CIE_XYZ_Color::value_type(1),
                                 newColor[1] / oldColor[1]);  // 受理確率
      film[oldIndex] +=
          oldColor * ((1 - a) / (oldColor[1] * c + kLargeStepProb));
      film[newIndex] +=
          newColor * ((b + a) / (newColor[1] * c + kLargeStepProb));

      if (::generate_random_number_for_single_thread() < a) {
        // 受理
        sample.accept();
        oldIndex = newIndex;
        oldColor = newColor;
      }
    }

    // 定期的な保存
    {
      std::clock_t const current_time = std::clock() - begin_time;
      std::clock_t const time = static_cast<std::clock_t>(
          current_time / double(CLOCKS_PER_SEC) + 0.5);

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %+6u[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60, spp);

      // save
      bool const is_finished = IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += nIntervalTime;

        ::_stprintf_s(filename, MAX_PATH,
                      _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                          _TEXT(" %06uspp %012lu avg=%g") _TEXT(".pfm"),
                      GetOutdirPath().c_str(), time / (60 * 60),
                      (time / 60) % 60, time % 60, spp,
                      scene.n_FindIntersection, image_contribution_luminance);

        film.Save(filename, spp);

        if (is_finished || (time >= nRenderingTime)) {
          FinishRendering();
          return;
        }

        begin_time = std::clock() - current_time;  // 時計を戻す
      }
    }
  }
}

void ImportanceDrivenPathTracerWithGoWithTheWinnersStrategy_MetropolisLightTransport::
    run() {
  ::_ftprintf_s(stderr,
                _TEXT(
                    "** Importance Driven Path Tracing with Go with the "
                    "Winners strategy & Metropolis Light Transport **\n"));

  Scene const &scene = Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  ImageFilm film(W, H);

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  // std::vector<std::mt19937_64> randoms(omp_get_max_threads());
  MSample sample;
  std::size_t oldIndex, newIndex;
  CIE_XYZ_Color oldColor, newColor;

  // 計測開始
  int const nRenderingTime = GetRenderingTime();
  int const nIntervalTime = GetIntervalTime();
  int saveTime = nIntervalTime;
  std::clock_t begin_time = std::clock();

  ::_ftprintf_s(stderr, _TEXT(">>>> Build Importance Map <<<<\n"));
  {
    Scene::GetInstance().CreateLightImportanceMap(
        ::get_instance_of_random_number_generator_for_single_thread());
  }

  // 初期化
  ::_ftprintf_s(stderr, _TEXT(">>>> Initialization <<<<\n"));
  float_t image_contribution_luminance = float_t(0);
  {
    static std::size_t const R = 10000;

    std::vector<MSample> samples(R);
    std::vector<CIE_XYZ_Color> colors(R);
    std::vector<std::size_t> indices(R);
    std::vector<float_t> cumulative_distribution_function(R + 1);
    cumulative_distribution_function[0] = float_t(0);

    for (std::size_t i = 0; i < R; ++i) {
      samples[i].init(true);

      float_t const x = sample.next();
      float_t const y = sample.next();
      indices[i] =
          static_cast<std::size_t>(y * H) * W + static_cast<std::size_t>(x * W);

      scene.EvaluteImportanceDrivenPathTracingWithGoWithTheWinnersStrategy(
          x, y, &colors[i], sample);

      cumulative_distribution_function[i + 1] =
          cumulative_distribution_function[i] + colors[i][1];
    }

    image_contribution_luminance = cumulative_distribution_function.back() / R;
    if (image_contribution_luminance <= float_t(0)) {
      FinishRendering();
      return;
    }

    {
      // float_t const xi = (c + ::generate_random_number_for_single_thread()) /
      // XYZ_CONFIG_kSystemCount;
      float_t const xi = ::generate_random_number_for_single_thread();
      float_t const value = xi * cumulative_distribution_function.back();
      std::vector<float_t>::const_iterator const it =
          std::upper_bound(cumulative_distribution_function.begin(),
                           cumulative_distribution_function.end(), value);

      assert(it != cumulative_distribution_function.begin());
      assert(it != cumulative_distribution_function.end());

      std::size_t const uIndex =
          it - cumulative_distribution_function.begin() - 1;
      sample = samples[uIndex];
      sample.accept();  // new を old に設定
      oldIndex = indices[uIndex];
      oldColor = colors[uIndex];
    }
  }
  float_t const c = hi::rcp(image_contribution_luminance);

  static float_t const kLargeStepProb =
      float_t(0.5);  // ラージステップを行う確率

  ::_ftprintf_s(stderr, _TEXT(">>>> Computation <<<<\n"));
  scene.n_FindIntersection = 0;
  for (std::size_t spp = 1;; ++spp) {
    //#pragma omp parallel for private(color)
    for (int n = 0, Q = W * H; n < Q; ++n) {
      bool const b = (::generate_random_number_for_single_thread() <
                      kLargeStepProb);  // 一様分布からサンプリングするか？
      sample.init(b);

      float_t const x = sample.next();
      float_t const y = sample.next();
      newIndex =
          static_cast<std::size_t>(y * H) * W + static_cast<std::size_t>(x * W);

      scene.EvaluteImportanceDrivenPathTracingWithGoWithTheWinnersStrategy(
          x, y, &newColor, sample);

      float_t const a = std::min(CIE_XYZ_Color::value_type(1),
                                 newColor[1] / oldColor[1]);  // 受理確率
      film[oldIndex] +=
          oldColor * ((1 - a) / (oldColor[1] * c + kLargeStepProb));
      film[newIndex] +=
          newColor * ((b + a) / (newColor[1] * c + kLargeStepProb));

      if (::generate_random_number_for_single_thread() < a) {
        // 受理
        sample.accept();
        oldIndex = newIndex;
        oldColor = newColor;
      }
    }

    // 定期的な保存
    {
      std::clock_t const current_time = std::clock() - begin_time;
      std::clock_t const time = static_cast<std::clock_t>(
          current_time / double(CLOCKS_PER_SEC) + 0.5);

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %+6u[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60, spp);

      // save
      bool const is_finished = IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += nIntervalTime;

        ::_stprintf_s(filename, MAX_PATH,
                      _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                          _TEXT(" %06uspp %012lu avg=%g") _TEXT(".pfm"),
                      GetOutdirPath().c_str(), time / (60 * 60),
                      (time / 60) % 60, time % 60, spp,
                      scene.n_FindIntersection, image_contribution_luminance);

        film.Save(filename, spp);

        if (is_finished || (time >= nRenderingTime)) {
          FinishRendering();
          return;
        }

        begin_time = std::clock() - current_time;  // 時計を戻す
      }
    }
  }
}

void BidirectionalPathTracer::run() {
  ::_ftprintf_s(stderr, _TEXT("** Bidirectional Path Tracing **\n"));

  Scene const &scene = Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  ImageFilm film(W, H);

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  // std::vector<std::mt19937_64> randoms(omp_get_max_threads());
  std::vector<pixel_descriptor_t> colors(XYZ_CONFIG_kMaxRandomWalkDepth);
  std::vector<PathVertex> lightsPathVertices(XYZ_CONFIG_kMaxRandomWalkDepth +
                                             2);
  std::vector<PathVertex> theEyePathVertices(XYZ_CONFIG_kMaxRandomWalkDepth +
                                             2);
  RSample sample;

  // 計測開始
  int const nRenderingTime = GetRenderingTime();
  int const nIntervalTime = GetIntervalTime();
  int saveTime = nIntervalTime;
  std::clock_t begin_time = std::clock();

  //
  // Stage1: Computation.
  //
  ::_ftprintf_s(stderr, _TEXT(">>>> Computation <<<<\n"));
  scene.n_FindIntersection = 0;
  for (std::size_t spp = 1;; ++spp) {
    //#pragma omp parallel for private(color)
    for (int y = 0; y < H; ++y) {
      // std::mt19937_64 & random = randoms[omp_get_thread_num()];
      for (int x = 0; x < W; ++x) {
        // 評価値を初期化する．
        colors[0].first = y * W + x;
        for (std::size_t i = 1; i < XYZ_CONFIG_kMaxRandomWalkDepth; ++i) {
          colors[i].first = ~0U;
        }

        scene.EvaluateBidirectionalPathTracing(
            (x + sample.next()) / W, (y + sample.next()) / H, &colors,
            lightsPathVertices, theEyePathVertices, sample);

        // 評価値を蓄える．
        for (std::size_t i = 0; i < XYZ_CONFIG_kMaxRandomWalkDepth; ++i) {
          if (~0U != colors[i].first) {
            film[colors[i].first] += colors[i].second;
          }
        }
      }
    }

    // 定期的な保存
    {
      std::clock_t const current_time = std::clock() - begin_time;
      std::clock_t const time = static_cast<std::clock_t>(
          current_time / double(CLOCKS_PER_SEC) + 0.5);

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %+6u[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60, spp);

      // save
      bool const is_finished = IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += nIntervalTime;

        ::_stprintf_s(
            filename, MAX_PATH, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                    _TEXT(" %06uspp %012lu") _TEXT(".pfm"),
            GetOutdirPath().c_str(), time / (60 * 60), (time / 60) % 60,
            time % 60, spp, scene.n_FindIntersection);

        film.Save(filename, spp);

        if (is_finished || (time >= nRenderingTime)) {
          FinishRendering();
          break;
        }

        begin_time = std::clock() - current_time;  // 時計を戻す
      }
    }
  }
}

void ImportanceDrivenBidirectionalPathTracer::run() {
  ::_ftprintf_s(stderr,
                _TEXT("** Importance Driven Bidirectional Path Tracing **\n"));

  Scene const &scene = Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  ImageFilm film(W, H);

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  // std::vector<std::mt19937_64> randoms(omp_get_max_threads());
  std::vector<pixel_descriptor_t> colors(XYZ_CONFIG_kMaxRandomWalkDepth);
  std::vector<PathVertex> lightsPathVertices(XYZ_CONFIG_kMaxRandomWalkDepth +
                                             2);
  std::vector<PathVertex> theEyePathVertices(XYZ_CONFIG_kMaxRandomWalkDepth +
                                             2);
  RSample sample;

  // 計測開始
  int const nRenderingTime = GetRenderingTime();
  int const nIntervalTime = GetIntervalTime();
  int saveTime = nIntervalTime;
  std::clock_t begin_time = std::clock();

  ::_ftprintf_s(stderr, _TEXT(">>>> Build Importance Map <<<<\n"));
  {
    Scene::GetInstance().CreateLightImportanceMap(
        get_instance_of_random_number_generator_for_single_thread());
  }

  ::_ftprintf_s(stderr, _TEXT(">>>> Computation <<<<\n"));
  scene.n_FindIntersection = 0;
  for (std::size_t spp = 1;; ++spp) {
    //#pragma omp parallel for private(color)
    for (int y = 0; y < H; ++y) {
      // std::mt19937_64 & random = randoms[omp_get_thread_num()];
      for (int x = 0; x < W; ++x) {
        // 評価値を初期化する．
        colors[0].first = y * W + x;
        for (std::size_t i = 1; i < XYZ_CONFIG_kMaxRandomWalkDepth; ++i) {
          colors[i].first = ~0U;
        }

        scene.EvaluateImportanceDrivenBidirectionalPathTracing(
            (x + sample.next()) / W, (y + sample.next()) / H, &colors,
            lightsPathVertices, theEyePathVertices, sample);

        // 評価値を蓄える．
        for (std::size_t i = 0; i < XYZ_CONFIG_kMaxRandomWalkDepth; ++i) {
          if (~0U != colors[i].first) {
            film[colors[i].first] += colors[i].second;
          }
        }
      }
    }

    // 定期的な保存
    {
      std::clock_t const current_time = std::clock() - begin_time;
      std::clock_t const time = static_cast<std::clock_t>(
          current_time / double(CLOCKS_PER_SEC) + 0.5);

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %+6u[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60, spp);

      // save
      bool const is_finished = IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += nIntervalTime;

        ::_stprintf_s(
            filename, MAX_PATH, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                    _TEXT(" %06uspp %012lu") _TEXT(".pfm"),
            GetOutdirPath().c_str(), time / (60 * 60), (time / 60) % 60,
            time % 60, spp, scene.n_FindIntersection);

        film.Save(filename, spp);

        if (is_finished || (time >= nRenderingTime)) {
          FinishRendering();
          break;
        }

        begin_time = std::clock() - current_time;  // 時計を戻す
      }
    }
  }
}

void BidirectionalPathTracer_MetropolisLightTransport::run() {
  ::_ftprintf_s(
      stderr,
      _TEXT("** Bidirectional Path Tracing & Metropolis Light Transport **\n"));

  Scene const &scene = Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  ImageFilm film(W, H);

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  // std::vector<std::mt19937_64> randoms(omp_get_max_threads());
  std::vector<pixel_descriptor_t> oldColors(XYZ_CONFIG_kMaxRandomWalkDepth);
  std::vector<pixel_descriptor_t> newColors(XYZ_CONFIG_kMaxRandomWalkDepth);
  float_t oldValue = 0;
  float_t newValue = 0;
  std::vector<PathVertex> lightsPathVertices(XYZ_CONFIG_kMaxRandomWalkDepth +
                                             2);
  std::vector<PathVertex> theEyePathVertices(XYZ_CONFIG_kMaxRandomWalkDepth +
                                             2);
  MSample sample;

  // 計測開始
  int const nRenderingTime = GetRenderingTime();
  int const nIntervalTime = GetIntervalTime();
  int saveTime = nIntervalTime;
  std::clock_t begin_time = std::clock();

  // 初期化
  ::_ftprintf_s(stderr, _TEXT(">>>> Initialization <<<<\n"));
  float_t image_contribution_luminance = float_t(0);
  {
    static std::size_t const R = 10000;

    std::vector<MSample> samples(R);
    std::vector<std::vector<pixel_descriptor_t>> colors(R);
    std::vector<float_t> cumulative_distribution_function(R + 1);
    cumulative_distribution_function[0] = float_t(0);

    for (std::size_t i = 0; i < R; ++i) {
      samples[i].init(true);

      colors[i].resize(XYZ_CONFIG_kMaxRandomWalkDepth);

      float_t const x = sample.next();
      float_t const y = sample.next();

      // 評価値を初期化する．
      colors[i][0].first =
          static_cast<std::size_t>(y * H) * W + static_cast<std::size_t>(x * W);
      for (std::size_t k = 1; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        colors[i][k].first = ~0U;
      }

      scene.EvaluateBidirectionalPathTracing(
          x, y, &colors[i], lightsPathVertices, theEyePathVertices, sample);

      // 評価値を蓄える．
      float_t evaluated_value = float_t();
      for (std::size_t k = 0; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        if (~0U != colors[i][k].first) {
          evaluated_value += colors[i][k].second[1];
        }
      }

      cumulative_distribution_function[i + 1] =
          cumulative_distribution_function[i] + evaluated_value;
    }

    image_contribution_luminance = cumulative_distribution_function.back() / R;
    if (image_contribution_luminance <= float_t(0)) {
      FinishRendering();
      return;
    }

    {
      // float_t const xi = (c + ::generate_random_number_for_single_thread()) /
      // XYZ_CONFIG_kSystemCount;
      float_t const xi = ::generate_random_number_for_single_thread();
      float_t const value = xi * cumulative_distribution_function.back();
      std::vector<float_t>::const_iterator const it =
          std::upper_bound(cumulative_distribution_function.begin(),
                           cumulative_distribution_function.end(), value);

      assert(it != cumulative_distribution_function.begin());
      assert(it != cumulative_distribution_function.end());

      std::size_t const uIndex =
          it - cumulative_distribution_function.begin() - 1;
      sample = samples[uIndex];
      sample.accept();  // new を old に設定
      std::swap(oldColors, colors[uIndex]);
      oldValue = cumulative_distribution_function[uIndex + 1] -
                 cumulative_distribution_function[uIndex];
    }
  }
  float_t const c = hi::rcp(image_contribution_luminance);

  static float_t const kLargeStepProb =
      float_t(0.5);  // ラージステップを行う確率

  ::_ftprintf_s(stderr, _TEXT(">>>> Computation <<<<\n"));
  scene.n_FindIntersection = 0;
  for (std::size_t spp = 1;; ++spp) {
    //#pragma omp parallel for private(color)
    for (int n = 0, Q = W * H; n < Q; ++n) {
      bool const b = (::generate_random_number_for_single_thread() <
                      kLargeStepProb);  // 一様分布からサンプリングするか？
      sample.init(b);

      float_t const x = sample.next();
      float_t const y = sample.next();

      // 評価値を初期化する．
      newColors[0].first =
          static_cast<std::size_t>(y * H) * W + static_cast<std::size_t>(x * W);
      for (std::size_t k = 1; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        newColors[k].first = ~0U;
      }

      scene.EvaluateBidirectionalPathTracing(
          x, y, &newColors, lightsPathVertices, theEyePathVertices, sample);

      newValue = float_t();
      for (std::size_t k = 0; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        if (~0U != newColors[k].first) {
          newValue += newColors[k].second[1];
        }
      }

      float_t const a = std::min(CIE_XYZ_Color::value_type(1),
                                 newValue / oldValue);  // 受理確率

      // 評価値を蓄える．
      float_t const oldScale = (1 - a) / (oldValue * c + kLargeStepProb);
      float_t const newScale = (b + a) / (newValue * c + kLargeStepProb);
      for (std::size_t k = 0; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        if (~0U != oldColors[k].first) {
          film[oldColors[k].first] += oldColors[k].second * oldScale;
        }
        if (~0U != newColors[k].first) {
          film[newColors[k].first] += newColors[k].second * newScale;
        }
      }

      if (::generate_random_number_for_single_thread() < a) {
        // 受理
        sample.accept();
        std::swap(oldColors, newColors);
        oldValue = newValue;
      }
    }

    // 定期的な保存
    {
      std::clock_t const current_time = std::clock() - begin_time;
      std::clock_t const time = static_cast<std::clock_t>(
          current_time / double(CLOCKS_PER_SEC) + 0.5);

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %+6u[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60, spp);

      // save
      bool const is_finished = IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += nIntervalTime;

        ::_stprintf_s(filename, MAX_PATH,
                      _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                          _TEXT(" %06uspp %012lu avg=%g") _TEXT(".pfm"),
                      GetOutdirPath().c_str(), time / (60 * 60),
                      (time / 60) % 60, time % 60, spp,
                      scene.n_FindIntersection, image_contribution_luminance);

        film.Save(filename, spp);

        if (is_finished || (time >= nRenderingTime)) {
          FinishRendering();
          return;
        }

        begin_time = std::clock() - current_time;  // 時計を戻す
      }
    }
  }
}

void ImportanceDrivenBidirectionalPathTracer_MetropolisLightTransport::run() {
  ::_ftprintf_s(stderr, _TEXT(
                            "** Importance Driven Bidirectional Path Tracing "
                            "& Metropolis Light Transport **\n"));

  Scene const &scene = Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  ImageFilm film(W, H);

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  // std::vector<std::mt19937_64> randoms(omp_get_max_threads());
  std::vector<pixel_descriptor_t> oldColors(XYZ_CONFIG_kMaxRandomWalkDepth);
  std::vector<pixel_descriptor_t> newColors(XYZ_CONFIG_kMaxRandomWalkDepth);
  float_t oldValue = 0;
  float_t newValue = 0;
  std::vector<PathVertex> lightsPathVertices(XYZ_CONFIG_kMaxRandomWalkDepth +
                                             2);
  std::vector<PathVertex> theEyePathVertices(XYZ_CONFIG_kMaxRandomWalkDepth +
                                             2);
  MSample sample;

  // 計測開始
  int const nRenderingTime = GetRenderingTime();
  int const nIntervalTime = GetIntervalTime();
  int saveTime = nIntervalTime;
  std::clock_t begin_time = std::clock();

  ::_ftprintf_s(stderr, _TEXT(">>>> Build Importance Map <<<<\n"));
  {
    Scene::GetInstance().CreateLightImportanceMap(
        get_instance_of_random_number_generator_for_single_thread());
  }

  // 初期化
  ::_ftprintf_s(stderr, _TEXT(">>>> Initialization <<<<\n"));
  float_t image_contribution_luminance = float_t(0);
  {
    static std::size_t const R = 10000;

    std::vector<MSample> samples(R);
    std::vector<std::vector<pixel_descriptor_t>> colors(R);
    std::vector<float_t> cumulative_distribution_function(R + 1);
    cumulative_distribution_function[0] = float_t(0);

    for (std::size_t i = 0; i < R; ++i) {
      samples[i].init(true);

      colors[i].resize(XYZ_CONFIG_kMaxRandomWalkDepth);

      float_t const x = sample.next();
      float_t const y = sample.next();

      // 評価値を初期化する．
      colors[i][0].first =
          static_cast<std::size_t>(y * H) * W + static_cast<std::size_t>(x * W);
      for (std::size_t k = 1; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        colors[i][k].first = ~0U;
      }

      scene.EvaluateBidirectionalPathTracing(
          x, y, &colors[i], lightsPathVertices, theEyePathVertices, sample);

      // 評価値を蓄える．
      float_t evaluated_value = float_t();
      for (std::size_t k = 0; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        if (~0U != colors[i][k].first) {
          evaluated_value += colors[i][k].second[1];
        }
      }

      cumulative_distribution_function[i + 1] =
          cumulative_distribution_function[i] + evaluated_value;
    }

    image_contribution_luminance = cumulative_distribution_function.back() / R;
    if (image_contribution_luminance <= float_t(0)) {
      FinishRendering();
      return;
    }

    {
      // float_t const xi = (c + ::generate_random_number_for_single_thread()) /
      // XYZ_CONFIG_kSystemCount;
      float_t const xi = ::generate_random_number_for_single_thread();
      float_t const value = xi * cumulative_distribution_function.back();
      std::vector<float_t>::const_iterator const it =
          std::upper_bound(cumulative_distribution_function.begin(),
                           cumulative_distribution_function.end(), value);

      assert(it != cumulative_distribution_function.begin());
      assert(it != cumulative_distribution_function.end());

      std::size_t const uIndex =
          it - cumulative_distribution_function.begin() - 1;
      sample = samples[uIndex];
      sample.accept();  // new を old に設定
      std::swap(oldColors, colors[uIndex]);
      oldValue = cumulative_distribution_function[uIndex + 1] -
                 cumulative_distribution_function[uIndex];
    }
  }
  float_t const c = hi::rcp(image_contribution_luminance);

  static float_t const kLargeStepProb =
      float_t(0.5);  // ラージステップを行う確率

  ::_ftprintf_s(stderr, _TEXT(">>>> Computation <<<<\n"));
  scene.n_FindIntersection = 0;
  for (std::size_t spp = 1;; ++spp) {
    //#pragma omp parallel for private(color)
    for (int n = 0, Q = W * H; n < Q; ++n) {
      bool const b = (::generate_random_number_for_single_thread() <
                      kLargeStepProb);  // 一様分布からサンプリングするか？
      sample.init(b);

      float_t const x = sample.next();
      float_t const y = sample.next();

      // 評価値を初期化する．
      newColors[0].first =
          static_cast<std::size_t>(y * H) * W + static_cast<std::size_t>(x * W);
      for (std::size_t k = 1; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        newColors[k].first = ~0U;
      }

      scene.EvaluateImportanceDrivenBidirectionalPathTracing(
          x, y, &newColors, lightsPathVertices, theEyePathVertices, sample);

      newValue = float_t();
      for (std::size_t k = 0; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        if (~0U != newColors[k].first) {
          newValue += newColors[k].second[1];
        }
      }

      float_t const a = std::min(CIE_XYZ_Color::value_type(1),
                                 newValue / oldValue);  // 受理確率

      // 評価値を蓄える．
      float_t const oldScale = (1 - a) / (oldValue * c + kLargeStepProb);
      float_t const newScale = (b + a) / (newValue * c + kLargeStepProb);
      for (std::size_t k = 0; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        if (~0U != oldColors[k].first) {
          film[oldColors[k].first] += oldColors[k].second * oldScale;
        }
        if (~0U != newColors[k].first) {
          film[newColors[k].first] += newColors[k].second * newScale;
        }
      }

      if (::generate_random_number_for_single_thread() < a) {
        // 受理
        sample.accept();
        std::swap(oldColors, newColors);
        oldValue = newValue;
      }
    }

    // 定期的な保存
    {
      std::clock_t const current_time = std::clock() - begin_time;
      std::clock_t const time = static_cast<std::clock_t>(
          current_time / double(CLOCKS_PER_SEC) + 0.5);

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %+6u[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60, spp);

      // save
      bool const is_finished = IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += nIntervalTime;

        ::_stprintf_s(filename, MAX_PATH,
                      _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                          _TEXT(" %06uspp %012lu avg=%g") _TEXT(".pfm"),
                      GetOutdirPath().c_str(), time / (60 * 60),
                      (time / 60) % 60, time % 60, spp,
                      scene.n_FindIntersection, image_contribution_luminance);

        film.Save(filename, spp);

        if (is_finished || (time >= nRenderingTime)) {
          FinishRendering();
          return;
        }

        begin_time = std::clock() - current_time;  // 時計を戻す
      }
    }
  }
}

void ResampledPathTracer::run() {
  ::_ftprintf_s(stderr, _TEXT("** Resampled Path Tracing **\n"));

  Scene const &scene = Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  ImageFilm film(W, H);

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  // std::vector<std::mt19937_64> randoms(omp_get_max_threads());
  RSample sample;

  std::vector<PathVertex> lightPathVertices(CONFIG_BLOCK_SIZE *
                                            CONFIG_BLOCK_SIZE);
  std::vector<Particle> ray(CONFIG_BLOCK_SIZE * CONFIG_BLOCK_SIZE);
  std::vector<Particle> tmp(CONFIG_BLOCK_SIZE * CONFIG_BLOCK_SIZE);
  std::vector<float_t> cdf(CONFIG_BLOCK_SIZE * CONFIG_BLOCK_SIZE + 1);

  // 計測開始
  int const nRenderingTime = GetRenderingTime();
  int const nIntervalTime = GetIntervalTime();
  int saveTime = nIntervalTime;
  std::clock_t begin_time = std::clock();

  //
  // Stage1: Computation.
  //
  std::vector<CIE_XYZ_Color> color(CONFIG_BLOCK_SIZE * CONFIG_BLOCK_SIZE);

  ::_ftprintf_s(stderr, _TEXT(">>>> Computation <<<<\n"));
  scene.n_FindIntersection = 0;
  for (std::size_t spp = 1;; ++spp) {
    //#pragma omp parallel for private(color)
    for (int y = 0; y < H; y += CONFIG_BLOCK_SIZE) {
      // std::mt19937_64 & random = randoms[omp_get_thread_num()];
      for (int x = 0; x < W; x += CONFIG_BLOCK_SIZE) {
        scene.EvaluteResampledPathTracing(x, y, &color, lightPathVertices, ray,
                                          tmp, cdf, sample);

        for (std::size_t i = 0, dy = 0; dy < CONFIG_BLOCK_SIZE; ++dy) {
          for (std::size_t dx = 0; dx < CONFIG_BLOCK_SIZE; ++dx, ++i) {
            film[(y + dy) * W + (x + dx)] += color[i];
          }
        }
      }
    }

    // 定期的な保存
    {
      std::clock_t const current_time = std::clock() - begin_time;
      std::clock_t const time = static_cast<std::clock_t>(
          current_time / double(CLOCKS_PER_SEC) + 0.5);

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %+6u[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60, spp);

      // save
      bool const is_finished = IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += nIntervalTime;

        ::_stprintf_s(
            filename, MAX_PATH, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                    _TEXT(" %06uspp %012lu") _TEXT(".pfm"),
            GetOutdirPath().c_str(), time / (60 * 60), (time / 60) % 60,
            time % 60, spp, scene.n_FindIntersection);

        film.Save(filename, spp);

        if (is_finished || (time >= nRenderingTime)) {
          FinishRendering();
          break;
        }

        begin_time = std::clock() - current_time;  // 時計を戻す
      }
    }
  }
}

void MetropolisLightTransportTestbed::run() {
  ::_ftprintf_s(stderr, _TEXT("** Metropolis Light Transport Testbed **\n"));
  static std::size_t const knInitialSample = 10000;  // 初期サンプル数

  Scene const &scene = Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();
  int const Q = W * H;

  ImageFilm film(W, H);    // 結果
  ImageFilm filmU1(W, H);  // u(x)
  ImageFilm filmP1(W, H);  // p(x)
  ImageFilm filmU2(W, H);  // u^2(x)
  ImageFilm filmP2(W, H);  // p^2(x)

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  // std::vector<std::mt19937_64> randoms(omp_get_max_threads());
  std::vector<pixel_descriptor_t> oldColors(XYZ_CONFIG_kMaxRandomWalkDepth);
  std::vector<pixel_descriptor_t> newColors(XYZ_CONFIG_kMaxRandomWalkDepth);
  float_t oldValue = 0;
  float_t newValue = 0;
  std::vector<PathVertex> lightsPathVertices(XYZ_CONFIG_kMaxRandomWalkDepth +
                                             2);
  std::vector<PathVertex> theEyePathVertices(XYZ_CONFIG_kMaxRandomWalkDepth +
                                             2);
  MSample sample;

  // 計測開始
  int const nRenderingTime = GetRenderingTime();
  int const nIntervalTime = GetIntervalTime();
  int saveTime = nIntervalTime;
  std::clock_t begin_time = std::clock();

  // 初期化
  ::_ftprintf_s(stderr, _TEXT(">>>> Initialization <<<<\n"));
  float_t image_contribution_luminance = float_t(0);
  {
    std::vector<MSample> samples(knInitialSample);
    std::vector<std::vector<pixel_descriptor_t>> colors(knInitialSample);
    std::vector<float_t> cumulative_distribution_function(knInitialSample + 1);
    cumulative_distribution_function[0] = float_t(0);

    for (std::size_t i = 0; i < knInitialSample; ++i) {
      samples[i].init(true);

      colors[i].resize(XYZ_CONFIG_kMaxRandomWalkDepth);

      float_t const x = sample.next();
      float_t const y = sample.next();

      // 評価値を初期化する．
      colors[i][0].first =
          static_cast<std::size_t>(y * H) * W + static_cast<std::size_t>(x * W);
      for (std::size_t k = 1; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        colors[i][k].first = ~0U;
      }

      scene.EvaluateBidirectionalPathTracing(
          x, y, &colors[i], lightsPathVertices, theEyePathVertices, sample);

      // 評価値を蓄える．
      float_t evaluated_value = float_t();
      for (std::size_t k = 0; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        if (~0U != colors[i][k].first) {
          evaluated_value += colors[i][k].second[1];
        }
      }

      cumulative_distribution_function[i + 1] =
          cumulative_distribution_function[i] + evaluated_value;
    }

    image_contribution_luminance = cumulative_distribution_function.back();
    if (image_contribution_luminance <= float_t(0)) {
      FinishRendering();
      return;
    }

    {
      // float_t const xi = (c + ::generate_random_number_for_single_thread()) /
      // XYZ_CONFIG_kSystemCount;
      float_t const xi = ::generate_random_number_for_single_thread();
      float_t const value = xi * cumulative_distribution_function.back();
      std::vector<float_t>::const_iterator const it =
          std::upper_bound(cumulative_distribution_function.begin(),
                           cumulative_distribution_function.end(), value);

      assert(it != cumulative_distribution_function.begin());
      assert(it != cumulative_distribution_function.end());

      std::size_t const uIndex =
          it - cumulative_distribution_function.begin() - 1;
      sample = samples[uIndex];
      sample.accept();  // new を old に設定
      std::swap(oldColors, colors[uIndex]);
      oldValue = cumulative_distribution_function[uIndex + 1] -
                 cumulative_distribution_function[uIndex];
    }
  }
  float_t const c =
      knInitialSample / image_contribution_luminance;  // 明るさの逆数

  static float_t const kLargeStepProb =
      float_t(0.5);  // ラージステップを行う確率

  ::_ftprintf_s(stderr, _TEXT(">>>> Computation <<<<\n"));
  scene.n_FindIntersection = 0;

  float_t c_true = 0;      // 真の明るさ
  std::size_t n_true = 0;  // 真のラージステップの回数

  for (std::size_t spp = 1;; ++spp) {
    //#pragma omp parallel for private(color)
    for (int n = 0; n < Q; ++n) {
      bool const b = (::generate_random_number_for_single_thread() <
                      kLargeStepProb);  // 一様分布からサンプリングするか？
      sample.init(b);

      float_t const x = sample.next();
      float_t const y = sample.next();

      // 評価値を初期化する．
      newColors[0].first =
          static_cast<std::size_t>(y * H) * W + static_cast<std::size_t>(x * W);
      for (std::size_t k = 1; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        newColors[k].first = ~0U;
      }

      scene.EvaluateBidirectionalPathTracing(
          x, y, &newColors, lightsPathVertices, theEyePathVertices, sample);

      newValue = float_t();
      for (std::size_t k = 0; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        if (~0U != newColors[k].first) {
          newValue += newColors[k].second[1];
        }
      }

      if (b) {
        c_true += newValue;
        ++n_true;
      }

      float_t const a = std::min(CIE_XYZ_Color::value_type(1),
                                 newValue / oldValue);  // 受理確率

      // 評価値を蓄える．
      float_t const oldScale = hi::rcp(oldValue * c + kLargeStepProb);
      float_t const newScale = hi::rcp(newValue * c + kLargeStepProb);
      float_t const oldScaleSq =
          hi::rcp(hi::square_of(oldValue * c) + hi::square_of(kLargeStepProb));
      float_t const newScaleSq =
          hi::rcp(hi::square_of(newValue * c) + hi::square_of(kLargeStepProb));
      for (std::size_t k = 0; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        if (~0U != oldColors[k].first) {
          filmP1[oldColors[k].first] +=
              oldColors[k].second * ((1 - a) * oldScale);
          filmP2[oldColors[k].first] +=
              oldColors[k].second * ((1 - a) * oldScaleSq * oldValue * c);
        }
        if (~0U != newColors[k].first) {
          filmP1[newColors[k].first] += newColors[k].second * ((a)*newScale);
          filmP2[newColors[k].first] +=
              newColors[k].second * ((a)*newScaleSq * newValue * c);
          if (b) {
            filmU1[newColors[k].first] += newColors[k].second * (newScale);
            filmU2[newColors[k].first] +=
                newColors[k].second * (newScaleSq * kLargeStepProb);
          }
        }
      }

      if (::generate_random_number_for_single_thread() < a) {
        // 受理
        sample.accept();
        std::swap(oldColors, newColors);
        oldValue = newValue;
      }
    }

    // 定期的な保存
    {
      std::clock_t const current_time = std::clock() - begin_time;
      std::clock_t const time = static_cast<std::clock_t>(
          current_time / double(CLOCKS_PER_SEC) + 0.5);

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %+6u[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60, spp);

      // save
      bool const is_finished = IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += nIntervalTime;

        // そのまま
        {
          ::_stprintf_s(filename, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                      _TEXT(" %06uspp %010lu") _TEXT("N.pfm"),
                        GetOutdirPath().c_str(), time / (60 * 60),
                        (time / 60) % 60, time % 60, spp,
                        scene.n_FindIntersection);

          for (int q = 0; q < Q; ++q) {
            film[q] = filmU1[q];
            film[q] += filmP1[q];
          }

          film.Save(filename, spp);
        }

        float_t const image_contribution_luminance_true =
            (image_contribution_luminance + c_true) /
            (knInitialSample + n_true);

        // バランス発見法
        {
          ::_stprintf_s(filename, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                      _TEXT(" %06uspp %010lu") _TEXT("B.pfm"),
                        GetOutdirPath().c_str(), time / (60 * 60),
                        (time / 60) % 60, time % 60, spp,
                        scene.n_FindIntersection);

          float_t const alpha = image_contribution_luminance_true * c;
          for (int q = 0; q < Q; ++q) {
            film[q] = filmU1[q];
            film[q] += filmP1[q] * alpha;
          }

          film.Save(filename, spp);
        }

        // 指数発見法
        {
          ::_stprintf_s(filename, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                      _TEXT(" %06uspp %010lu") _TEXT("P.pfm"),
                        GetOutdirPath().c_str(), time / (60 * 60),
                        (time / 60) % 60, time % 60, spp,
                        scene.n_FindIntersection);

          float_t const alpha = image_contribution_luminance_true * c;
          for (int q = 0; q < Q; ++q) {
            film[q] = filmU2[q];
            film[q] += filmP2[q] * alpha;
          }

          film.Save(filename, spp);
        }

        ::_ftprintf_s(stdout, _T("I, %g, %g\n"),
                      image_contribution_luminance / knInitialSample,
                      image_contribution_luminance_true);

        if (is_finished || (time >= nRenderingTime)) {
          FinishRendering();
          return;
        }

        begin_time = std::clock() - current_time;  // 時計を戻す
      }
    }
  }
}

void PathTracer_PopulationLightTransportTestbed::run() {
  ::_ftprintf_s(
      stderr,
      _TEXT("** Path Tracing & Metropolis Light Transport Testbed **\n"));
  static std::size_t const knInitialSample = 10000;  // 初期サンプル数
  static std::size_t const knPDFSample = 2;  // g(x)/I からのサンプル数

  Scene const &scene = Scene::GetInstance();
  scene.n_FindIntersection = 0;

  int const kImageWidth = scene.GetWidth();
  int const kImageHeight = scene.GetHeight();
  int const knImagePixel = kImageWidth * kImageHeight;

  ImageFilm film(kImageWidth, kImageHeight);   // 結果
  ImageFilm filmU(kImageWidth, kImageHeight);  // u(x)
  ImageFilm filmP(kImageWidth, kImageHeight);  // p(x)

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  // std::vector<std::mt19937_64> randoms(omp_get_max_threads());

  // 計測開始
  int const nRenderingTime = GetRenderingTime();
  int const nIntervalTime = GetIntervalTime();
  int saveTime = nIntervalTime;
  std::clock_t begin_time = std::clock();

  // 乱数生成器
  std::vector<MSample> old_primary_samples(knInitialSample);
  std::vector<MSample> new_primary_samples(knInitialSample);

  // 寄与ベクトル
  std::vector<pixel_descriptor_t> old_contribution_vectors(knInitialSample);
  std::vector<pixel_descriptor_t> new_contribution_vectors(knInitialSample);

  // 累積密度分布(p_{k}からp_{k+1}のサンプリング用)
  std::vector<float_t> cumulative_mass_function(knInitialSample + 1);

  // Population Annealing を反復する
  float_t image_contribution_luminance = float_t();
  for (std::size_t n = 1;; ++n) {
    //
    // p_{k=0}(x)=u(x) から X_{i} をサンプリング
    //
    cumulative_mass_function[0] = float_t();
    for (std::size_t i = 0; i < knInitialSample; ++i) {
      new_primary_samples[i].init(true);

      // 画像平面上の位置を特定
      float_t const x = new_primary_samples[i].next();
      float_t const y = new_primary_samples[i].next();

      // 画素の位置(j)を初期化
      new_contribution_vectors[i].first =
          static_cast<std::size_t>(y * kImageHeight) * kImageWidth +
          static_cast<std::size_t>(x * kImageWidth);

      // シーンをサンプリング
      scene.EvalutePathTracing(x, y, &new_contribution_vectors[i].second,
                               new_primary_samples[i]);

      // 累積密度関数を更新
      cumulative_mass_function[i + 1] =
          cumulative_mass_function[i] + new_contribution_vectors[i].second[1];
    }

    if (cumulative_mass_function.back() <= float_t()) {
      continue;  // probably empty scene
    }

    // 正確な I の計算のために寄与を累積
    image_contribution_luminance += cumulative_mass_function.back();

    // 近似的な I_k を計算
    float_t const image_contribution_luminance_k =
        image_contribution_luminance / (knInitialSample * n);

    // X_{i}~p_{k=0}(x)=u(x) の寄与をバッファに格納
    for (std::size_t i = 0; i < knInitialSample; ++i) {
      // Multiple Importance Sampling 用の重みを計算
      float_t const heuristic_weight =
          1 / (knPDFSample * new_contribution_vectors[i].second[1] /
                   image_contribution_luminance_k +
               1);

      // 寄与をバッファに貯蓄
      filmU[new_contribution_vectors[i].first] +=
          new_contribution_vectors[i].second * heuristic_weight;
    }

    // 再サンプリング
    {
      typedef std::vector<float_t>::const_iterator cmf_const_iterator;
      cmf_const_iterator const begin = cumulative_mass_function.begin();
      cmf_const_iterator const end = cumulative_mass_function.end();
      cmf_const_iterator it = begin;
      for (std::size_t i = 0; i < knInitialSample; ++i) {
        float_t const xi = (i + ::generate_random_number_for_single_thread()) /
                           knInitialSample;
        float_t const value = xi * cumulative_mass_function.back();

        it = std::upper_bound(it, end, value);  // binary search

        assert(it != begin);
        assert(it != end);

        std::size_t const resampled_index = it - begin - 1;
        old_primary_samples[i] = new_primary_samples[resampled_index];
        old_contribution_vectors[i] = new_contribution_vectors[resampled_index];
      }
      new_primary_samples.swap(old_primary_samples);  // old <=> new を入れ替え
    }

    //
    // p_{k=1}(x)=g(x)/I から X_{i} をサンプリング
    //
    for (std::size_t i = 0; i < knInitialSample; ++i) {
      // short length Metropolis-Hastings sampling
      pixel_descriptor_t &old_contribution_vector = old_contribution_vectors[i];
      pixel_descriptor_t &new_contribution_vector = new_contribution_vectors[i];

      for (std::size_t m = 0;;) {
        new_primary_samples[i].init(false);

        // 画像平面上の位置を特定
        float_t const x = new_primary_samples[i].next();
        float_t const y = new_primary_samples[i].next();

        // 画素の位置(j)を初期化
        new_contribution_vector.first =
            static_cast<std::size_t>(y * kImageHeight) * kImageWidth +
            static_cast<std::size_t>(x * kImageWidth);

        // シーンをサンプリング
        scene.EvalutePathTracing(x, y, &new_contribution_vector.second,
                                 new_primary_samples[i]);

        // Metropolis-Hastings Ratio から受理確率を計算する
        float_t const a =
            std::min(float_t(1), new_contribution_vector.second[1] /
                                     old_contribution_vector.second[1]);

        // Multiple Importance Sampling 用の重みを計算
        //   u(x) に対する重みに 1/image_contribution_luminance_k を掛けたもの
        float_t const new_heuristic_weight =
            knPDFSample / (knPDFSample * new_contribution_vector.second[1] +
                           image_contribution_luminance_k);
        float_t const old_heuristic_weight =
            knPDFSample / (knPDFSample * old_contribution_vector.second[1] +
                           image_contribution_luminance_k);

        filmP[new_contribution_vector.first] +=
            new_contribution_vector.second * ((a)*new_heuristic_weight);
        filmP[old_contribution_vector.first] +=
            old_contribution_vector.second * ((1 - a) * old_heuristic_weight);

        if (++m >= knPDFSample) {
          break;
        }

        if (::generate_random_number_for_single_thread() < a) {
          // 受理
          new_primary_samples[i].accept();
          old_contribution_vector = new_contribution_vector;
        }
      }
    }

    //
    // 定期的な保存
    //
    {
      std::clock_t const current_time = std::clock() - begin_time;
      std::clock_t const time = static_cast<std::clock_t>(
          current_time / double(CLOCKS_PER_SEC) + 0.5);
      std::size_t const the_number_of_samples = knInitialSample * n;
      float_t const samples_per_pixel =
          static_cast<float_t>(the_number_of_samples) / knImagePixel;

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %g[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60,
                    samples_per_pixel);

      // save
      bool const is_finished = IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += nIntervalTime;

        {
          ::_stprintf_s(filename, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                      _TEXT(" %gspp %010lu") _TEXT(".pfm"),
                        GetOutdirPath().c_str(), time / (60 * 60),
                        (time / 60) % 60, time % 60, samples_per_pixel,
                        scene.n_FindIntersection);

          float_t const alpha =
              image_contribution_luminance / the_number_of_samples;
          for (int q = 0; q < knImagePixel; ++q) {
            film[q] = filmU[q];
            film[q] += filmP[q] * alpha;
          }

          film.Save(filename, samples_per_pixel);
        }

        // u(x)
        {
          ::_stprintf_s(filename, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                      _TEXT(" %gspp %010lu") _TEXT("u.pfm"),
                        GetOutdirPath().c_str(), time / (60 * 60),
                        (time / 60) % 60, time % 60, samples_per_pixel,
                        scene.n_FindIntersection);

          filmU.Save(filename, samples_per_pixel);
        }

        // p(0)
        {
          ::_stprintf_s(filename, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                      _TEXT(" %gspp %010lu") _TEXT("p.pfm"),
                        GetOutdirPath().c_str(), time / (60 * 60),
                        (time / 60) % 60, time % 60, samples_per_pixel,
                        scene.n_FindIntersection);

          filmP.Save(filename, samples_per_pixel);
        }

        //::_ftprintf_s(stdout, _T("I, %g, %g\n"),
        // image_contribution_luminance_k, alpha);

        if (is_finished || (time >= nRenderingTime)) {
          FinishRendering();
          return;
        }

        begin_time = std::clock() - current_time;  // 時計を戻す
      }
    }
  }
}

void BidirectionalPathTracer_PopulationLightTransportTestbed::run() {
  ::_ftprintf_s(stderr, _TEXT(
                            "** Bidirectional Path Tracing & Metropolis "
                            "Light Transport Testbed **\n"));
  static std::size_t const knInitialSample = 10000;  // 初期サンプル数
  static std::size_t const knPDFSample = 2;  // g(x)/I からのサンプル数

  Scene const &scene = Scene::GetInstance();
  scene.n_FindIntersection = 0;

  int const kImageWidth = scene.GetWidth();
  int const kImageHeight = scene.GetHeight();
  int const knImagePixel = kImageWidth * kImageHeight;

  ImageFilm film(kImageWidth, kImageHeight);   // 結果
  ImageFilm filmU(kImageWidth, kImageHeight);  // u(x)
  ImageFilm filmP(kImageWidth, kImageHeight);  // p(x)

  std::tchar_t filename[MAX_PATH];
  filename[0] = _TEXT('\0');

  // スレッド固有のデータ
  // NOTE: 粒度の小さいものはomp privateとして扱う
  // std::vector<std::mt19937_64> randoms(omp_get_max_threads());
  std::vector<PathVertex> lightsPathVertices(XYZ_CONFIG_kMaxRandomWalkDepth +
                                             2);
  std::vector<PathVertex> theEyePathVertices(XYZ_CONFIG_kMaxRandomWalkDepth +
                                             2);

  // 計測開始
  int const nRenderingTime = GetRenderingTime();
  int const nIntervalTime = GetIntervalTime();
  int saveTime = nIntervalTime;
  std::clock_t begin_time = std::clock();

  // 乱数生成器
  std::vector<MSample> old_primary_samples(knInitialSample);
  std::vector<MSample> new_primary_samples(knInitialSample);

  // 寄与ベクトル
  std::vector<std::vector<pixel_descriptor_t>> old_contribution_vectors(
      knInitialSample);
  std::vector<std::vector<pixel_descriptor_t>> new_contribution_vectors(
      knInitialSample);
  for (std::size_t i = 0; i < knInitialSample; ++i) {
    // 最大経路長を設定
    old_contribution_vectors[i].resize(XYZ_CONFIG_kMaxRandomWalkDepth);
    new_contribution_vectors[i].resize(XYZ_CONFIG_kMaxRandomWalkDepth);
  }

  // 寄与値
  std::vector<float_t> old_contribution_values(knInitialSample);
  std::vector<float_t> new_contribution_values(knInitialSample);

  // 累積密度分布(p_{k}からp_{k+1}のサンプリング用)
  std::vector<float_t> cumulative_mass_function(knInitialSample + 1);

  // Population Annealing を反復する
  float_t image_contribution_luminance = float_t();
  for (std::size_t n = 1;; ++n) {
    //
    // p_{k=0}(x)=u(x) から X_{i} をサンプリング
    //
    cumulative_mass_function[0] = float_t(0);
    for (std::size_t i = 0; i < knInitialSample; ++i) {
      new_primary_samples[i].init(true);

      // 画像平面上の位置を特定
      float_t const x = new_primary_samples[i].next();
      float_t const y = new_primary_samples[i].next();

      // 画素の位置(j)を初期化
      new_contribution_vectors[i][0].first =
          static_cast<std::size_t>(y * kImageHeight) * kImageWidth +
          static_cast<std::size_t>(x * kImageWidth);
      for (std::size_t k = 1; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        new_contribution_vectors[i][k].first = ~0U;
      }

      // シーンをサンプリング
      scene.EvaluateBidirectionalPathTracing(
          x, y, &new_contribution_vectors[i], lightsPathVertices,
          theEyePathVertices, new_primary_samples[i]);

      // 寄与(g(x))を集約
      float_t sampled_value = float_t();
      for (std::size_t k = 0; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        if (~0U != new_contribution_vectors[i][k].first) {
          sampled_value += new_contribution_vectors[i][k].second[1];
        }
      }
      new_contribution_values[i] = sampled_value;

      // 累積密度関数を更新
      cumulative_mass_function[i + 1] =
          cumulative_mass_function[i] + sampled_value;
    }

    if (cumulative_mass_function.back() <= float_t()) {
      continue;  // probably empty scene
    }

    // 正確な I の計算のために寄与を累積
    image_contribution_luminance += cumulative_mass_function.back();

    // 近似的な I を計算
    float_t const image_contribution_luminance_k =
        image_contribution_luminance / (knInitialSample * n);

    // X_{i}~p_{k=0}(x)=u(x) の寄与をバッファに格納
    for (std::size_t i = 0; i < knInitialSample; ++i) {
      // Multiple Importance Sampling 用の重みを計算
      float_t const sampled_value = new_contribution_values[i];
      float_t const heuristic_weight =
          1 /
          (knPDFSample * sampled_value / image_contribution_luminance_k + 1);

      for (std::size_t k = 0; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
        if (~0U != new_contribution_vectors[i][k].first) {
          filmU[new_contribution_vectors[i][k].first] +=
              new_contribution_vectors[i][k].second * heuristic_weight;
        }
      }
    }

    // 再サンプリング
    {
      typedef std::vector<float_t>::const_iterator cmf_const_iterator;
      cmf_const_iterator const begin = cumulative_mass_function.begin();
      cmf_const_iterator const end = cumulative_mass_function.end();
      cmf_const_iterator it = begin;
      for (std::size_t i = 0; i < knInitialSample; ++i) {
        float_t const xi = (i + ::generate_random_number_for_single_thread()) /
                           knInitialSample;
        float_t const value = xi * cumulative_mass_function.back();

        it = std::upper_bound(it, end, value);  // binary search

        assert(it != begin);
        assert(it != end);

        std::size_t const resampled_index = it - begin - 1;
        old_primary_samples[i] = new_primary_samples[resampled_index];
        old_contribution_vectors[i] = new_contribution_vectors[resampled_index];
        old_contribution_values[i] = new_contribution_values[resampled_index];
      }

      // old <=> new を入れ替え
      new_primary_samples.swap(old_primary_samples);  // 過去の値は参照しない
      // new_contribution_vectors = old_contribution_vectors; //
      // 過去の値しか参照しない
      // new_contribution_values.swap(old_contribution_values); //
      // 過去の値は参照しない
    }

    //
    // p_{k=1}(x)=g(x)/I から X_{i} をサンプリング
    //
    for (std::size_t i = 0; i < knInitialSample; ++i) {
      // short length Metropolis-Hastings sampling
      std::vector<pixel_descriptor_t> &old_contribution_vector =
          old_contribution_vectors[i];
      std::vector<pixel_descriptor_t> &new_contribution_vector =
          new_contribution_vectors[i];

      for (std::size_t m = 0;;) {
        new_primary_samples[i].init(false);

        // 画像平面上の位置を特定
        float_t const x = new_primary_samples[i].next();
        float_t const y = new_primary_samples[i].next();

        // 画素の位置(j)を初期化
        new_contribution_vector[0].first =
            static_cast<std::size_t>(y * kImageHeight) * kImageWidth +
            static_cast<std::size_t>(x * kImageWidth);
        for (std::size_t k = 1; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
          new_contribution_vector[k].first = ~0U;
        }

        // シーンをサンプリング
        scene.EvaluateBidirectionalPathTracing(
            x, y, &new_contribution_vector, lightsPathVertices,
            theEyePathVertices, new_primary_samples[i]);

        // 寄与(g(x))を集約
        float_t sampled_value = float_t();
        for (std::size_t k = 0; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
          if (~0U != new_contribution_vector[k].first) {
            sampled_value += new_contribution_vector[k].second[1];
          }
        }
        new_contribution_values[i] = sampled_value;

        // Metropolis-Hastings Ratio から受理確率を計算する
        float_t const a = std::min(float_t(1), new_contribution_values[i] /
                                                   old_contribution_values[i]);

        // Multiple Importance Sampling 用の重みを計算
        // u(x) に対する重みに 1/image_contribution_luminance_k を掛けたもの
        float_t const new_heuristic_weight =
            knPDFSample / (knPDFSample * new_contribution_values[i] +
                           image_contribution_luminance_k);
        float_t const old_heuristic_weight =
            knPDFSample / (knPDFSample * old_contribution_values[i] +
                           image_contribution_luminance_k);

        for (std::size_t k = 0; k < XYZ_CONFIG_kMaxRandomWalkDepth; ++k) {
          if (~0U != new_contribution_vector[k].first) {
            filmP[new_contribution_vector[k].first] +=
                new_contribution_vector[k].second * ((a)*new_heuristic_weight);
          }
          if (~0U != old_contribution_vector[k].first) {
            filmP[old_contribution_vector[k].first] +=
                old_contribution_vector[k].second *
                ((1 - a) * old_heuristic_weight);
          }
        }

        if (++m >= knPDFSample) {
          break;
        }

        if (::generate_random_number_for_single_thread() < a) {
          // 受理
          new_primary_samples[i].accept();
          old_contribution_vector.swap(new_contribution_vector);
          old_contribution_values[i] = new_contribution_values[i];
        }
      }
    }

    //
    // 定期的な保存
    //
    {
      std::clock_t const current_time = std::clock() - begin_time;
      std::clock_t const time = static_cast<std::clock_t>(
          current_time / double(CLOCKS_PER_SEC) + 0.5);
      std::size_t const the_number_of_samples = knInitialSample * n;
      float_t const samples_per_pixel =
          static_cast<float_t>(the_number_of_samples) / knImagePixel;

      // sample per. pixel
      ::_ftprintf_s(stderr, _TEXT("\r  %02u時間%02u分%02u秒 %g[spp]"),
                    time / (60 * 60), (time / 60) % 60, time % 60,
                    samples_per_pixel);

      // save
      bool const is_finished = IsFinished();
      if (is_finished || (time >= saveTime)) {
        saveTime += nIntervalTime;

        {
          ::_stprintf_s(filename, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                      _TEXT(" %gspp %010lu") _TEXT(".pfm"),
                        GetOutdirPath().c_str(), time / (60 * 60),
                        (time / 60) % 60, time % 60, samples_per_pixel,
                        scene.n_FindIntersection);

          float_t const alpha =
              image_contribution_luminance / the_number_of_samples;
          for (int q = 0; q < knImagePixel; ++q) {
            film[q] = filmU[q];
            film[q] += filmP[q] * alpha;
          }

          film.Save(filename, samples_per_pixel);
        }
        // u(x)
        {
          ::_stprintf_s(filename, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                      _TEXT(" %gspp %010lu") _TEXT("u.pfm"),
                        GetOutdirPath().c_str(), time / (60 * 60),
                        (time / 60) % 60, time % 60, samples_per_pixel,
                        scene.n_FindIntersection);

          filmU.Save(filename, samples_per_pixel);
        }

        // p(0)
        {
          ::_stprintf_s(filename, _TEXT("%s") _TEXT("%02dh%02dm%02ds")
                                      _TEXT(" %gspp %010lu") _TEXT("p.pfm"),
                        GetOutdirPath().c_str(), time / (60 * 60),
                        (time / 60) % 60, time % 60, samples_per_pixel,
                        scene.n_FindIntersection);

          filmP.Save(filename, samples_per_pixel);
        }

        //::_ftprintf_s(stdout, _T("I, %g, %g\n"),
        // image_contribution_luminance_k, alpha);

        if (is_finished || (time >= nRenderingTime)) {
          FinishRendering();
          return;
        }

        begin_time = std::clock() - current_time;  // 時計を戻す
      }
    }
  }
}
