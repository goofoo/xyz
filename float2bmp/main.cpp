// main.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

#include "stdafx.hpp"

#pragma comment(lib, "Rpcrt4")  // for UuidCreate(), UuidToString()

namespace {
bool load_pfm(std::tstring const &filename, std::vector<float> &floatmap,
              int *const pWidth, int *const pHeight) {
  std::ifstream fin(filename.c_str(), std::ios::in | std::ios::binary);

  if (!fin) {
    return false;
  }

  std::string str;
  std::getline(fin, str);
  if ("PF" != str) {
    return false;
  }

  std::getline(fin, str);
  {
    std::istringstream in(str);

    in >> *pWidth >> *pHeight;

    if ((*pWidth <= 0) || (pHeight <= 0)) {
      return false;
    }
  }

  std::getline(fin, str);
  {
    std::istringstream in(str);

    int a;
    in >> a;

    if (-1 != a) {
      return false;
    }
  }

  floatmap.resize(*pWidth * *pHeight * 3);
  fin.read(reinterpret_cast<char *>(&floatmap[0]),
           static_cast<std::streamsize>(floatmap.size()) * sizeof(float));

  return true;
}

bool save_bitmap(std::tstring const &filename, BITMAPFILEHEADER const &bmh,
                 BITMAPINFOHEADER const &bmi,
                 std::vector<hi::ubyte> const &bitmap) {
  std::ofstream fout(filename.c_str(), std::ios::out | std::ios::binary);

  if (!fout) {
    return false;
  }

  fout.write(reinterpret_cast<char const *>(&bmh), sizeof(bmh));
  fout.write(reinterpret_cast<char const *>(&bmi), sizeof(bmi));
  fout.write(reinterpret_cast<char const *>(&bitmap[0]),
             static_cast<std::streamsize>(bitmap.size()) * sizeof(hi::ubyte));

  return true;
}

inline hi::ubyte clamp(float const value) {
  if (value <= 0) {
    return 0;
  }
  if (value >= 1) {
    return 255;
  }
  return static_cast<hi::ubyte>(value * 256);  // 量子化
}

inline hi::ubyte tonemap(float const light,
                         float const white,  // 白色にマッピングしたい明るさ
                         float const exposure,  // 露出
                         float const dst_gamma = 2.2f,
                         float const src_gamma = 1.0f) {
  if (light <= 0) {
    return 0;
  }

// トーンマッピング
#if 0
    // なし
    UNREFERENCED_PARAMETER(white);
    float const toonmapped_value = exposure * light;
#elif 1
  // フィルム モデル
  UNREFERENCED_PARAMETER(white);
  float const toonmapped_value = 1.0f - std::exp(-exposure * light);
#elif 0
  // Reinhard モデル(white=∞)
  UNREFERENCED_PARAMETER(white);
  float const toonmapped_value = exposure * light / (1 + light);
#else
  // Reinhard モデル
  float const toonmapped_value =
      exposure * light * (1 + light / hi::square_of(white)) / (1 + light);
#endif

// ガンマ補正
#if 0
    // なし
    float const gammacorreted_value = toonmapped_value;
#else
  // あり
  float const gamma = src_gamma / dst_gamma;
  float const gammacorreted_value = std::pow(toonmapped_value, gamma);
#endif

  return clamp(gammacorreted_value);
}
}

int _tmain(int argc, std::tchar_t *argv[]) {
  if (argc < 2) {
    std::tcerr << _TEXT("usage: ") << argv[0]
               << _TEXT(" width height dirname [reference]") << std::endl;
    std::getchar();
    return 1;
  }

  float const exposure =
      (argc < 3) ? 1.0f : static_cast<float>(::_tstof(argv[2]));

  // 保存ファイル名
  std::tchar_t filename[MAX_PATH];

  // 保存ビットマップのヘッダ
  BITMAPFILEHEADER bmh;
  BITMAPINFOHEADER bmi;

  ::memset(&bmh, 0, sizeof(bmh));
  bmh.bfType = 'B' | ('M' << 8);
  bmh.bfOffBits = sizeof(bmh) + sizeof(bmi);

  ::memset(&bmi, 0, sizeof(bmi));
  bmi.biSize = sizeof(bmi);
  bmi.biPlanes = 1;
  bmi.biBitCount = 24;
  bmi.biCompression = BI_RGB;

  // 保存ビットマップのデータ
  std::vector<hi::ubyte> bitmap;

  // 対象ディレクトリ
  std::tstring dirname(argv[1]);
  std::replace(dirname.begin(), dirname.end(), _TEXT('/'), _TEXT('\\'));
  if (_TEXT('\\') != dirname[dirname.size() - 1]) {
    dirname += _TEXT("\\");
  }
  ::_ftprintf_s(stdout, _T("%s\n"), dirname.c_str());

  // ファイルの列挙
  std::vector<std::tstring> files;
  hi::enum_files(dirname, _TEXT(".pfm"), files);

  std::vector<float> floatmap;

  {
    // 順に処理する
    for (std::size_t i = 0, size = files.size(); i < size; ++i) {
      // 読み込む
      ::_ftprintf_s(stdout, _T("%s\n"), files[i].c_str());
      files[i].insert(0, dirname);

      int width;
      int height;
      if (!::load_pfm(files[i].c_str(), floatmap, &width, &height)) {
        continue;
      }
      if (width & 3) {
        std::tcerr << _TEXT("画像の横幅は四の倍数でなければなりません．")
                   << std::endl;
        continue;
      }
      int const count = width * height * 3;

      float max_value = 0;
#if 1
      for (int p = 0; p < count; p += 3) {
        float const luminance = 0.6069f * floatmap[p + 0] +
                                0.1735f * floatmap[p + 1] +
                                0.2003f * floatmap[p + 2];
        if (max_value < luminance) {
          max_value = luminance;
        }
      }
#else
      for (int y = height / 5, ymax = height * 4 / 5; y < ymax; ++y) {
        for (int x = width / 5, xmax = width * 4 / 5; x < xmax; ++x) {
          int const p = x + y * width;
          float const luminance = 0.6069f * floatmap[p + 0] +
                                  0.1735f * floatmap[p + 1] +
                                  0.2003f * floatmap[p + 2];
          if (max_value < luminance) {
            max_value = luminance;
          }
        }
      }
#endif

      bitmap.resize(count);
      for (int p = 0; p < count; p += 3) {
        bitmap[p + 2] = ::tonemap(floatmap[p + 0], max_value, exposure);
        bitmap[p + 1] = ::tonemap(floatmap[p + 1], max_value, exposure);
        bitmap[p + 0] = ::tonemap(floatmap[p + 2], max_value, exposure);
      }

      // 連番で保存
      bmi.biWidth = width;
      bmi.biHeight = height;
      bmh.bfSize = bmh.bfOffBits + count * sizeof(hi::ubyte);
      ::_stprintf_s(filename, MAX_PATH, _TEXT("%s%06u.bmp"), dirname.c_str(),
                    i + 1);
      ::save_bitmap(filename, bmh, bmi, bitmap);
    }
  }

  return 0;
}
#if 0
namespace
{
  bool load_pfm(
    std::tstring const & filename,
    std::vector<float> & floatmap,
    int const width, int const height)
  {
    std::ifstream fin(filename.c_str(), std::ios::in | std::ios::binary);

    if (!fin) {
      return false;
    }

    std::string str;
    std::getline(fin, str);
    if ("PF" != str)
    {
      return false;
    }

    std::getline(fin, str);
    {
      std::istringstream in(str);

      int w, h;
      in >> w >> h;

      if ((width != w) || (height != h))
      {
        return false;
      }
    }

    std::getline(fin, str);
    {
      std::istringstream in(str);

      int a;
      in >> a;

      if (-1 != a)
      {
        return false;
      }
    }

    fin.read(reinterpret_cast<char*>(&floatmap[0]),
      static_cast<std::streamsize>(floatmap.size()) * sizeof(float));

    return true;
  }

  bool load_float(
    std::tstring const & filename,
    std::vector<float> & floatmap)
  {
    std::ifstream fin(filename.c_str(), std::ios::in | std::ios::binary);

    if (!fin)
    {
      return false;
    }

    fin.read(reinterpret_cast<char*>(&floatmap[0]),
      static_cast<std::streamsize>(floatmap.size()) * sizeof(float));

    return true;
  }

  bool save_bitmap(
    std::tstring const & filename,
    BITMAPFILEHEADER const & bmh,
    BITMAPINFOHEADER const & bmi,
    std::vector<hi::ubyte> const & bitmap)
  {
    std::ofstream fout(filename.c_str(), std::ios::out | std::ios::binary);

    if (!fout)
    {
      return false;
    }

    fout.write(reinterpret_cast<char const *>(&bmh), sizeof(bmh));
    fout.write(reinterpret_cast<char const *>(&bmi), sizeof(bmi));
    fout.write(reinterpret_cast<char const *>(&bitmap[0]),
      static_cast<std::streamsize>(bitmap.size()) * sizeof(hi::ubyte));

    return true;
  }

  void hsv2rgb(float const hsv[], float rgb[])
  {
    float h = (hsv[0] - std::floor(hsv[0])) * 6;
    float s = hsv[1];
    float v = hsv[2];

    int   i = static_cast<int>(h);
    float f = h - i;

    float p = v * (1 - s        );
    float q = v * (1 - s * ( f));
    float t = v * (1 - s * (1-f));

    switch (i)
    {
    case 6:
    case 0: rgb[0] = v; rgb[1] = t; rgb[2] = p; break;
    case 1: rgb[0] = q; rgb[1] = v; rgb[2] = p; break;
    case 2: rgb[0] = p; rgb[1] = v; rgb[2] = t; break;
    case 3: rgb[0] = p; rgb[1] = q; rgb[2] = v; break;
    case 4: rgb[0] = t; rgb[1] = p; rgb[2] = v; break;
    case 5: rgb[0] = v; rgb[1] = p; rgb[2] = q; break;
    }
  }

  inline hi::ubyte clamp(float const value)
  {
    if (value <= 0)
    {
      return 0;
    }
    if (value >= 1)
    {
      return 255;
    }

    return static_cast<hi::ubyte>(value * 256); // 量子化
  }

  inline hi::ubyte tonemap(
    float const light,
    float const white, // 白色にマッピングしたい明るさ
    float const exposure = 3.0f, // 露出
    float const dst_gamma = 2.2f,
    float const src_gamma = 1.0f)
  {
    if (light <= 0)
    {
      return 0;
    }
 
    // トーンマッピング
#if 0
    // なし
    float const toonmapped_value = exposure * light;
#elif 1
    // フィルム モデル
    float const toonmapped_value = 1 - std::exp(-exposure * light);
#elif 1
    // Reinhard モデル(white=∞)
    float const toonmapped_value = exposure * light / (1 + light);
#else
    // Reinhard モデル
    float const toonmapped_value = exposure * light * (1 + light / hi::square_of(white)) / (1 + light);
#endif

    // ガンマ補正
#if 0
    // なし
    float const gammacorreted_value = toonmapped_value;
#else
    // あり
    float const gamma = src_gamma / dst_gamma;
    float const gammacorreted_value = std::pow(toonmapped_value, gamma);
#endif

    return clamp(gammacorreted_value);
  }
}

int _tmain(int argc, std::tchar_t * argv[])
{
  if (argc < 4)
  {
    std::tcerr << _TEXT("usage: ") << argv[0]
      << _TEXT(" width height dirname [reference]") << std::endl;
    std::getchar();
    return 1;
  }

  int const width  = ::_ttoi(argv[1]);
  int const height = ::_ttoi(argv[2]);
  int const count  = width * height * 3;

  if (width & 3)
  {
    std::tcerr << _TEXT("画像の横幅は四の倍数でなければなりません．") << std::endl;
    return 1;
  }

  // 保存ファイル名
  std::tchar_t filename[MAX_PATH];

  // 保存ビットマップのヘッダ
  BITMAPFILEHEADER bmh;
  BITMAPINFOHEADER bmi;

  ::memset(&bmh, 0, sizeof(bmh));
  bmh.bfType = 'B' | ('M' << 8);
  bmh.bfOffBits = sizeof(bmh) + sizeof(bmi);
  bmh.bfSize = bmh.bfOffBits + count * sizeof(hi::ubyte);

  ::memset(&bmi, 0, sizeof(bmi));
  bmi.biSize = sizeof(bmi);
  bmi.biWidth = width;
  bmi.biHeight = height;
  bmi.biPlanes = 1;
  bmi.biBitCount = 24;
  bmi.biCompression = BI_RGB;

  // 保存ビットマップのデータ
  std::vector<hi::ubyte> bitmap(count);

  // 対象ディレクトリ
  std::tstring dirname(argv[3]);
  std::replace(dirname.begin(), dirname.end(), _TEXT('/'), _TEXT('\\'));
  if (_TEXT('\\') != dirname[dirname.size()-1])
  {
    dirname += _TEXT("\\");
  }
  ::_ftprintf_s(stdout, _T("%s\n"), dirname.c_str());

  // ファイルの列挙
  std::vector<std::tstring> files;
  hi::enum_files(dirname, _TEXT(".pfm"), files);

  std::vector<float> testimage(count);

  if (argc > 4)
  {
    // リファレンス
    std::vector<float> reference(count);
    if (!load_pfm(argv[4], reference, width, height))
    {
      ::_ftprintf_s(stdout, _T("リファレンス(%s)のロードに失敗しました．\n"), argv[4]);
      return 1;
    }
#if 0
    // 信号の実行値を求める
    long double signal = 0;
    for (int p = 0; p < count; p += 3)
    {
      hi::basic_vector3<double> t;
      hi::CCIR601_1_RGB2XYZ(reference.begin() + p, t);

      long double const error = hi::length_squared(t);
      signal += error;
    }
    signal /= count * 3;
    signal  = std::sqrt(signal);

    ::_ftprintf_s(stdout, _TEXT("signal = %g"), signal);
#endif
#if 1
    // 実際の計算
    std::tstring name;
    {
      UUID uuid;
      while (::UuidCreate(&uuid) != RPC_S_OK); // UUIDの生成

      RPC_WSTR guid = NULL;
      while (::UuidToString(&uuid, &guid) != RPC_S_OK); // 文字列化(GUID化)
    
      name = reinterpret_cast<wchar_t>(guid);

      ::RpcStringFree(&guid);
    }

    std::tstring logfile(dirname + name + _T(".csv"));
    std::ofstream logger(logfile.c_str());
    if (!logger)
    {
      ::_ftprintf_s(stdout, _T("ログファイル(%s)のロードに失敗しました．\n"), logfile.c_str());
      return 1;
    }

    logger << "RMSE" << std::endl;
    // 順に処理する
    for (std::size_t i = 0, size = files.size(); i < size; ++i)
    {
      // 読み込む
      ::_ftprintf_s(stdout, _T("%s\n"), files[i].c_str());
      files[i].insert(0, dirname);

      if (!::load_pfm(files[i].c_str(), testimage, width, height))
      {
        continue;
      }

      long double rmse = 0;
      for (int p = 0; p < count; p += 3)
      {
        hi::basic_vector3<double> s;
        hi::CCIR601_1_RGB2XYZ(testimage.begin() + p, s);

        hi::basic_vector3<double> t;
        hi::CCIR601_1_RGB2XYZ(reference.begin() + p, t);

        hi::basic_vector3<double> e(s - t);

        // XYZ色空間で考える
        long double const error = hi::length_squared(e);
        rmse += error;

        // dB表現に直す
        long double const error_level = 10 * std::log10(error / 3);

        // min_level から max_level までの表現に直す
        long double const min_level = -50;
        long double const max_level =  1;
        long double const nomalized_level =
          std::min(std::max(
            (error_level - min_level) /
            (  max_level - min_level), 0.0L), 1.0L);

        // HSV色空間で表現する
        float hsv[3];
        hsv[0] = static_cast<float>((1 - nomalized_level) * 2 / 3);
        hsv[1] = 1;
        hsv[2] = 1;

        // RGB色空間に変換する
        float rgb[3];
        hsv2rgb(hsv, rgb);
        bitmap[p+2] = ::tonemap(rgb[0], 1);
        bitmap[p+1] = ::tonemap(rgb[1], 1);
        bitmap[p+0] = ::tonemap(rgb[2], 1);
      }
      rmse /= count * 3;
      rmse  = std::sqrt(rmse);

      // 誤差のログを出力
      logger << rmse << "," << 20 * std::log10(rmse) << std::endl;

      // 連番で保存
      ::_stprintf_s(filename, MAX_PATH,
        _TEXT("%s%06u RMSE=%g.bmp"),
        dirname.c_str(), i + 1, rmse);
      ::save_bitmap(filename, bmh, bmi, bitmap);
    }
#endif
  }
  else
  {
    // 順に処理する
    for (std::size_t i = 0, size = files.size(); i < size; ++i)
    {
      // 読み込む
      ::_ftprintf_s(stdout, _T("%s\n"), files[i].c_str());
      files[i].insert(0, dirname);

      if (!::load_pfm(files[i].c_str(), testimage, width, height))
      {
        continue;
      }

      float max_value = 0;
      for (int p = 0; p < count; p += 3)
      {
        float const luminance =
          0.6069f * testimage[p+0] +
          0.1735f * testimage[p+1] +
          0.2003f * testimage[p+2];
        if (max_value < luminance) {
          max_value = luminance;
        }
      }

      for (int p = 0; p < count; p += 3) {
        bitmap[p+2] = ::tonemap(testimage[p+0], max_value);
        bitmap[p+1] = ::tonemap(testimage[p+1], max_value);
        bitmap[p+0] = ::tonemap(testimage[p+2], max_value);
      }

      // 連番で保存
      ::_stprintf_s(filename, MAX_PATH, _TEXT("%s%06u.bmp"), dirname.c_str(), i + 1);
      ::save_bitmap(filename, bmh, bmi, bitmap);
    }
  }

  return 0;
}

/*
    logger << "絶対誤差,,,相対誤差,," << std::endl;
    logger << "L1ノルム,L2ノルム,最大値ノルム,L1ノルム,L2ノルム,最大値ノルム" << std::endl;

    // 順に処理する
    for (std::size_t i = 0, size = files.size(); i < size; ++i)
    {
      // 読み込む
      ::_ftprintf_s(stdout, _T("%s\n"), files[i].c_str());
      files[i].insert(0, dirname);

      if (!::load_pfm(files[i].c_str(), testimage, width, height))
      {
        continue;
      }

      // 絶対誤差
      long double absolute_L1 = 0;
      long double absolute_L2 = 0;
      long double absolute_Lmax = 0;

      for (int p = 0; p < count; p += 3)
      {
        hi::basic_vector3<double> s;
        hi::CCIR601_1_RGB2XYZ(testimage.begin() + p, s);

        hi::basic_vector3<double> t;
        hi::CCIR601_1_RGB2XYZ(reference.begin() + p, t);

        // XYZ色空間で考える
        long double const absolute_error = hi::absolute_error(s, t);

        // 絶対誤差を蓄積する
        absolute_L1 += absolute_error;
        absolute_L2 += hi::square_of(absolute_error);
        hi::max(absolute_Lmax, absolute_error);

        // HSV色空間で表現する
        float hsv[3];
        hsv[0] = static_cast<float>((absolute_error > 1) ? 0 : (1 - absolute_error) * 2 / 3);
        hsv[1] = 1;
        hsv[2] = 1;

        // RGB色空間に変換する
        float rgb[3];
        hsv2rgb(hsv, rgb);
        bitmap[p+2] = ::tonemap(rgb[0], 1);
        bitmap[p+1] = ::tonemap(rgb[1], 1);
        bitmap[p+0] = ::tonemap(rgb[2], 1);
      }
      absolute_L2 = std::sqrt(absolute_L2);

      // 連番で保存
      ::_stprintf_s(filename, MAX_PATH, _TEXT("%sABSOLUTE %06u L1=%g L2=%g Lmax=%g.bmp"),
        dirname.c_str(), i + 1, absolute_L1, absolute_L2, absolute_Lmax);
      ::save_bitmap(filename, bmh, bmi, bitmap);

      // 相対誤差
      long double relative_L1 = 0;
      long double relative_L2 = 0;
      long double relative_Lmax = 0;

      for (int p = 0; p < count; p += 3)
      {
        hi::basic_vector3<double> s;
        hi::CCIR601_1_RGB2XYZ(testimage.begin() + p, s);

        hi::basic_vector3<double> t;
        hi::CCIR601_1_RGB2XYZ(reference.begin() + p, t );

        // XYZ色空間で考える
        long double const relative_error = hi::relative_error(s, t);

        // 相対誤差を蓄積する
        relative_L1 += relative_error;
        relative_L2 += hi::square_of(relative_error);
        hi::max(relative_Lmax, relative_error);

        // HSV色空間で表現する
        float const hValue = static_cast<float>(relative_error * 0.5);
        float hsv[3];
        hsv[0] = (hValue > 1) ? 0 : (1 - hValue) * 2 / 3;
        hsv[1] = 1;
        hsv[2] = 1;

        // RGB色空間に変換する
        float rgb[3];
        hsv2rgb(hsv, rgb);
        bitmap[p+2] = ::clamp(rgb[0]);
        bitmap[p+1] = ::clamp(rgb[1]);
        bitmap[p+0] = ::clamp(rgb[2]);
      }
      relative_L2 = std::sqrt(relative_L2);

      // 連番で保存
      ::_stprintf_s(filename, MAX_PATH, _TEXT("%sRELATIVE %06u L1=%g L2=%g Lmax=%g.bmp"),
        dirname.c_str(), i + 1, relative_L1, relative_L2, relative_Lmax);
      ::save_bitmap(filename, bmh, bmi, bitmap);
    
      // 誤差のログを出力
      logger
        << absolute_L1 << "," << absolute_L2 << "," << absolute_Lmax << ","
        << relative_L1 << "," << relative_L2 << "," << relative_Lmax << std::endl;
    }
*/
#endif
