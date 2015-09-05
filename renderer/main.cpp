#include "core/Scene.hpp"
#include "core/shell.hpp"
#include "core/Camera.hpp"
#include "core/Film.hpp"
#include <clocale>

#pragma comment(lib, "gdi32")
#pragma comment(lib, "user32")

namespace {
bool CorrectWindowGamma(HDC hDC) {
  if (!hDC) {
    std::tcerr << "can not get dc." << std::endl;
    return false;
  }

  WORD gammaRamp[3][256];

  for (std::size_t i = 0; i < 256; ++i) {
    gammaRamp[0][i] = gammaRamp[1][i] = gammaRamp[2][i] =
        static_cast<WORD>(std::pow(i / 255.0, 1.0) * 65535);
  }

  return !!::SetDeviceGammaRamp(hDC, gammaRamp);
}
}

int _tmain(int argc, std::tchar_t *argv[]) {
#ifndef NDEBUG
  ::_CrtSetDbgFlag(_CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF);
#endif NDEBUG

  // ロケールの設定
  {
    // C++のグローバルロケールの設定
    std::locale::global(std::locale(""));

    // Cのロケールの設定
    ::_tsetlocale(LC_ALL, _T(""));

    // すでに作成されているオブジェクトのロケールを変更
    std::cin.imbue(std::locale(""));
    std::wcin.imbue(std::locale(""));
    std::cout.imbue(std::locale(""));
    std::wcout.imbue(std::locale(""));
    std::cerr.imbue(std::locale(""));
    std::wcerr.imbue(std::locale(""));
  }

  // バージョンの表示
  std::tcerr << RENDERER_NAME << _TEXT(" (build ") << RENDERER_BUILD
             << _TEXT(", ") _TEXT(__TIME__) _TEXT(")") << std::endl;
  std::tcerr << _TEXT("  max threads: ") << omp_get_max_threads() << std::endl;

#if 0
  {
    int const n = omp_get_max_threads();
#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
      HANDLE hThread;
      if (::DuplicateHandle(
        ::GetCurrentProcess(),
        ::GetCurrentThread(),
        ::GetCurrentProcess(),
        &hThread, 0, FALSE,
        DUPLICATE_SAME_ACCESS))
      {
        ::_ftprintf_s(stderr, _TEXT("  CPU[%d] <= Thread[%p]\n"), i, hThread);

        ::SetThreadAffinityMask(hThread, 1 << (i%n));
        ::SetThreadPriority(hThread, THREAD_PRIORITY_HIGHEST);

        ::CloseHandle(hThread);
      }
    }
  }
#endif

  if (argc > 1) {
#if 0
    // パラメータのパース
    std::map<std::tstring, std::tstring> params;

    for (int i = 1; i < argc; ++i)
    {
      if (_TEXT('/') != argv[i][0])
      {
        continue;
      }

      std::tstring param(&argv[i][1]);

      std::size_t pos = param.find(_TEXT('='));
      if (std::tstring::npos == pos)
      {
        continue;
      }

      std::tstring name = param.substr(0, pos);
      hi::trim(name);

      std::tstring value = param.substr(pos + 1);
      hi::trim(value);

      params.insert(std::make_pair(name, value));
    }

    std::map<std::tstring, std::tstring>::iterator it = params.find(_TEXT("file"));
    std::tstring const filename((params.end() != it) ? (*it).second : _TEXT(""));
#else
    std::tstring filename(RENDIR_PATH);
    filename += argv[1];
    filename += _TEXT(".ren");
#endif
    std::tifstream fin(filename.c_str());
    if (!!fin) {
      return tgir::Shell(fin, false);
    }
  }

  std::vector<char *> c_argv;
  for (int i = 0; i < argc; ++i) {
    int size =
        ::WideCharToMultiByte(CP_ACP, 0, argv[i], -1, NULL, 0, NULL, NULL);
    c_argv.push_back(nullptr);
    c_argv.back() = new char[size];
    ::WideCharToMultiByte(CP_ACP, 0, argv[i], -1, c_argv.back(), size, NULL,
                          NULL);
  }

  ::glutInit(&argc, &c_argv[0]);
  ::glutInitWindowPosition(0, 0);
  ::glutInitWindowSize(tgir::Scene::GetInstance().GetWidth(),
                       tgir::Scene::GetInstance().GetHeight());
  ::glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  ::glutCreateWindow(c_argv[0]);

  // ガンマ補正の設定
  if (!::CorrectWindowGamma(::wglGetCurrentDC())) {
    std::tcerr << _TEXT("failed correct display gamma.") << std::endl;
  }

  ::glutDisplayFunc(tgir::Display);
  ::glutReshapeFunc(tgir::Resize);
  ::glutKeyboardFunc(tgir::Keyboard);

  ::glClearColor(0.0, 0.0, 0.0, 1.0);
  ::glLineWidth(2.0);
  ::glEnable(GL_DEPTH_TEST);
#ifdef CONFIG_BACKFACE_CULLING
  ::glEnable(GL_CULL_FACE);
  ::glCullFace(GL_FRONT);  // CullCounterClockwise
#else
  ::glDisable(GL_CULL_FACE);
#endif CONFIG_BACKFACE_CULLING
  ::glEnable(GL_LIGHTING);
  ::glutMainLoop();

  return 0;
}

#if 0
  int w = ::_ttoi(argv[1]);
  int h = ::_ttoi(argv[2]);
  int s = w * h * 3;

  std::vector<float> img1(s);
  std::vector<float> img2(s);

  {
    std::ifstream in(argv[3], std::ios::in | std::ios::binary);
    if (!in)
    {
      std::tcerr << _TEXT("Usage: ") << argv[0] << " width height refarence filenames" << std::endl;
      return 1;
    }
    in.read(reinterpret_cast<char*>(&img1[0]), sizeof(float) * s);
  }

  std::vector<std::tstring> filenames;
  for (int i = 4; i < argc; ++i)
  {
    filenames.push_back(std::tstring(argv[i]));
  }
  std::sort(filenames.Begin(), filenames.End());

  for each (std::tstring const & filename in filenames)
  {
    {
      std::ifstream in(filename.c_str(), std::ios::in | std::ios::binary);
      if (!in)
      {
        continue;
      }
      in.read(reinterpret_cast<char*>(&img2[0]), sizeof(float) * s);
    }

    float absolute_L1 = 0;
    float absolute_L2 = 0;
    float absolute_Lmax = 0;
    float relative_L1 = 0;
    float relative_L2 = 0;
    float relative_Lmax = 0;
    float tmp[3];
    for (int i = 0; i < s; i += 3)
    {
      tgir::CCIR601_1_RGB2XYZ(img1.Begin()+i, tmp); float a = tmp[1];
      tgir::CCIR601_1_RGB2XYZ(img2.Begin()+i, tmp); float b = tmp[1];

      float const absolute_error = std::abs(b - a);
      absolute_L1 += absolute_error;
      absolute_L2 += hi::square_of(absolute_error);
      hi::max(absolute_Lmax, absolute_error);

      if (a > hi::rcp<float>(256))
      {
        tgir::Real const relative_error = absolute_error / a;
        relative_L1 += relative_error;
        relative_L2 += hi::square_of(relative_error);
        hi::max(relative_Lmax, relative_error);
      }
    }
    absolute_L2 = std::sqrt(absolute_L2);
    relative_L2 = std::sqrt(relative_L2);

    std::tcerr
      << absolute_L1 << _TEXT(",") << absolute_L2 << _TEXT(",") << absolute_Lmax << _TEXT(",")
      << relative_L1 << _TEXT(",") << relative_L2 << _TEXT(",") << relative_Lmax << std::endl;
  }
#endif
/**
  {
    SYSTEM_INFO system_info;
    ::GetSystemInfo(&system_info);

    std::tcerr << _TEXT("  processor architecture: ");
    switch (system_info.wProcessorArchitecture)
    {
    case PROCESSOR_ARCHITECTURE_INTEL:
      std::tcerr
        << _TEXT("x86 Family ") << system_info.wProcessorLevel
        << _TEXT(" Model ")     << ((system_info.wProcessorRevision >> 8) &
  0xFF) // 80386と80486では厳密には異なる
        << _TEXT(" Stepping ")  << (system_info.wProcessorRevision & 0xFF); //
  同上
      break;
    case PROCESSOR_ARCHITECTURE_MIPS:
      if (4 == system_info.wProcessorLevel)
      {
        std::tcerr << _TEXT("MIPS R4000");
      }
      else
      {
        std::tcerr << _TEXT("MIPS ProcessorLebel ") <<
  system_info.wProcessorLevel;
      }
      break;
    case PROCESSOR_ARCHITECTURE_ALPHA:
      switch (system_info.wProcessorLevel)
      {
      case 21064:
        std::tcerr << _TEXT("Alpha 21064");
        break;
      case 21066:
        std::tcerr << _TEXT("Alpha 21066");
        break;
      case 21164:
        std::tcerr << _TEXT("Alpha 21164");
        break;
      default:
        std::tcerr << _TEXT("Alpha ProcessorLebel ") <<
  system_info.wProcessorLevel;
        break;
      }
      break;
    case PROCESSOR_ARCHITECTURE_PPC:
      switch (system_info.wProcessorLevel)
      {
      case 1:
        std::tcerr << _TEXT("PPC 601");
        break;
      case 3:
        std::tcerr << _TEXT("PPC 603");
        break;
      case 4:
        std::tcerr << _TEXT("PPC 604");
        break;
      case 6:
        std::tcerr << _TEXT("PPC 603+");
        break;
      case 9:
        std::tcerr << _TEXT("PPC 604+");
        break;
      case 20:
        std::tcerr << _TEXT("PPC 620");
        break;
      default:
        std::tcerr << _TEXT("PPC ProcessorLebel ") <<
  system_info.wProcessorLevel;
        break;
      }
      break;
    }
    std::tcerr << std::endl;

    std::tcerr << _TEXT("  processor type: ");
    switch (system_info.dwProcessorType)
    {
    case PROCESSOR_INTEL_486:
      std::tcerr << _TEXT("Intel 486 Processor");
      break;
    case PROCESSOR_INTEL_PENTIUM:
      std::tcerr << _TEXT("Intel Pentium Processor");
      break;
    case PROCESSOR_ARCHITECTURE_UNKNOWN:
      std::tcerr << _TEXT("Unknown Processor");
      break;
    }
    std::tcerr << std::endl;
    std::tcerr << _TEXT("  the number of processors: ")
      << system_info.dwNumberOfProcessors << std::endl;
  }
*/