#include "core/scene.hpp"
#include "core/shell.hpp"
#include "core/thinlenscamera.hpp"
#include "core/imagefilm.hpp"
#include <clocale>

#pragma comment(lib, "gdi32")
#pragma comment(lib, "user32")

namespace {
bool CorrectWindowGamma(HDC hDC) {
  if (!hDC) {
    ::_ftprintf_s(stderr, _TEXT("can not get dc.\n"));
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
  ::_ftprintf_s(stderr, _TEXT("%s (build %s, %s)\n"), RENDERER_NAME,
                RENDERER_BUILD, _TEXT(__TIME__));
  ::_ftprintf_s(stderr, _TEXT("max threads: %d\n"), omp_get_max_threads());

  if (argc > 1) {
    std::tstring filename(RENDIR_PATH);
    filename += argv[1];
    filename += _TEXT(".ren");

    std::tifstream fin(filename.c_str());
    if (!!fin) {
      return xyz::Shell(fin, false);
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
  ::glutInitWindowSize(xyz::Scene::GetInstance().GetWidth(),
                       xyz::Scene::GetInstance().GetHeight());
  ::glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  ::glutCreateWindow(c_argv[0]);

  // ガンマ補正の設定
  if (!::CorrectWindowGamma(::wglGetCurrentDC())) {
    ::_ftprintf_s(stderr, _TEXT("failed correct display gamma.\n"));
  }

  ::glutDisplayFunc(xyz::Display);
  ::glutReshapeFunc(xyz::Resize);
  ::glutKeyboardFunc(xyz::Keyboard);

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
