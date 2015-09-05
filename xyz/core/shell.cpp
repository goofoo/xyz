#define PRAGMATIZE_GLOBAL_OBJECT
#include "shell.hpp"
#include "scene.hpp"
#include "../integrator/integrator.hpp"

namespace {
// Scene Definition
std::tstring s_scene_name;
std::tstring s_outdir_path;

// System Setting Parameters
int s_rendering_time;
int s_interval_time;
int s_algorithm_type;
bool volatile s_is_rendering;

std::auto_ptr<hi::thread> s_kicker;

int stream2time(std::tistream &in) {
  int time = 0;

  while (!!in) {
    int val;
    in >> val;
    if (!in) {
      break;
    }

    std::tstring unit;
    in >> unit;
    if ((_TEXT("d") == unit) || (_TEXT("day") == unit) ||
        (_TEXT("days") == unit)) {
      time += val * (60 * 60 * 24);
    } else if ((_TEXT("h") == unit) || (_TEXT("hr") == unit) ||
               (_TEXT("hour") == unit) || (_TEXT("hours") == unit)) {
      time += val * (60 * 60);
    } else if ((_TEXT("m") == unit) || (_TEXT("min") == unit) ||
               (_TEXT("minute") == unit) || (_TEXT("minutes") == unit)) {
      time += val * 60;
    } else if ((_TEXT("s") == unit) || (_TEXT("sec") == unit) ||
               (_TEXT("second") == unit) || (_TEXT("seconds") == unit)) {
      time += val;
    }
  }

  return time;
}

std::tstring time2string(int time) {
  int days = time / (60 * 60 * 24);
  int hours = (time / (60 * 60)) % 24;
  int minutes = (time / 60) % 60;
  int seconds = time % 60;

  std::tostringstream out;

  if (days > 0) {
    out << days << ((days > 1) ? _TEXT(" days") : _TEXT(" day"));
  }

  if (hours > 0) {
    out << ((days > 0) ? _TEXT(" ") : _TEXT("")) << hours
        << ((hours > 1) ? _TEXT(" hours") : _TEXT(" hour"));
  }

  if (minutes > 0) {
    out << (((days > 0) || (hours > 0)) ? _TEXT(" ") : _TEXT("")) << minutes
        << ((minutes > 1) ? _TEXT(" minutes") : _TEXT(" minute"));
  }

  if (seconds > 0) {
    out << (((days > 0) || (hours > 0) || (minutes > 0)) ? _TEXT(" ")
                                                         : _TEXT("")) << seconds
        << ((seconds > 1) ? _TEXT(" seconds") : _TEXT(" second"));
  }

  if ((0 == days) && (0 == hours) && (0 == minutes) && (0 == seconds)) {
    out << _TEXT(" 0 s");
  }

  return out.str();
}
}

namespace xyz {
std::tstring const &GetOutdirPath() { return s_outdir_path; }

int GetRenderingTime() { return s_rendering_time; }

int GetIntervalTime() { return s_interval_time; }

bool IsFinished() { return !s_is_rendering; }

void FinishRendering() { s_is_rendering = false; }
}

namespace xyz {
void Display() {
  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  { Scene::GetInstance().Render(); }
  ::glutSwapBuffers();
}

void Resize(int w, int h) {
  Scene const &scene = Scene::GetInstance();
  if ((w != scene.GetWidth()) || (h != scene.GetHeight())) {
    ::glutReshapeWindow(scene.GetWidth(), scene.GetHeight());
  } else {
    ::glViewport(0, 0, w, h);
  }
}

void Keyboard(unsigned char key, int, int) {
  switch (key) {
  // カメラ位置の移動
  case 'a':
    if (!s_is_rendering)
      Scene::GetInstance().ThinLensCamera().MoveBy(float3_t(-0.05, 0, 0));
    break;
  case 'd':
    if (!s_is_rendering)
      Scene::GetInstance().ThinLensCamera().MoveBy(float3_t(+0.05, 0, 0));
    break;
  case 's':
    if (!s_is_rendering)
      Scene::GetInstance().ThinLensCamera().MoveBy(float3_t(0, -0.05, 0));
    break;
  case 'w':
    if (!s_is_rendering)
      Scene::GetInstance().ThinLensCamera().MoveBy(float3_t(0, +0.05, 0));
    break;
  case 'e':
    if (!s_is_rendering)
      Scene::GetInstance().ThinLensCamera().MoveBy(float3_t(0, 0, -0.05));
    break;
  case 'q':
    if (!s_is_rendering)
      Scene::GetInstance().ThinLensCamera().MoveBy(float3_t(0, 0, +0.05));
    break;

  // 注視点の移動
  case 'j':
    if (!s_is_rendering)
      Scene::GetInstance().ThinLensCamera().LookBy(float3_t(-0.05, 0, 0));
    break;
  case 'l':
    if (!s_is_rendering)
      Scene::GetInstance().ThinLensCamera().LookBy(float3_t(+0.05, 0, 0));
    break;
  case 'k':
    if (!s_is_rendering)
      Scene::GetInstance().ThinLensCamera().LookBy(float3_t(0, -0.05, 0));
    break;
  case 'i':
    if (!s_is_rendering)
      Scene::GetInstance().ThinLensCamera().LookBy(float3_t(0, +0.05, 0));
    break;
  case 'o':
    if (!s_is_rendering)
      Scene::GetInstance().ThinLensCamera().LookBy(float3_t(0, 0, -0.05));
    break;
  case 'u':
    if (!s_is_rendering)
      Scene::GetInstance().ThinLensCamera().LookBy(float3_t(0, 0, +0.05));
    break;

  // 表示方法の切り替え
  case '1':
    ::glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
    break;
  case '2':
    ::glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    break;
  case '3':
    ::glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    break;

  // コマンドシェル
  case 'x':
    Shell(std::tcin);
    break;

  // アプリケーションの終了
  case '\033':
    std::exit(0);
    return;

  default:
    break;
  }

  ::glutPostRedisplay();
}
}

namespace xyz {
enum interpret_result {
  SHELL_EXIT,
  SHELL_QUIT,
  SHELL_CONTINUE,
  SHELL_UNKOWN,
};

interpret_result Interpret(std::tstring &command, bool const wait_shell) {
  std::tistringstream sin(command);
  sin >> command;

  if (_TEXT("quit") == command) {
    return SHELL_QUIT;
  }

  if (_TEXT("exit") == command) {
    return SHELL_EXIT;
  }

  if (_TEXT("abort") == command) {
    FinishRendering();
    return SHELL_QUIT;
  }

  if (!IsFinished()) {
    return SHELL_CONTINUE;
  }

  if (_TEXT("set") == command) {
    sin >> command;
    if (_TEXT("rep") == command) {
      sin >> s_algorithm_type;
      std::tcerr << _TEXT("rep (") << s_algorithm_type << _TEXT(")");
    } else if (_TEXT("rt") == command) {
      s_rendering_time = ::stream2time(sin);
      if (s_interval_time > s_rendering_time) {
        s_interval_time = s_rendering_time;
      }
      std::tcerr << _TEXT("it (") << s_interval_time << _TEXT("), ")
                 << _TEXT("rt (") << s_rendering_time << _TEXT(")");
    } else if (_TEXT("it") == command) {
      s_interval_time = ::stream2time(sin);
      if (s_interval_time > s_rendering_time) {
        s_rendering_time = s_interval_time;
      }
      std::tcerr << _TEXT("it (") << s_interval_time << _TEXT("), ")
                 << _TEXT("rt (") << s_rendering_time << _TEXT(")");
    } else if (_TEXT("eye") == command) {
      float3_t eye;
      sin >> eye;
      Scene::GetInstance().ThinLensCamera().MoveTo(eye);

      std::tcerr << _TEXT("eye (") << eye[0] << _TEXT(", ") << eye[1]
                 << _TEXT(", ") << eye[2] << _TEXT(")");
    } else if (_TEXT("at") == command) {
      float3_t at;
      sin >> at;
      Scene::GetInstance().ThinLensCamera().LookAt(at);

      std::tcerr << _TEXT("at (") << at[0] << _TEXT(", ") << at[1]
                 << _TEXT(", ") << at[2] << _TEXT(")");
    } else if (_TEXT("dim") == command) {
      int width;
      int height;

      sin >> width >> height;

      if (width <= 0) {
        width = 1;
      }

      if (height <= 0) {
        height = 1;
      }

      Scene::GetInstance().SetDimension(width, height);

      std::tcerr << _TEXT("dim (") << width << _TEXT(", ") << height
                 << _TEXT(")");
    }

    return SHELL_CONTINUE;
  }

  if (_TEXT("load") == command) {
    sin >> s_scene_name;

    std::tstring filename(MQODIR_PATH);
    filename += s_scene_name;
    filename += _TEXT(".mqo");

    // load scene
    std::tifstream in(filename.c_str());
    if (!!in) {
      std::tcerr << std::endl
                 << _TEXT("load scene (") << filename << _TEXT(") ... ")
                 << std::endl;
      {
        std::clock_t begin_times = std::clock();
        Scene::GetInstance().Load(in);  // {sponza, winosi, scene1, scene2}

        std::clock_t times = std::clock() - begin_times;
        times = static_cast<std::clock_t>(times / double(CLOCKS_PER_SEC) + 0.5);
        ::_ftprintf_s(stderr, _TEXT("done: %02u時間%02u分%02u秒"),
                      std::size_t(times / (60 * 60)),
                      std::size_t(times / 60) % 60, std::size_t(times) % 60);
      }
      std::tcerr << std::endl;

      // buil kd-tree
      std::tcerr << _TEXT("build kd-tree ... ");
      {
        std::clock_t begin_times = std::clock();
        Scene::GetInstance().Build();

        std::clock_t times = std::clock() - begin_times;
        times = static_cast<std::clock_t>(times / double(CLOCKS_PER_SEC) + 0.5);
        ::_ftprintf_s(stderr, _TEXT("done: %02u時間%02u分%02u秒"),
                      std::size_t(times / (60 * 60)),
                      std::size_t(times / 60) % 60, std::size_t(times) % 60);
      }
    }

    return SHELL_CONTINUE;
  }

  if (_TEXT("save") == command) {
    sin >> command;

    hi::mkdir(RENDIR_PATH);

    std::tstring filename(RENDIR_PATH);
    filename += command;
    filename += _TEXT(".ren");

    std::wofstream fout(filename.c_str());
    if (!!fout) {
      Scene const &scene = Scene::GetInstance();
      float3_t const eye = scene.ThinLensCamera().Eye();
      float3_t const at = scene.ThinLensCamera().At();

      fout << _TEXT("load ") << s_scene_name << std::endl;
      fout << std::endl;
      fout << _TEXT("set eye ") << eye[0] << _TEXT(' ') << eye[1] << _TEXT(' ')
           << eye[2] << std::endl;
      fout << _TEXT("set at  ") << at[0] << _TEXT(' ') << at[1] << _TEXT(' ')
           << at[2] << std::endl;
      fout << std::endl;
      fout << _TEXT("set dim ") << scene.GetWidth() << _TEXT(' ')
           << scene.GetHeight() << std::endl;
      fout << _TEXT("set it  ") << ::time2string(s_interval_time) << std::endl;
      fout << _TEXT("set rt  ") << ::time2string(s_rendering_time) << std::endl;
      fout << std::endl;
      fout << _TEXT("set rep ") << s_algorithm_type << std::endl;
      fout << _TEXT("ren") << std::endl;
    }

    return SHELL_CONTINUE;
  }

  if (_TEXT("ren") == command) {
    if (s_rendering_time <= 0) {
      return SHELL_CONTINUE;
    }

    // 出力フォルダを生成
    {
      std::tchar_t filename[MAX_PATH];

      std::tchar_t modulepath[MAX_PATH];
      ::GetModuleFileName(NULL, modulepath, MAX_PATH);

      std::tstring modulename;
      hi::path2name(std::tstring(modulepath), modulename);

      std::time_t time;
      std::time(&time);
      std::tm t;
      ::localtime_s(&t, &time);

      ::_stprintf_s(
          filename, MAX_PATH,
          OUTDIR_PATH _TEXT("%s/%s %04dx%04d %04d-%02d-%02dT%02d.%02d.%02d "),
          modulename.c_str(), s_scene_name.c_str(),
          Scene::GetInstance().GetWidth(), Scene::GetInstance().GetHeight(),
          1900 + t.tm_year, 1 + t.tm_mon, t.tm_mday, t.tm_hour, t.tm_min,
          t.tm_sec);

      s_outdir_path = filename;
    }

    //
    // rendering
    //
    switch (s_algorithm_type) {
    case 0:
      s_kicker = std::auto_ptr<hi::thread>(new PathTracer());
      s_outdir_path += _TEXT("PT/");
      break;
    case 1:
      s_kicker = std::auto_ptr<hi::thread>(new ImportanceDrivenPathTracer());
      s_outdir_path += _TEXT("IDPT/");
      break;
    case 2:
      s_kicker =
          std::auto_ptr<hi::thread>(new PathTracer_MetropolisLightTransport());
      s_outdir_path += _TEXT("PT_MLT/");
      break;
    case 3:
      s_kicker = std::auto_ptr<hi::thread>(
          new ImportanceDrivenPathTracer_MetropolisLightTransport());
      s_outdir_path += _TEXT("IDPT_MLT/");
      break;

    case 4:
      s_kicker = std::auto_ptr<hi::thread>(
          new PathTracerWithGoWithTheWinnersStrategy());
      s_outdir_path += _TEXT("PTwGW/");
      break;
    case 5:
      s_kicker = std::auto_ptr<hi::thread>(
          new ImportanceDrivenPathTracerWithGoWithTheWinnersStrategy());
      s_outdir_path += _TEXT("IDPTwGW/");
      break;
    case 6:
      s_kicker = std::auto_ptr<hi::thread>(
          new PathTracerWithGoWithTheWinnersStrategy_MetropolisLightTransport());
      s_outdir_path += _TEXT("PTwGW_MLT/");
      break;
    case 7:
      s_kicker = std::auto_ptr<hi::thread>(
          new ImportanceDrivenPathTracerWithGoWithTheWinnersStrategy_MetropolisLightTransport());
      s_outdir_path += _TEXT("IDPTwGW_MLT/");
      break;

    case 8:
      s_kicker = std::auto_ptr<hi::thread>(new BidirectionalPathTracer());
      s_outdir_path += _TEXT("BPT/");
      break;
    case 9:
      s_kicker = std::auto_ptr<hi::thread>(
          new ImportanceDrivenBidirectionalPathTracer());
      s_outdir_path += _TEXT("IDBPT/");
      break;
    case 10:
      s_kicker = std::auto_ptr<hi::thread>(
          new BidirectionalPathTracer_MetropolisLightTransport());
      s_outdir_path += _TEXT("BPT_MLT/");
      break;
    case 11:
      s_kicker = std::auto_ptr<hi::thread>(
          new ImportanceDrivenBidirectionalPathTracer_MetropolisLightTransport());
      s_outdir_path += _TEXT("IDBPT_MLT/");
      break;

    case 12:
      s_kicker =
          std::auto_ptr<hi::thread>(new MetropolisLightTransportTestbed());
      s_outdir_path += _TEXT("MLT/");
      break;

    case 13:
      s_kicker = std::auto_ptr<hi::thread>(
          new PathTracer_PopulationLightTransportTestbed());
      s_outdir_path += _TEXT("PT_PLT/");
      break;
    case 14:
      s_kicker = std::auto_ptr<hi::thread>(
          new BidirectionalPathTracer_PopulationLightTransportTestbed());
      s_outdir_path += _TEXT("BPT_PLT/");
      break;

    default:
      return SHELL_CONTINUE;
    }

    if (s_algorithm_type >= 0) {
      std::tcerr << _TEXT("mkdir: `") << s_outdir_path << _TEXT("'")
                 << std::endl;
      if (hi::mkdir(s_outdir_path.c_str())) {
        s_is_rendering = true;
        s_kicker->start();
      } else {
        std::tcerr << _TEXT("出力ディレクトリ(") << s_outdir_path
                   << _TEXT(")の作成に失敗しました．") << std::endl;
        return SHELL_CONTINUE;
      }
    } else {
      s_is_rendering = true;
      s_kicker->start();
    }
    if (wait_shell) {
      s_kicker->join();
      s_is_rendering = false;
    }

    return SHELL_QUIT;
  }

  return SHELL_UNKOWN;
}

int Shell(std::tistream &in, bool const gui) {
  std::tstring line;

  for (;;) {
    std::tcerr << std::endl << _TEXT("[") << RENDERER_NAME << ("]> ");

    do {
      if (!in) {
        std::tcerr << _TEXT("Bye Bye") << std::endl;
        std::exit(0);
      }

      std::getline(in, line);
      hi::trim(line);

      if (_TEXT("=begin") == line) {
        do {
          if (!in) {
            std::tcerr << _TEXT("bye-bye") << std::endl;
            std::exit(0);
          }

          std::getline(in, line);
          hi::trim(line);
        } while (_TEXT("=end") != line);

        line.clear();
      }
    } while (line.empty() || hi::starts_with(line, _TEXT("rem ")));

    switch (Interpret(line, !gui)) {
    case SHELL_EXIT:  // exit system
      std::exit(0);
      break;
    case SHELL_QUIT:  // exit shell
      if (gui) {
        return 0;
      }
    case SHELL_CONTINUE:  // continue
      break;
    default:  // undefined
      std::tcerr << _TEXT("undefined") << std::endl;
      break;
    }

    std::tcerr << std::endl;
  }

  return 0;
}

}  // end of namespace xyz

#ifdef IMPORTANCE_DRIVEN
//
// init the photon map
//
std::tcerr << _TEXT("init the photon map") << std::endl;
{
  s_scene.BeginEmitPhotons(20000);  // とりあえず2万個蓄える

  std::tcerr << _TEXT("emit photons ... ");
  s_scene.EmitPhotons(20000);  // とりあえず2万個射出する
  std::tcerr << _TEXT("done.") << std::endl;

  std::tcerr << _TEXT("build the photon map ... ");
  s_scene.EndEmitPhotons();  // フォトンマップの構築
  std::tcerr << _TEXT("done.") << std::endl;
}
s_scene.ShowPhotonMapInfo();
#endif
