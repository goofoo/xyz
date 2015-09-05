#include "Scene.hpp"
#include "Bsdf/Lambert.hpp"
#include "Bsdf/Mirror.hpp"
#include "Bsdf/Glass.hpp"
#include "Bsdf/Metal.hpp"
#include "Bsdf/Light.hpp"

namespace {
struct MoqMaterial {
  MoqMaterial()
      : name(),
        shader(3),
        dif(0.8),
        amb(0.6),
        emi(0),
        spc(0),
        power(5),
        tex(),
        aplane(),
        bump() {
    col[0] = col[1] = col[2] = col[3] = 1;
  }

  std::tstring name;
  std::size_t shader;
  float col[4];
  float dif;
  float amb;
  float emi;
  float spc;
  float power;
  std::tstring tex;
  std::tstring aplane;
  std::tstring bump;
};

struct MoqFace {
  inline MoqFace()
      : v0(0), v1(0), v2(0), m(0), s0(0), t0(0), s1(1), t1(0), s2(0), t2(1) {}

  int v0, v1, v2;
  int m;
  float s0, t0;
  float s1, t1;
  float s2, t2;
};

void getstring(std::tistream &in, std::tstring &str) {
  std::getline(in, str, _TEXT('\"'));
  std::getline(in, str, _TEXT('\"'));
}

void skip_blanket(std::tistream &in) {
  std::size_t indent = 0;
  std::tstring line;

  do {
    std::getline(in, line);
    hi::trim(line);
    if (_TEXT('{') == line[line.length() - 1]) {
      ++indent;
    } else if (_TEXT('}') == line[line.length() - 1]) {
      --indent;
    }
  } while (indent > 0);
}

void LoadMqoMaterial(std::tistream &in, std::size_t const size,
                     std::vector<tgir::Bsdf const *> &bsdfs) {
  std::tstring line;
  std::tstring tag;
  std::tstring data;

  for (std::size_t i = 0; i < size; ++i) {
    MoqMaterial mat;

    std::getline(in, line);
    std::tistringstream lin(line);  // line input
    getstring(lin, mat.name);       // get name

    while (!lin.eof()) {
      std::getline(lin, tag, _TEXT('('));
      hi::trim(tag);

      if (tag.length() == 0) {
        break;
      }

      if (_TEXT("shader") == tag) {
        std::getline(lin, data, _TEXT(')'));
        std::tistringstream din(data);  // data input
        din >> mat.shader;
      } else if (_TEXT("col") == tag) {
        std::getline(lin, data, _TEXT(')'));
        std::tistringstream din(data);  // data input
        din >> mat.col[0] >> mat.col[1] >> mat.col[2] >> mat.col[3];
      } else if (_TEXT("dif") == tag) {
        std::getline(lin, data, _TEXT(')'));
        std::tistringstream din(data);  // data input
        din >> mat.dif;
      } else if (_TEXT("amb") == tag) {
        std::getline(lin, data, _TEXT(')'));
        std::tistringstream din(data);  // data input
        din >> mat.amb;
      } else if (_TEXT("emi") == tag) {
        std::getline(lin, data, _TEXT(')'));
        std::tistringstream din(data);  // data input
        din >> mat.emi;
      } else if (_TEXT("spc") == tag) {
        std::getline(lin, data, _TEXT(')'));
        std::tistringstream din(data);  // data input
        din >> mat.spc;
      } else if (_TEXT("power") == tag) {
        std::getline(lin, data, _TEXT(')'));
        std::tistringstream din(data);  // data input
        din >> mat.power;
      } else if (_TEXT("tex") == tag) {
        getstring(lin, mat.tex);
        std::getline(in, data, _TEXT(')'));
      } else if (_TEXT("aplane") == tag) {
        getstring(lin, mat.aplane);
        std::getline(in, data, _TEXT(')'));
      } else if (_TEXT("bump") == tag) {
        getstring(lin, mat.bump);
        std::getline(in, data, _TEXT(')'));
      }
    }

    // マテリアルの生成
    bsdfs.push_back(nullptr);

    if (1 == mat.emi) {
      bsdfs.back() = new tgir::Light();
    } else if (1 == mat.col[3]) {
      if ((1 == mat.spc) && (mat.power > 0)) {
        if (100 == mat.power) {
          // 鏡面
          bsdfs.back() = new tgir::Mirror();
        } else {
          // 金属面
          bsdfs.back() =
              (mat.name.find(_TEXT("mcc-dark-skin")) != std::tstring::npos) ? new tgir::Metal(
                                                                                  mat.power,
                                                                                  tgir::SpectrumData::
                                                                                      MccDarkSkin())
                                                                            : (mat.name.find(
                                                                                   _TEXT(
                                                                                       "mcc-light-skin")) !=
                                                                               std::tstring::npos)
                                                                                  ? new tgir::Metal(
                                                                                        mat.power,
                                                                                        tgir::SpectrumData::
                                                                                            MccLightSkin())
                                                                                  : (mat.name.find(_TEXT(
                                                                                         "mcc-blue-sky")) !=
                                                                                     std::
                                                                                         tstring::npos)
                                                                                        ? new tgir::
                                                                                              Metal(mat.power,
                                                                                                    tgir::SpectrumData::
                                                                                                        MccBlueSky())
                                                                                        : (mat.name.find(_TEXT(
                                                                                               "mcc-foliage")) != std::tstring::npos)
                                                                                              ? new tgir::
                                                                                                    Metal(mat.power, tgir::SpectrumData::MccFoliage())
                                                                                              : (mat.name.find(_TEXT(
                                                                                                     "mcc-blue-flower")) != std::tstring::
                                                                                                                                npos)
                                                                                                    ? new tgir::Metal(mat.power, tgir::SpectrumData::MccBlueFlower())
                                                                                                    : (mat.name.find(_TEXT(
                                                                                                           "mcc-bluish-"
                                                                                                           "green")) != std::tstring::
                                                                                                                            npos)
                                                                                                          ? new tgir::
                                                                                                                Metal(mat.power, tgir::
                                                                                                                                     SpectrumData::MccBluishGreen())
                                                                                                          : (mat.name.find(_TEXT(
                                                                                                                 "mcc-"
                                                                                                                 "orang"
                                                                                                                 "e")) !=
                                                                                                             std::
                                                                                                                 tstring::
                                                                                                                     npos)
                                                                                                                ? new tgir::Metal(
                                                                                                                      mat.power,
                                                                                                                      tgir::SpectrumData::
                                                                                                                          MccOrange())
                                                                                                                : (
                                                                                                                      mat.name
                                                                                                                          .find(
                                                                                                                              _TEXT(
                                                                                                                                  "mcc"
                                                                                                                                  "-pu"
                                                                                                                                  "rpl"
                                                                                                                                  "ish"
                                                                                                                                  "-bl"
                                                                                                                                  "u"
                                                                                                                                  "e")) !=
                                                                                                                      std::tstring::
                                                                                                                          npos)
                                                                                                                      ? new tgir::Metal(
                                                                                                                            mat.power,
                                                                                                                            tgir::SpectrumData::
                                                                                                                                MccPurplishBlue())
                                                                                                                      : (mat.name.find(_TEXT(
                                                                                                                             "m"
                                                                                                                             "c"
                                                                                                                             "c"
                                                                                                                             "-"
                                                                                                                             "m"
                                                                                                                             "o"
                                                                                                                             "d"
                                                                                                                             "e"
                                                                                                                             "r"
                                                                                                                             "a"
                                                                                                                             "t"
                                                                                                                             "e"
                                                                                                                             "-"
                                                                                                                             "r"
                                                                                                                             "e"
                                                                                                                             "d")) !=
                                                                                                                         std::tstring::
                                                                                                                             npos)
                                                                                                                            ? new tgir::Metal(mat.power, tgir::SpectrumData::
                                                                                                                                                             MccModerateRed())
                                                                                                                            : (mat.name.find(
                                                                                                                                   _TEXT(
                                                                                                                                       "mcc-purple")) !=
                                                                                                                               std::tstring::
                                                                                                                                   npos)
                                                                                                                                  ? new tgir::Metal(
                                                                                                                                        mat.power,
                                                                                                                                        tgir::SpectrumData::
                                                                                                                                            MccPurple())
                                                                                                                                  : (mat.name.find(
                                                                                                                                         _TEXT(
                                                                                                                                             "mcc-yellow-green")) !=
                                                                                                                                     std::
                                                                                                                                         tstring::
                                                                                                                                             npos)
                                                                                                                                        ? new tgir::Metal(
                                                                                                                                              mat.power,
                                                                                                                                              tgir::SpectrumData::
                                                                                                                                                  MccYellowGreen())
                                                                                                                                        : (mat.name.find(
                                                                                                                                               _TEXT(
                                                                                                                                                   "mcc-orange-yellow")) !=
                                                                                                                                           std::tstring::
                                                                                                                                               npos)
                                                                                                                                              ? new tgir::Metal(
                                                                                                                                                    mat.power,
                                                                                                                                                    tgir::SpectrumData::
                                                                                                                                                        MccOrangeYellow())
                                                                                                                                              : (mat
                                                                                                                                                     .name.find(_TEXT(
                                                                                                                                                         "mcc-blue")) !=
                                                                                                                                                 std::tstring::
                                                                                                                                                     npos)
                                                                                                                                                    ? new tgir::Metal(
                                                                                                                                                          mat.power,
                                                                                                                                                          tgir::SpectrumData::
                                                                                                                                                              MccBlue())
                                                                                                                                                    : (mat.name.find(
                                                                                                                                                           _TEXT(
                                                                                                                                                               "mcc-green")) != std::
                                                                                                                                                                                    tstring::
                                                                                                                                                                                        npos)
                                                                                                                                                          ? new tgir::Metal(
                                                                                                                                                                mat.power,
                                                                                                                                                                tgir::SpectrumData::
                                                                                                                                                                    MccGreen())
                                                                                                                                                          : (mat.name.find(
                                                                                                                                                                 _TEXT(
                                                                                                                                                                     "mcc-red")) !=
                                                                                                                                                             std::tstring::
                                                                                                                                                                 npos)
                                                                                                                                                                ? new tgir::Metal(mat
                                                                                                                                                                                      .power,
                                                                                                                                                                                  tgir::SpectrumData::
                                                                                                                                                                                      MccRed())
                                                                                                                                                                : (mat.name.find(
                                                                                                                                                                       _TEXT(
                                                                                                                                                                           "mcc-yellow")) !=
                                                                                                                                                                   std::tstring::
                                                                                                                                                                       npos)
                                                                                                                                                                      ? new tgir::
                                                                                                                                                                            Metal(mat.power,
                                                                                                                                                                                  tgir::SpectrumData::MccYellow())
                                                                                                                                                                      : (mat.name.find(
                                                                                                                                                                             _TEXT(
                                                                                                                                                                                 "mcc-magenta")) !=
                                                                                                                                                                         std::tstring::
                                                                                                                                                                             npos)
                                                                                                                                                                            ? new tgir::Metal(mat.power, tgir::SpectrumData::MccMagenta())
                                                                                                                                                                            : (mat.name.find(
                                                                                                                                                                                   _TEXT(
                                                                                                                                                                                       "mcc-cyan")) !=
                                                                                                                                                                               std::tstring::
                                                                                                                                                                                   npos)
                                                                                                                                                                                  ? new tgir::
                                                                                                                                                                                        Metal(
                                                                                                                                                                                            mat.power, tgir::
                                                                                                                                                                                                           SpectrumData::MccCyan())
                                                                                                                                                                                  : (mat.name.find(
                                                                                                                                                                                         _TEXT(
                                                                                                                                                                                             "mcc-white")) !=
                                                                                                                                                                                     std::tstring::
                                                                                                                                                                                         npos)
                                                                                                                                                                                        ? new tgir::Metal(mat.power, tgir::
                                                                                                                                                                                                                         SpectrumData::
                                                                                                                                                                                                                             MccWhite())
                                                                                                                                                                                        : (mat.name.find(
                                                                                                                                                                                               _TEXT(
                                                                                                                                                                                                   "mcc-neutral-80")) !=
                                                                                                                                                                                           std::tstring::
                                                                                                                                                                                               npos)
                                                                                                                                                                                              ? new tgir::Metal(
                                                                                                                                                                                                    mat.power,
                                                                                                                                                                                                    tgir::SpectrumData::
                                                                                                                                                                                                        MccNeutral80())
                                                                                                                                                                                              : (mat.name.find(
                                                                                                                                                                                                     _TEXT(
                                                                                                                                                                                                         "mcc-neutral-65")) !=
                                                                                                                                                                                                 std::tstring::
                                                                                                                                                                                                     npos)
                                                                                                                                                                                                    ? new tgir::
                                                                                                                                                                                                          Metal(
                                                                                                                                                                                                              mat.power,
                                                                                                                                                                                                              tgir::SpectrumData::
                                                                                                                                                                                                                  MccNeutral65())
                                                                                                                                                                                                    : (mat.name.find(_TEXT(
                                                                                                                                                                                                           "mcc-neutral-50")) != std::tstring::
                                                                                                                                                                                                                                     npos)
                                                                                                                                                                                                          ? new tgir::Metal(
                                                                                                                                                                                                                mat.power,
                                                                                                                                                                                                                tgir::SpectrumData::
                                                                                                                                                                                                                    MccNeutral50())
                                                                                                                                                                                                          : (mat.name.find(
                                                                                                                                                                                                                 _TEXT(
                                                                                                                                                                                                                     "mcc-neutral-35")) !=
                                                                                                                                                                                                             std::tstring::
                                                                                                                                                                                                                 npos)
                                                                                                                                                                                                                ? new tgir::Metal(
                                                                                                                                                                                                                      mat.power,
                                                                                                                                                                                                                      tgir::SpectrumData::
                                                                                                                                                                                                                          MccNeutral35())
                                                                                                                                                                                                                : (mat.name.find(
                                                                                                                                                                                                                       _TEXT(
                                                                                                                                                                                                                           "mcc-black")) !=
                                                                                                                                                                                                                   std::tstring::
                                                                                                                                                                                                                       npos)
                                                                                                                                                                                                                      ? new tgir::Metal(
                                                                                                                                                                                                                            mat.power, tgir::
                                                                                                                                                                                                                                           SpectrumData::MccBlack())
                                                                                                                                                                                                                      : /* otherwise */ nullptr;

          if (!bsdfs.back()) {
            // ガンマ補正をかけて，スペクトルに変換する
            tgir::Vector3 const rgb(std::pow(mat.col[0] * mat.spc, 2.2f),
                                    std::pow(mat.col[1] * mat.spc, 2.2f),
                                    std::pow(mat.col[2] * mat.spc, 2.2f));
            std::tcerr << _TEXT("glossy(") << mat.col[0] << _TEXT(" ")
                       << mat.col[1] << _TEXT(" ") << mat.col[2] << _TEXT(")")
                       << std::endl;

            tgir::SpectrumData spd;
            hi::rgb2spd(rgb, spd);
            hi::clamp(spd);

            bsdfs.back() = new tgir::Metal(mat.power, spd);
          }
        }
      } else {
        // 完全拡散反射面
        bsdfs
            .back() = (mat.name.find(_TEXT("mcc-dark-skin")) !=
                       std::tstring::npos)
                          ? new tgir::Lambert(tgir::SpectrumData::MccDarkSkin())
                          : (mat.name.find(_TEXT("mcc-light-skin")) !=
                             std::tstring::npos)
                                ? new tgir::Lambert(
                                      tgir::SpectrumData::MccLightSkin())
                                : (mat.name.find(_TEXT("mcc-blue-sky")) !=
                                   std::tstring::npos)
                                      ? new tgir::Lambert(
                                            tgir::SpectrumData::MccBlueSky())
                                      : (mat.name.find(_TEXT("mcc-foliage")) !=
                                         std::tstring::npos)
                                            ? new tgir::Lambert(
                                                  tgir::SpectrumData::
                                                      MccFoliage())
                                            : (mat.name.find(
                                                   _TEXT("mcc-blue-flower")) !=
                                               std::tstring::npos)
                                                  ? new tgir::Lambert(
                                                        tgir::SpectrumData::
                                                            MccBlueFlower())
                                                  : (mat.name.find(_TEXT(
                                                         "mcc-bluish-green")) !=
                                                     std::tstring::npos)
                                                        ? new tgir::Lambert(
                                                              tgir::SpectrumData::
                                                                  MccBluishGreen())
                                                        : (mat.name.find(_TEXT(
                                                               "mcc-orange")) !=
                                                           std::tstring::npos)
                                                              ? new tgir::Lambert(
                                                                    tgir::SpectrumData::
                                                                        MccOrange())
                                                              : (mat.name.find(
                                                                     _TEXT(
                                                                         "mcc-"
                                                                         "purpl"
                                                                         "ish-"
                                                                         "blu"
                                                                         "e")) !=
                                                                 std::
                                                                     tstring::
                                                                         npos)
                                                                    ? new tgir::Lambert(
                                                                          tgir::SpectrumData::
                                                                              MccPurplishBlue())
                                                                    : (mat.name.find(_TEXT(
                                                                           "mcc"
                                                                           "-mo"
                                                                           "der"
                                                                           "ate"
                                                                           "-re"
                                                                           "d")) !=
                                                                       std::
                                                                           tstring::
                                                                               npos)
                                                                          ? new tgir::Lambert(
                                                                                tgir::SpectrumData::
                                                                                    MccModerateRed())
                                                                          : (mat.name.find(
                                                                                 _TEXT(
                                                                                     "mcc-purple")) !=
                                                                             std::
                                                                                 tstring::
                                                                                     npos)
                                                                                ? new tgir::Lambert(
                                                                                      tgir::SpectrumData::
                                                                                          MccPurple())
                                                                                : (mat.name.find(
                                                                                       _TEXT(
                                                                                           "mcc-yellow-green")) !=
                                                                                   std::
                                                                                       tstring::
                                                                                           npos)
                                                                                      ? new tgir::Lambert(
                                                                                            tgir::SpectrumData::
                                                                                                MccYellowGreen())
                                                                                      : (mat.name.find(
                                                                                             _TEXT(
                                                                                                 "mcc-orange-yellow")) !=
                                                                                         std::
                                                                                             tstring::
                                                                                                 npos)
                                                                                            ? new tgir::Lambert(
                                                                                                  tgir::SpectrumData::
                                                                                                      MccOrangeYellow())
                                                                                            : (mat.name.find(
                                                                                                   _TEXT(
                                                                                                       "mcc-blue")) !=
                                                                                               std::
                                                                                                   tstring::
                                                                                                       npos)
                                                                                                  ? new tgir::Lambert(
                                                                                                        tgir::SpectrumData::
                                                                                                            MccBlue())
                                                                                                  : (mat.name.find(
                                                                                                         _TEXT(
                                                                                                             "mcc-green")) !=
                                                                                                     std::
                                                                                                         tstring::
                                                                                                             npos)
                                                                                                        ? new tgir::
                                                                                                              Lambert(tgir::
                                                                                                                          SpectrumData::MccGreen())
                                                                                                        : (mat.name.find(
                                                                                                               _TEXT(
                                                                                                                   "mcc-red")) !=
                                                                                                           std::
                                                                                                               tstring::
                                                                                                                   npos)
                                                                                                              ? new tgir::Lambert(
                                                                                                                    tgir::SpectrumData::
                                                                                                                        MccRed())
                                                                                                              : (mat.name.find(
                                                                                                                     _TEXT(
                                                                                                                         "mcc-yellow")) !=
                                                                                                                 std::
                                                                                                                     tstring::
                                                                                                                         npos)
                                                                                                                    ? new tgir::Lambert(
                                                                                                                          tgir::SpectrumData::
                                                                                                                              MccYellow())
                                                                                                                    : (mat.name.find(
                                                                                                                           _TEXT(
                                                                                                                               "mcc-magenta")) !=
                                                                                                                       std::tstring::
                                                                                                                           npos)
                                                                                                                          ? new tgir::Lambert(
                                                                                                                                tgir::SpectrumData::
                                                                                                                                    MccMagenta())
                                                                                                                          : (mat.name.find(
                                                                                                                                 _TEXT(
                                                                                                                                     "mcc-cyan")) !=
                                                                                                                             std::tstring::
                                                                                                                                 npos)
                                                                                                                                ? new tgir::Lambert(
                                                                                                                                      tgir::SpectrumData::
                                                                                                                                          MccCyan())
                                                                                                                                : (mat.name.find(
                                                                                                                                       _TEXT(
                                                                                                                                           "mcc-white")) !=
                                                                                                                                   std::tstring::
                                                                                                                                       npos)
                                                                                                                                      ? new tgir::Lambert(tgir::SpectrumData::
                                                                                                                                                              MccWhite())
                                                                                                                                      : (mat.name.find(
                                                                                                                                             _TEXT(
                                                                                                                                                 "mcc-neutral-80")) !=
                                                                                                                                         std::tstring::
                                                                                                                                             npos)
                                                                                                                                            ? new tgir::Lambert(
                                                                                                                                                  tgir::SpectrumData::
                                                                                                                                                      MccNeutral80())
                                                                                                                                            : (mat.name.find(_TEXT(
                                                                                                                                                   "mcc-neutral-65")) != std::tstring::
                                                                                                                                                                             npos)
                                                                                                                                                  ? new tgir::Lambert(
                                                                                                                                                        tgir::SpectrumData::
                                                                                                                                                            MccNeutral65())
                                                                                                                                                  : (mat.name.find(
                                                                                                                                                         _TEXT(
                                                                                                                                                             "mcc-neutral-50")) !=
                                                                                                                                                     std::tstring::
                                                                                                                                                         npos)
                                                                                                                                                        ? new tgir::Lambert(tgir::
                                                                                                                                                                                SpectrumData::MccNeutral50())
                                                                                                                                                        : (mat.name.find(
                                                                                                                                                               _TEXT(
                                                                                                                                                                   "mcc-neutral-35")) !=
                                                                                                                                                           std::tstring::
                                                                                                                                                               npos)
                                                                                                                                                              ? new tgir::Lambert(
                                                                                                                                                                    tgir::SpectrumData::
                                                                                                                                                                        MccNeutral35())
                                                                                                                                                              : (mat.name
                                                                                                                                                                     .find(
                                                                                                                                                                         _TEXT(
                                                                                                                                                                             "mcc-black")) !=
                                                                                                                                                                 std::tstring::
                                                                                                                                                                     npos)
                                                                                                                                                                    ? new tgir::Lambert(
                                                                                                                                                                          tgir::SpectrumData::
                                                                                                                                                                              MccBlack())
                                                                                                                                                                    : /* otherwise */ nullptr;

        if (!bsdfs.back()) {
          // ガンマ補正をかけて，スペクトルに変換する
          tgir::Vector3 const rgb(std::pow(mat.col[0] * mat.dif, 2.2f),
                                  std::pow(mat.col[1] * mat.dif, 2.2f),
                                  std::pow(mat.col[2] * mat.dif, 2.2f));
          std::tcerr << _TEXT("diffuse(") << mat.col[0] << _TEXT(" ")
                     << mat.col[1] << _TEXT(" ") << mat.col[2] << _TEXT(")")
                     << std::endl;

          tgir::SpectrumData spd;
          hi::rgb2spd(rgb, spd);
          hi::clamp(spd);

          bsdfs.back() = new tgir::Lambert(spd);
        }
      }
    } else {
      // ガラス
      bsdfs.back() = new tgir::Glass();
    }
  }
}

tgir::Real LoadMqoObject(std::tistream &in,
                         std::vector<tgir::Bsdf const *> const &bsdfs,
                         std::vector<tgir::Triangle::Vector3> &virtices,
                         std::vector<tgir::Triangle const *> &triangles,
                         std::vector<tgir::Triangle const *> &lights,
                         std::vector<tgir::Real> &cdf) {
  std::size_t const vn = virtices.size();
  std::tstring line;
  std::tstring tag;

  int vertex_count = 0;
  int triangle_count = 0;

  int visible = 0;
  int shading = 1;

  std::vector<MoqFace> data;

  std::tcerr << _TEXT("  load object ... ") << std::endl;
  for (;;) {
    in >> line;
    if (_TEXT("visible") == line) {
      in >> visible;
    } else if (_TEXT("shading") == line) {
      in >> shading;
    } else if (_TEXT("vertex") == line) {
      if (!visible) {
        skip_blanket(in);
        continue;
      }

      in >> vertex_count;

      virtices.reserve(virtices.size() + vertex_count);

      std::getline(in, line);  // {

      std::tcerr << _TEXT("    load vertices (") << vertex_count
                 << _TEXT("): 0%\r");
      int log_interval = std::max(1, vertex_count / 100);
      for (int i = 0; i < vertex_count; ++i) {
        std::getline(in, line);
        std::tistringstream sin(line);

        tgir::Triangle::Vector3 p;
        sin >> p;
        p *= tgir::Real(0.01);  // cm単位からm単位へ変換
        virtices.push_back(p);

        if (i % log_interval == 0) {
          std::tcerr << _TEXT("    load vertices (") << vertex_count
                     << _TEXT("): ") << (100 * (i + 1) / vertex_count)
                     << _TEXT("%\r");
          ;
        }
      }
      std::tcerr << _TEXT("    load vertices (") << vertex_count
                 << _TEXT("): 100%\r") << std::endl;

      std::getline(in, line);  // }
    } else if (_TEXT("face") == line) {
      if (!visible) {
        skip_blanket(in);
        continue;
      }

      in >> triangle_count;

      data.resize(triangle_count);

      std::getline(in, line);  // {

      std::tcerr << _TEXT("    load faces (") << triangle_count
                 << _TEXT("): 0%\r");
      int log_interval = std::max(1, triangle_count / 100);
      for (int i = 0; i < triangle_count; ++i) {
        std::getline(in, line);

        // 括弧()を空白に置き換え
        std::replace_if(line.begin(), line.end(),
                        std::bind1st(std::equal_to<std::tchar_t>(), _TEXT('(')),
                        _TEXT(' '));
        std::replace_if(line.begin(), line.end(),
                        std::bind1st(std::equal_to<std::tchar_t>(), _TEXT(')')),
                        _TEXT(' '));

        std::size_t n;
        std::tistringstream sin(line);
        sin >> n;

        if (3 == n) {
          for (std::size_t j = 0; j < 3; ++j) {
            sin >> tag;
            if (_TEXT("V") == tag) {
              sin >> data[i].v0 >> data[i].v1 >> data[i].v2;
            } else if (_TEXT("M") == tag) {
              sin >> data[i].m;
            } else if (_TEXT("UV") == tag) {
              sin >> data[i].s0 >> data[i].t0 >> data[i].s1 >> data[i].t1 >>
                  data[i].s2 >> data[i].t2;
            }
          }
        }

        if (i % log_interval == 0) {
          std::tcerr << _TEXT("    load faces (") << triangle_count
                     << _TEXT("): ") << (100 * (i + 1) / triangle_count)
                     << _TEXT("%\r");
        }
      }
      std::tcerr << _TEXT("    load faces (") << triangle_count
                 << _TEXT("): 100%") << std::endl;

      std::getline(in, line);  // }
    } else if (_TEXT("}") == line) {
      break;
    } else {
      skip_blanket(in);
    }
  }

  if (!visible) {
    std::tcerr << _TEXT("    unvisible.") << std::endl;
    return 0;
  }

  tgir::Real rAllSurfaceArea = 0;

  std::tcerr << _TEXT("    compute triangle basis ... ");
  triangles.reserve(triangles.size() + triangle_count);
  if (shading) {
    std::vector<tgir::Vector3> normals(vertex_count, tgir::Vector3(0));

    for (int i = 0; i < triangle_count; ++i) {
#if 1
      tgir::Vector3 const ea =
          virtices[vn + data[i].v1] - virtices[vn + data[i].v2];
      tgir::Vector3 const eb =
          virtices[vn + data[i].v2] - virtices[vn + data[i].v0];
      tgir::Vector3 const ec =
          virtices[vn + data[i].v0] - virtices[vn + data[i].v1];

      tgir::Real const a = hi::length(ea);
      tgir::Real const b = hi::length(eb);
      tgir::Real const c = hi::length(ec);

      tgir::Real const aa =
          std::max(0.0, 1 - M_1_PI * std::acos(hi::dot(eb, ec) / (b * c)));
      tgir::Real const ab =
          std::max(0.0, 1 - M_1_PI * std::acos(hi::dot(ec, ea) / (c * a)));
      tgir::Real const ac =
          std::max(0.0, 1 - M_1_PI * std::acos(hi::dot(ea, eb) / (a * b)));

#if 1
      tgir::Vector3 const n = hi::normalize(tgir::NormalVector(ec, -eb));

      normals[data[i].v0] += n * aa;
      normals[data[i].v1] += n * ab;
      normals[data[i].v2] += n * ac;
#else
      tgir::Vector3 const n = tgir::NormalVector(ec, -eb);

      normals[data[i].v0] += n * (aa / (aa + ab + ac));
      normals[data[i].v1] += n * (ab / (aa + ab + ac));
      normals[data[i].v2] += n * (ac / (aa + ab + ac));
#endif

#else
      tgir::Vector3 const e1 =
          virtices[vn + data[i].v1] - virtices[vn + data[i].v0];
      tgir::Vector3 const e2 =
          virtices[vn + data[i].v2] - virtices[vn + data[i].v0];

      tgir::Vector3 const n = hi::normalize(tgir::NormalVector(e1, e2));

      normals[data[i].v0] += n;
      normals[data[i].v1] += n;
      normals[data[i].v2] += n;
#endif
    }
#if 1
    for (int i = 0; i < vertex_count; ++i) {
      normals[i] = hi::normalize(normals[i]);
    }
#endif
    for (int i = 0; i < triangle_count; ++i) {
      triangles.push_back(nullptr);
      triangles.back() = new tgir::Triangle(
          data[i].m, virtices[data[i].v0 + vn], virtices[data[i].v1 + vn],
          virtices[data[i].v2 + vn], normals[data[i].v0], normals[data[i].v1],
          normals[data[i].v2], data[i].s0, data[i].t0, data[i].s1, data[i].t1,
          data[i].s2, data[i].t2);

      tgir::Real const rSurfaceArea = triangles.back()->area();
      rAllSurfaceArea += rSurfaceArea;

      if (bsdfs[data[i].m]->What() == tgir::Bsdf::LIGHT) {
#ifndef NDEBUG
        std::tcerr << _TEXT("Triangle[") << i << _TEXT("] as a light")
                   << std::endl;
#endif
        lights.push_back(triangles.back());
        cdf.push_back(cdf.back() + rSurfaceArea);
      }
    }
  } else {
    for (int i = 0; i < triangle_count; ++i) {
      triangles.push_back(nullptr);
      triangles.back() = new tgir::Triangle(
          data[i].m, virtices[vn + data[i].v0], virtices[vn + data[i].v1],
          virtices[vn + data[i].v2], data[i].s0, data[i].t0, data[i].s1,
          data[i].t1, data[i].s2, data[i].t2);

      tgir::Real const rSurfaceArea = triangles.back()->area();
      rAllSurfaceArea += rSurfaceArea;

      if (bsdfs[data[i].m]->What() == tgir::Bsdf::LIGHT) {
        lights.push_back(triangles.back());
        cdf.push_back(cdf.back() + rSurfaceArea);
      }
    }
  }
  std::tcerr << _TEXT("done.") << std::endl;    // compute Triangle basis
  std::tcerr << _TEXT("  done.") << std::endl;  // load object

  return rAllSurfaceArea;
}

}  // end of unnamed namespace

namespace tgir {

void Scene::Load(std::tistream &in) {
  Clear();

  std::tstring line;
  std::tstring tag;

  std::getline(in, line);
  if (_TEXT("Metasequoia Document") != line) {
    std::tcerr << _TEXT("MQOファイルでない") << std::endl;
    return;
  }

  std::getline(in, line);
  if (_TEXT("Format Text Ver 1.0") != line) {
    std::tcerr << _TEXT("対応していないフォーマット") << std::endl;
    return;
  }

  std::vector<tgir::Triangle::Vector3> vertices;

  while (!in.eof()) {
    in >> tag;
    if (_TEXT("Scene") == tag) {
      tgir::Vector3 position(0, 0, 400);
      tgir::Vector3 lookat(0, 100, 0);
      tgir::Real head = 0;
      tgir::Real pitch = 0;

      std::size_t indent = 0;
      do {
        std::getline(in, line);
        hi::trim(line);

        if (_TEXT('{') == line[line.length() - 1]) {
          ++indent;
        } else if (_TEXT('}') == line[line.length() - 1]) {
          --indent;
        } else {
          std::tistringstream din(line);
          din >> tag;
          if (_TEXT("pos") == tag) {
            din >> position;
          } else if (_TEXT("lookat") == tag) {
            din >> lookat;
          } else if (_TEXT("head") == tag) {
            din >> head;
          } else if (_TEXT("pich") == tag) {
            din >> pitch;
          }
        }
      } while (indent > 0);

      // pich変換(Y-Z)
      {
        tgir::Real const sin = std::sin(pitch);
        tgir::Real const cos = std::cos(pitch);
        tgir::Vector3 const v(position);
        position[1] = cos * v[1] + sin * v[2];
        position[2] = sin * v[1] - cos * v[2];
      }

      // head変換(Z-X)
      {
        head -= M_PI;
        tgir::Real const sin = std::sin(head);
        tgir::Real const cos = std::cos(head);
        tgir::Vector3 const v(position);
        position[2] = cos * v[2] + sin * v[0];
        position[0] = sin * v[2] - cos * v[0];
      }

      // 単位変換 (cm -> m)
      position *= 0.01;
      lookat *= 0.01;

      camera_.MoveTo(lookat + position);
      camera_.LookAt(lookat);
    }
    if (_TEXT("Material") == tag) {
      std::size_t materials;
      in >> materials;

      if (materials > 0) {
        std::getline(in, line);  // {
        LoadMqoMaterial(in, materials, bsdfs_);
        std::getline(in, line);  // }
      } else {
        bsdfs_.push_back(nullptr);
        bsdfs_.back() = new tgir::Lambert(tgir::SpectrumData::MccWhite());
        skip_blanket(in);
      }
    } else if (_TEXT("Object") == tag) {
      std::getline(in, line);  // {
      rSurfaceArea_ +=
          LoadMqoObject(in, bsdfs_, vertices, triangles_, lights_, lightsCdf_);
      std::getline(in, line);  // \n
    } else if (_TEXT("Eof") == tag) {
      break;
    }
  }

#ifndef NDEBUG
  std::size_t minMaterialIndex = std::numeric_limits<std::size_t>::max();
  std::size_t maxMaterialIndex = std::numeric_limits<std::size_t>::min();
  for (std::size_t i = 0, size = triangles_.size(); i < size; ++i) {
    minMaterialIndex = std::min(minMaterialIndex, triangles_[i]->bsdf());
    maxMaterialIndex = std::max(maxMaterialIndex, triangles_[i]->bsdf());
  }
  ::_ftprintf_s(stdout, _TEXT("material(%d) = [%u, %u]\n"), bsdfs_.size(),
                minMaterialIndex, maxMaterialIndex);
#endif NDEBUG

  std::tcerr << _TEXT("  SceneArea=") << rSurfaceArea_ << std::endl;
  std::tcerr << _TEXT("  LightArea=") << lightsCdf_.back() << std::endl;

  //
  // create display list
  //
  if (~0 == displayList_) {
    displayList_ = ::glGenLists(1);
  }

  ::glNewList(displayList_, GL_COMPILE);
  {
#if 0
      // 光源の設定
      if (!lights_.empty())
      {
        std::mt19937_64 random;
        tgir::Vector3 P, T, B, N, D;

        tgir::Vector3 xyz;
        hi::spd2xyz(tgir::GetLightPower(), xyz);
        xyz *= M_1_PI / (2 * 8);

        hi::basic_vector3<float> light_radiance;
        hi::CCIR601_1_XYZ2RGB(xyz.data(), light_radiance.data());

        for (int i = 0; i < 8; ++i)
        {
          tgir::Real u = std::uniform_real_distribution<tgir::Real>()(random);
          u *= lightsCdf_.back();
          std::vector<tgir::Real>::const_iterator const it =
            std::upper_bound(lightsCdf_.Begin(), lightsCdf_.End(), u);
          assert(it != lightsCdf_.Begin());
          assert(it != lightsCdf_.End());

          std::size_t n = it - lightsCdf_.Begin() - 1;
          lights_[n]->Sample(
            tgir::Vector2(
              std::uniform_real_distribution<tgir::Real>()(random),
              std::uniform_real_distribution<tgir::Real>()(random)),
              P, T, N, B, D);

          float const light_position [] = { float(P[0]), float(P[1]), float(P[2]), 1 };
          hi::basic_vector3<float> const light_direction(N);

          hi::glLight(GL_LIGHT0 + i, GL_POSITION, light_position);
        //hi::glLight(GL_LIGHT0 + i, GL_AMBIENT , light_ambient );
          hi::glLight(GL_LIGHT0 + i, GL_DIFFUSE , light_radiance.data());
          hi::glLight(GL_LIGHT0 + i, GL_SPECULAR, light_radiance.data());
          hi::glLight(GL_LIGHT0 + i, GL_CONSTANT_ATTENUATION , 0.0f);
          hi::glLight(GL_LIGHT0 + i, GL_LINEAR_ATTENUATION   , 0.0f);
          hi::glLight(GL_LIGHT0 + i, GL_QUADRATIC_ATTENUATION, 1.0f);
          hi::glLight(GL_LIGHT0 + i, GL_SPOT_DIRECTION, light_direction.data());
          hi::glLight(GL_LIGHT0 + i, GL_SPOT_CUTOFF   , 88.0f);
          hi::glLight(GL_LIGHT0 + i, GL_SPOT_EXPONENT ,  1.0f);
          ::glEnable (GL_LIGHT0 + i);
        }
      }
#else
    {
      GLfloat const diffuse[4] = {0.8f, 0.8f, 0.8f, 1.0f};
      GLfloat const ambient[4] = {0.2f, 0.2f, 0.2f, 1.0f};
      ::glEnable(GL_LIGHT0);
      hi::glLight(GL_LIGHT0, GL_DIFFUSE, diffuse);
      hi::glLight(GL_LIGHT0, GL_AMBIENT, ambient);
    }
#endif
      for
        each(tgir::Triangle const *t in triangles_) {
          bsdfs_[t->bsdf()]->Render();
          t->Render();
        }
  }
  ::glEndList();

#ifdef LIGHT_PROBE
  // 背景画像のロード
  {
    FILE *fp;
    if (::_tfopen_s(&fp, RESOURCE_PATH _TEXT("light_probe.float"),
                    _TEXT("rb"))) {
      return;
    }
    std::fread(&light_probe_[0], sizeof(float), 3 * LIGHT_PROBE * LIGHT_PROBE,
               fp);
    std::fclose(fp);
  }
#endif
}

}  // end of namespace tgir
