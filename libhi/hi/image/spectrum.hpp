#pragma once

#include <cstddef>
#include <cassert>

namespace hi {
/**
 * スペクトルデータ (400-700nmの波長を10nm間隔で31サンプル)
 */
//[TODO]冗長すぎるのでどうにかしたい
class basic_spectrum {
 public:
  typedef float value_type;
  static std::size_t const Size = 31;

 public:
  basic_spectrum();
  basic_spectrum(basic_spectrum const &rhs);

  explicit basic_spectrum(value_type const a);
  explicit basic_spectrum(
      value_type const _00, value_type const _01, value_type const _02,
      value_type const _03, value_type const _04, value_type const _05,
      value_type const _06, value_type const _07, value_type const _08,
      value_type const _09, value_type const _10, value_type const _11,
      value_type const _12, value_type const _13, value_type const _14,
      value_type const _15, value_type const _16, value_type const _17,
      value_type const _18, value_type const _19, value_type const _20,
      value_type const _21, value_type const _22, value_type const _23,
      value_type const _24, value_type const _25, value_type const _26,
      value_type const _27, value_type const _28, value_type const _29,
      value_type const _30);

  basic_spectrum &operator=(basic_spectrum const &rhs);
  basic_spectrum &operator-=(basic_spectrum const &rhs);
  basic_spectrum &operator+=(basic_spectrum const &rhs);
  basic_spectrum &operator*=(basic_spectrum const &rhs);
  basic_spectrum &operator/=(basic_spectrum const &rhs);

  basic_spectrum operator+(basic_spectrum const &rhs) const;
  basic_spectrum operator-(basic_spectrum const &rhs) const;
  basic_spectrum operator*(basic_spectrum const &rhs) const;
  basic_spectrum operator/(basic_spectrum const &rhs) const;

  basic_spectrum &operator=(value_type const a);
  basic_spectrum &operator+=(value_type const a);
  basic_spectrum &operator*=(value_type const a);
  basic_spectrum &operator/=(value_type const a);

  basic_spectrum operator+(value_type const a) const;
  basic_spectrum operator*(value_type const a) const;
  basic_spectrum operator/(value_type const a) const;

  inline value_type *data() { return data_; }
  inline value_type const *data() const { return data_; }
  inline value_type &operator()(std::size_t const i) {
    assert(i < Size);
    return data_[i];
  }
  inline value_type operator()(std::size_t const i) const {
    assert(i < Size);
    return data_[i];
  }
  inline value_type &operator[](std::size_t const i) { return operator()(i); }
  inline value_type operator[](std::size_t const i) const {
    return operator()(i);
  }

  value_type Sample(value_type u) const;

  //
  // lights [W]
  //
  static basic_spectrum const &LightSourceA();
  static basic_spectrum const &LightSourceB();
  static basic_spectrum const &LightSourceC();
  static basic_spectrum const &LightD65();
  static basic_spectrum const &LightD100();
  static basic_spectrum const &LightStandardWarmWhite();
  static basic_spectrum const &LightWhite();
  static basic_spectrum const &LightStandardWhite();
  static basic_spectrum const &LightDaylight();
  static basic_spectrum const &LightWarmWhiteDeLuxe();
  static basic_spectrum const &LightSoftWhite();
  static basic_spectrum const &LightCoolWhiteDeLuxe();
  static basic_spectrum const &LightMercuryArcLamp();

  //
  // Macbeth Color Checker
  //
  static basic_spectrum const &MccDarkSkin();
  static basic_spectrum const &MccLightSkin();
  static basic_spectrum const &MccBlueSky();
  static basic_spectrum const &MccFoliage();
  static basic_spectrum const &MccBlueFlower();
  static basic_spectrum const &MccBluishGreen();
  static basic_spectrum const &MccOrange();
  static basic_spectrum const &MccPurplishBlue();
  static basic_spectrum const &MccModerateRed();
  static basic_spectrum const &MccPurple();
  static basic_spectrum const &MccYellowGreen();
  static basic_spectrum const &MccOrangeYellow();
  static basic_spectrum const &MccBlue();
  static basic_spectrum const &MccGreen();
  static basic_spectrum const &MccRed();
  static basic_spectrum const &MccYellow();
  static basic_spectrum const &MccMagenta();
  static basic_spectrum const &MccCyan();
  static basic_spectrum const &MccWhite();
  static basic_spectrum const &MccNeutral80();
  static basic_spectrum const &MccNeutral65();
  static basic_spectrum const &MccNeutral50();
  static basic_spectrum const &MccNeutral35();
  static basic_spectrum const &MccBlack();

 private:
  value_type data_[Size];
};
}  // end of namespace

namespace hi {
template <typename T>
class basic_vector3;

typedef hi::basic_vector3<basic_spectrum::value_type> spectrum_vector;

basic_spectrum::value_type dot(basic_spectrum const &a,
                               basic_spectrum const &b);
void clamp(basic_spectrum &spd, basic_spectrum::value_type const &max = 1,
           basic_spectrum::value_type const &min = 0);

void spd2xyz(basic_spectrum const &, spectrum_vector &);
void spd2xyz(std::size_t const, basic_spectrum::value_type const,
             spectrum_vector &);

basic_spectrum::value_type XyzCmf(std::size_t const i);

void rgb2spd(spectrum_vector const &, basic_spectrum &);

}  // end of namespace tgir
