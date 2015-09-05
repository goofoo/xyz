#include <hi/image/spectrum.hpp>
#include <hi/image/color_space.hpp>
#include <hi/math.hpp>

namespace hi {
basic_spectrum::basic_spectrum() {}

basic_spectrum::basic_spectrum(
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
    value_type const _30) {
  data_[0] = _00;
  data_[1] = _01;
  data_[2] = _02;
  data_[3] = _03;
  data_[4] = _04;
  data_[5] = _05;
  data_[6] = _06;
  data_[7] = _07;
  data_[8] = _08;
  data_[9] = _09;
  data_[10] = _10;
  data_[11] = _11;
  data_[12] = _12;
  data_[13] = _13;
  data_[14] = _14;
  data_[15] = _15;
  data_[16] = _16;
  data_[17] = _17;
  data_[18] = _18;
  data_[19] = _19;
  data_[20] = _20;
  data_[21] = _21;
  data_[22] = _22;
  data_[23] = _23;
  data_[24] = _24;
  data_[25] = _25;
  data_[26] = _26;
  data_[27] = _27;
  data_[28] = _28;
  data_[29] = _29;
  data_[30] = _30;
}

basic_spectrum::basic_spectrum(value_type const a) { (*this) = a; }

basic_spectrum &basic_spectrum::operator=(value_type const a) {
  for (std::size_t i = 0; i < Size; ++i) {
    data_[i] = a;
  }
  return *this;
}

basic_spectrum::basic_spectrum(basic_spectrum const &rhs) { (*this) = rhs; }

basic_spectrum &basic_spectrum::operator=(basic_spectrum const &rhs) {
  for (std::size_t i = 0; i < Size; ++i) {
    (*this)[i] = rhs[i];
  }
  return *this;
}

basic_spectrum &basic_spectrum::operator+=(value_type const a) {
  for (std::size_t i = 0; i < Size; ++i) {
    (*this)[i] += a;
  }
  return *this;
}

basic_spectrum basic_spectrum::operator+(value_type const a) const {
  return basic_spectrum(*this) += a;
}

basic_spectrum &basic_spectrum::operator+=(basic_spectrum const &rhs) {
  for (std::size_t i = 0; i < Size; ++i) {
    (*this)[i] += rhs[i];
  }
  return *this;
}

basic_spectrum basic_spectrum::operator+(basic_spectrum const &rhs) const {
  return basic_spectrum(*this) += rhs;
}

basic_spectrum &basic_spectrum::operator-=(basic_spectrum const &rhs) {
  for (std::size_t i = 0; i < Size; ++i) {
    (*this)[i] -= rhs[i];
  }
  return *this;
}

basic_spectrum basic_spectrum::operator-(basic_spectrum const &rhs) const {
  return basic_spectrum(*this) -= rhs;
}

basic_spectrum &basic_spectrum::operator*=(value_type const a) {
  for (std::size_t i = 0; i < Size; ++i) {
    (*this)[i] *= a;
  }
  return *this;
}

basic_spectrum basic_spectrum::operator*(value_type const a) const {
  return basic_spectrum(*this) *= a;
}

basic_spectrum &basic_spectrum::operator*=(basic_spectrum const &rhs) {
  for (std::size_t i = 0; i < Size; ++i) {
    (*this)[i] *= rhs[i];
  }
  return *this;
}

basic_spectrum basic_spectrum::operator*(basic_spectrum const &rhs) const {
  return basic_spectrum(*this) *= rhs;
}

basic_spectrum &basic_spectrum::operator/=(value_type const a) {
  for (std::size_t i = 0; i < Size; ++i) {
    (*this)[i] /= a;
  }
  return *this;
}

basic_spectrum basic_spectrum::operator/(value_type const a) const {
  return basic_spectrum(*this) /= a;
}

basic_spectrum &basic_spectrum::operator/=(basic_spectrum const &rhs) {
  for (std::size_t i = 0; i < Size; ++i) {
    (*this)[i] /= rhs[i];
  }
  return *this;
}

basic_spectrum basic_spectrum::operator/(basic_spectrum const &rhs) const {
  return basic_spectrum(*this) /= rhs;
}

basic_spectrum::value_type basic_spectrum::Sample(value_type u) const {
  u *= Size;
  u -= 0.5;

  if (u <= 0) {
    return (*this)[0];
  } else if (u >= Size - 1) {
    return (*this)[Size - 1];
  }

  std::size_t const k = static_cast<std::size_t>(u);
  u -= k;

  value_type const a = (*this)[k + 0];
  value_type const b = (*this)[k + 1];

  return a + (b - a) * u;
}

}  // end of namespace

//
// lights
//
namespace hi {
basic_spectrum const &basic_spectrum::LightSourceA() {
  static basic_spectrum const spd(
      14.71, 17.78, 21.00, 24.67, 28.70, 33.09, 37.82, 42.87, 48.25, 53.91,
      59.86, 66.06, 72.50, 79.13, 85.95, 92.91, 100.00, 107.18, 114.44, 121.73,
      129.04, 136.34, 143.62, 150.83, 157.98, 165.03, 171.96, 178.77, 185.43,
      191.93, 198.26);
  return spd;
}

basic_spectrum const &basic_spectrum::LightSourceB() {
  static basic_spectrum const spd(
      41.3, 52.1, 63.2, 73.1, 80.8, 85.4, 88.3, 92.0, 95.2, 96.5, 94.2, 90.7,
      89.5, 92.2, 96.9, 101.0, 102.8, 102.6, 101.0, 99.2, 98.0, 98.5, 99.7,
      101.0, 102.2, 103.9, 105.0, 104.9, 103.9, 101.6, 99.1);
  return spd;
}

basic_spectrum const &basic_spectrum::LightSourceC() {
  static basic_spectrum const spd(
      63.3, 80.6, 98.1, 112.4, 121.5, 124.0, 123.1, 123.8, 123.9, 120.7, 112.1,
      102.3, 96.9, 98.0, 102.1, 105.2, 105.3, 102.3, 97.8, 93.2, 89.7, 88.4,
      88.1, 88.0, 87.8, 88.2, 87.9, 86.3, 84.0, 80.2, 76.3);
  return spd;
}

basic_spectrum const &basic_spectrum::LightD65() {
  static basic_spectrum const spd(
      82.8, 91.5, 93.4, 86.7, 104.9, 117.0, 117.8, 114.9, 115.9, 108.8, 109.4,
      107.8, 104.8, 107.7, 104.4, 104.0, 100.0, 96.3, 95.8, 88.7, 90.0, 89.6,
      87.7, 83.3, 83.7, 80.0, 80.2, 82.3, 78.3, 69.7, 71.6);
  return spd;
}

basic_spectrum const &basic_spectrum::LightD100() {
  static basic_spectrum const spd(
      138.8, 151.5, 150.3, 134.6, 151.8, 162.8, 159.4, 150.3, 146.9, 134.4,
      130.0, 124.6, 115.6, 115.3, 109.7, 106.4, 100.0, 94.3, 91.4, 84.8, 83.5,
      81.6, 78.0, 72.6, 71.6, 68.3, 67.3, 67.1, 63.8, 56.7, 57.3);
  return spd;
}

basic_spectrum const &basic_spectrum::LightStandardWarmWhite() {
  static basic_spectrum const spd(
      7.9, 10.1, 12.1, 14.3, 16.5, 18.3, 19.5, 20.4, 20.2, 20.7, 20.0, 20.1,
      23.2, 31.8, 47.3, 71.3, 100.0, 121.2, 135.8, 135.8, 122.4, 102.4, 82.0,
      63.7, 48.1, 34.3, 25.4, 18.1, 12.4, 7.6, 4.0);
  return spd;
}

basic_spectrum const &basic_spectrum::LightWhite() {
  static basic_spectrum const spd(
      11.2, 14.9, 18.1, 21.4, 24.7, 28.0, 29.6, 31.2, 29.9, 31.3, 30.2, 29.6,
      32.8, 41.1, 55.3, 76.8, 100.0, 118.3, 127.7, 124.8, 110.8, 90.4, 74.5,
      57.2, 44.5, 33.0, 24.0, 16.5, 10.5, 6.0, 2.5);
  return spd;
}

basic_spectrum const &basic_spectrum::LightStandardWhite() {
  static basic_spectrum const spd(
      17.7, 23.1, 27.4, 33.4, 38.7, 44.3, 46.7, 50.2, 49.7, 50.0, 48.4, 45.7,
      47.4, 53.2, 64.7, 81.6, 100.0, 115.0, 120.8, 116.6, 103.1, 84.4, 69.8,
      52.7, 40.5, 30.5, 22.0, 15.5, 9.5, 5.0, 2.0);
  return spd;
}

basic_spectrum const &basic_spectrum::LightDaylight() {
  static basic_spectrum const spd(
      31.3, 42.1, 51.1, 61.5, 71.0, 81.4, 85.0, 89.1, 86.9, 86.0, 82.7, 75.5,
      74.2, 76.2, 80.6, 90.2, 100.0, 108.2, 108.8, 101.3, 88.0, 72.6, 58.1,
      46.8, 37.5, 27.9, 20.9, 14.2, 9.0, 4.1, 1.4);
  return spd;
}

basic_spectrum const &basic_spectrum::LightWarmWhiteDeLuxe() {
  static basic_spectrum const spd(
      13.1, 12.3, 12.3, 13.4, 14.5, 16.2, 17.1, 17.2, 17.2, 21.0, 27.4, 46.0,
      75.0, 87.5, 89.7, 92.4, 100.0, 112.7, 130.0, 150.0, 162.2, 169.3, 164.5,
      152.0, 133.0, 109.6, 88.0, 64.4, 42.5, 25.4, 14.0);
  return spd;
}

basic_spectrum const &basic_spectrum::LightSoftWhite() {
  static basic_spectrum const spd(
      21.8, 25.0, 28.4, 33.6, 38.6, 45.0, 47.2, 49.7, 49.7, 51.8, 53.0, 58.3,
      70.0, 76.3, 81.7, 88.8, 100.0, 114.0, 131.3, 150.0, 161.9, 166.8, 160.0,
      143.1, 125.0, 102.3, 80.9, 59.0, 39.0, 24.0, 13.5);
  return spd;
}

basic_spectrum const &basic_spectrum::LightCoolWhiteDeLuxe() {
  static basic_spectrum const spd(
      22.0, 25.6, 29.7, 35.0, 40.3, 46.5, 49.0, 51.5, 51.5, 55.6, 61.0, 79.0,
      100.0, 106.3, 102.2, 98.0, 100.0, 106.2, 117.3, 130.4, 138.5, 141.3,
      133.7, 120.7, 107.0, 87.7, 70.0, 52.7, 37.0, 23.0, 13.0);
  return spd;
}

basic_spectrum const &basic_spectrum::LightMercuryArcLamp() {
  static basic_spectrum const spd(
      46.2, 59.5, 46.2, 45.7, 71.9, 40.9, 30.2, 26.6, 11.7, 26.6, 22.2, 18.7,
      19.5, 23.9, 51.5, 72.8, 51.5, 45.3, 58.6, 47.1, 19.8, 24.9, 23.9, 23.1,
      23.2, 23.3, 23.1, 21.2, 17.7, 21.3, 25.9);
  return spd;
}
}  // end of namespace tgir

//
// Macbeth Color Checker
//
namespace hi {
basic_spectrum const &basic_spectrum::MccDarkSkin() {
  static basic_spectrum const spd(
      0.05, 0.05, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05,
      0.06, 0.06, 0.06, 0.07, 0.08, 0.10, 0.11, 0.13, 0.13, 0.14, 0.15, 0.16,
      0.18, 0.19, 0.21, 0.23, 0.25, 0.27, 0.29);
  return spd;
}

basic_spectrum const &basic_spectrum::MccLightSkin() {
  static basic_spectrum const spd(
      0.21, 0.21, 0.22, 0.23, 0.24, 0.26, 0.29, 0.31, 0.33, 0.34, 0.35, 0.34,
      0.31, 0.30, 0.31, 0.31, 0.32, 0.38, 0.47, 0.53, 0.57, 0.60, 0.62, 0.64,
      0.66, 0.68, 0.70, 0.73, 0.77, 0.80, 0.83);
  return spd;
}

basic_spectrum const &basic_spectrum::MccBlueSky() {
  static basic_spectrum const spd(
      0.33, 0.35, 0.35, 0.35, 0.34, 0.34, 0.32, 0.30, 0.29, 0.27, 0.25, 0.23,
      0.22, 0.20, 0.20, 0.18, 0.17, 0.17, 0.16, 0.15, 0.15, 0.14, 0.13, 0.12,
      0.12, 0.11, 0.11, 0.11, 0.10, 0.10, 0.10);
  return spd;
}

basic_spectrum const &basic_spectrum::MccFoliage() {
  static basic_spectrum const spd(
      0.03, 0.03, 0.03, 0.03, 0.03, 0.04, 0.04, 0.04, 0.05, 0.05, 0.08, 0.13,
      0.17, 0.17, 0.14, 0.12, 0.11, 0.10, 0.09, 0.08, 0.08, 0.08, 0.08, 0.08,
      0.08, 0.08, 0.08, 0.10, 0.13, 0.16, 0.18);
  return spd;
}

basic_spectrum const &basic_spectrum::MccBlueFlower() {
  static basic_spectrum const spd(
      0.43, 0.46, 0.46, 0.46, 0.45, 0.44, 0.42, 0.40, 0.38, 0.36, 0.33, 0.28,
      0.24, 0.22, 0.22, 0.21, 0.20, 0.21, 0.23, 0.24, 0.25, 0.25, 0.25, 0.27,
      0.33, 0.40, 0.47, 0.52, 0.55, 0.56, 0.57);
  return spd;
}

basic_spectrum const &basic_spectrum::MccBluishGreen() {
  static basic_spectrum const spd(
      0.32, 0.34, 0.35, 0.37, 0.38, 0.42, 0.48, 0.54, 0.59, 0.61, 0.60, 0.58,
      0.57, 0.53, 0.49, 0.44, 0.40, 0.36, 0.31, 0.27, 0.24, 0.23, 0.23, 0.22,
      0.21, 0.22, 0.22, 0.24, 0.25, 0.26, 0.26);
  return spd;
}

basic_spectrum const &basic_spectrum::MccOrange() {
  static basic_spectrum const spd(
      0.05, 0.04, 0.05, 0.04, 0.04, 0.05, 0.05, 0.05, 0.05, 0.06, 0.07, 0.09,
      0.12, 0.18, 0.23, 0.30, 0.39, 0.46, 0.52, 0.55, 0.56, 0.58, 0.59, 0.60,
      0.61, 0.61, 0.61, 0.62, 0.62, 0.63, 0.64);
  return spd;
}

basic_spectrum const &basic_spectrum::MccPurplishBlue() {
  static basic_spectrum const spd(
      0.30, 0.33, 0.35, 0.37, 0.38, 0.37, 0.34, 0.30, 0.24, 0.20, 0.16, 0.14,
      0.12, 0.11, 0.10, 0.09, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.09,
      0.10, 0.11, 0.11, 0.11, 0.10, 0.10, 0.11);
  return spd;
}

basic_spectrum const &basic_spectrum::MccModerateRed() {
  static basic_spectrum const spd(
      0.15, 0.15, 0.14, 0.14, 0.13, 0.13, 0.13, 0.13, 0.11, 0.11, 0.10, 0.10,
      0.09, 0.10, 0.10, 0.11, 0.12, 0.17, 0.29, 0.42, 0.51, 0.56, 0.59, 0.60,
      0.60, 0.60, 0.60, 0.60, 0.60, 0.61, 0.61);
  return spd;
}

basic_spectrum const &basic_spectrum::MccPurple() {
  static basic_spectrum const spd(
      0.20, 0.19, 0.18, 0.16, 0.13, 0.11, 0.09, 0.08, 0.07, 0.06, 0.05, 0.05,
      0.05, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.06, 0.08, 0.11, 0.14,
      0.17, 0.19, 0.23, 0.26, 0.31, 0.36, 0.40);
  return spd;
}

basic_spectrum const &basic_spectrum::MccYellowGreen() {
  static basic_spectrum const spd(
      0.05, 0.05, 0.05, 0.06, 0.07, 0.08, 0.10, 0.13, 0.19, 0.27, 0.38, 0.47,
      0.53, 0.54, 0.52, 0.50, 0.48, 0.45, 0.41, 0.37, 0.35, 0.34, 0.33, 0.33,
      0.32, 0.32, 0.33, 0.35, 0.35, 0.37, 0.38);
  return spd;
}

basic_spectrum const &basic_spectrum::MccOrangeYellow() {
  static basic_spectrum const spd(
      0.05, 0.05, 0.05, 0.05, 0.05, 0.06, 0.07, 0.08, 0.10, 0.12, 0.14, 0.18,
      0.26, 0.35, 0.44, 0.50, 0.55, 0.59, 0.61, 0.62, 0.62, 0.63, 0.64, 0.65,
      0.65, 0.66, 0.66, 0.66, 0.67, 0.68, 0.68);
  return spd;
}

basic_spectrum const &basic_spectrum::MccBlue() {
  static basic_spectrum const spd(
      0.20, 0.23, 0.26, 0.29, 0.32, 0.32, 0.28, 0.22, 0.16, 0.12, 0.09, 0.07,
      0.06, 0.05, 0.05, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04,
      0.04, 0.05, 0.05, 0.05, 0.05, 0.05, 0.06);
  return spd;
}

basic_spectrum const &basic_spectrum::MccGreen() {
  static basic_spectrum const spd(
      0.06, 0.06, 0.06, 0.07, 0.07, 0.08, 0.10, 0.12, 0.15, 0.18, 0.23, 0.31,
      0.36, 0.36, 0.33, 0.29, 0.25, 0.21, 0.17, 0.13, 0.10, 0.09, 0.08, 0.08,
      0.07, 0.07, 0.07, 0.08, 0.09, 0.09, 0.09);
  return spd;
}

basic_spectrum const &basic_spectrum::MccRed() {
  static basic_spectrum const spd(
      0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.04, 0.04, 0.04, 0.04, 0.04,
      0.04, 0.04, 0.04, 0.05, 0.05, 0.07, 0.11, 0.20, 0.33, 0.49, 0.61, 0.67,
      0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77);
  return spd;
}

basic_spectrum const &basic_spectrum::MccYellow() {
  static basic_spectrum const spd(
      0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.06, 0.10, 0.16, 0.26, 0.37, 0.47,
      0.56, 0.61, 0.64, 0.67, 0.70, 0.73, 0.75, 0.75, 0.75, 0.76, 0.77, 0.78,
      0.78, 0.78, 0.79, 0.79, 0.79, 0.80, 0.80);
  return spd;
}

basic_spectrum const &basic_spectrum::MccMagenta() {
  static basic_spectrum const spd(
      0.36, 0.37, 0.36, 0.34, 0.31, 0.28, 0.25, 0.22, 0.20, 0.18, 0.15, 0.13,
      0.11, 0.10, 0.11, 0.10, 0.11, 0.14, 0.22, 0.31, 0.41, 0.52, 0.60, 0.67,
      0.70, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77);
  return spd;
}

basic_spectrum const &basic_spectrum::MccCyan() {
  static basic_spectrum const spd(
      0.26, 0.27, 0.30, 0.33, 0.35, 0.39, 0.43, 0.44, 0.45, 0.43, 0.39, 0.34,
      0.29, 0.24, 0.19, 0.15, 0.12, 0.10, 0.09, 0.08, 0.07, 0.07, 0.07, 0.07,
      0.07, 0.07, 0.08, 0.07, 0.07, 0.07, 0.07);
  return spd;
}

basic_spectrum const &basic_spectrum::MccWhite() {
  static basic_spectrum const spd(
      0.67, 0.82, 0.87, 0.88, 0.88, 0.90, 0.90, 0.90, 0.91, 0.92, 0.91, 0.91,
      0.92, 0.92, 0.91, 0.91, 0.92, 0.93, 0.93, 0.92, 0.91, 0.92, 0.92, 0.92,
      0.92, 0.92, 0.91, 0.92, 0.92, 0.93, 0.92);
  return spd;
}

basic_spectrum const &basic_spectrum::MccNeutral80() {
  static basic_spectrum const spd(
      0.53, 0.58, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60,
      0.61, 0.61, 0.60, 0.60, 0.61, 0.61, 0.61, 0.61, 0.60, 0.60, 0.60, 0.60,
      0.60, 0.59, 0.59, 0.59, 0.59, 0.59, 0.59);
  return spd;
}

basic_spectrum const &basic_spectrum::MccNeutral65() {
  static basic_spectrum const spd(
      0.34, 0.35, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36,
      0.37, 0.37, 0.36, 0.36, 0.37, 0.37, 0.37, 0.36, 0.36, 0.36, 0.36, 0.36,
      0.36, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35);
  return spd;
}

basic_spectrum const &basic_spectrum::MccNeutral50() {
  static basic_spectrum const spd(
      0.19, 0.19, 0.19, 0.20, 0.20, 0.20, 0.20, 0.19, 0.20, 0.20, 0.20, 0.20,
      0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20,
      0.20, 0.20, 0.19, 0.19, 0.20, 0.20, 0.19);
  return spd;
}

basic_spectrum const &basic_spectrum::MccNeutral35() {
  static basic_spectrum const spd(
      0.08, 0.08, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09,
      0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09,
      0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09);
  return spd;
}

basic_spectrum const &basic_spectrum::MccBlack() {
  static basic_spectrum const spd(
      0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,
      0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,
      0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03);
  return spd;
}
}
// end of namespace hi

namespace hi {
basic_spectrum::value_type dot(basic_spectrum const &a,
                               basic_spectrum const &b) {
  basic_spectrum::value_type sum = basic_spectrum::value_type();
  for (std::size_t i = 0; i < basic_spectrum::Size; ++i) {
    sum += a[i] * b[i];
  }
  return sum;
}

void clamp(basic_spectrum &spd, basic_spectrum::value_type const &max,
           basic_spectrum::value_type const &min) {
  for (std::size_t i = 0; i < basic_spectrum::Size; ++i) {
    if (spd[i] < min) {
      spd[i] = min;
    } else if (spd[i] > max) {
      spd[i] = max;
    }
  }
}
}
// end of namespace hi

namespace hi {
/// <summary>
/// CIE XYZ color matching functions - Andrew S. Glassner, `Principles of
/// Digital Image Synthesis', 1995.
/// </summary>
basic_spectrum const XYZ_CMF[3] = {
    basic_spectrum(
        0.014 / 9.80328, 0.044 / 9.80328, 0.134 / 9.80328, 0.284 / 9.80328,
        0.348 / 9.80328, 0.336 / 9.80328, 0.291 / 9.80328, 0.195 / 9.80328,
        0.096 / 9.80328, 0.032 / 9.80328, 0.005 / 9.80328, 0.009 / 9.80328,
        0.063 / 9.80328, 0.166 / 9.80328, 0.290 / 9.80328, 0.433 / 9.80328,
        0.595 / 9.80328, 0.762 / 9.80328, 0.916 / 9.80328, 1.026 / 9.80328,
        1.062 / 9.80328, 1.003 / 9.80328, 0.854 / 9.80328, 0.642 / 9.80328,
        0.448 / 9.80328, 0.284 / 9.80328, 0.165 / 9.80328, 0.087 / 9.80328,
        0.047 / 9.80328, 0.023 / 9.80328, 0.011 / 9.80328),
    basic_spectrum(
        0.000 / 9.80328, 0.001 / 9.80328, 0.004 / 9.80328, 0.012 / 9.80328,
        0.023 / 9.80328, 0.038 / 9.80328, 0.060 / 9.80328, 0.091 / 9.80328,
        0.139 / 9.80328, 0.208 / 9.80328, 0.323 / 9.80328, 0.503 / 9.80328,
        0.710 / 9.80328, 0.862 / 9.80328, 0.954 / 9.80328, 0.995 / 9.80328,
        0.995 / 9.80328, 0.952 / 9.80328, 0.870 / 9.80328, 0.757 / 9.80328,
        0.631 / 9.80328, 0.503 / 9.80328, 0.381 / 9.80328, 0.265 / 9.80328,
        0.175 / 9.80328, 0.107 / 9.80328, 0.061 / 9.80328, 0.032 / 9.80328,
        0.017 / 9.80328, 0.008 / 9.80328, 0.004 / 9.80328),
    basic_spectrum(
        0.068 / 9.80328, 0.207 / 9.80328, 0.646 / 9.80328, 1.386 / 9.80328,
        1.747 / 9.80328, 1.772 / 9.80328, 1.669 / 9.80328, 1.288 / 9.80328,
        0.813 / 9.80328, 0.465 / 9.80328, 0.272 / 9.80328, 0.158 / 9.80328,
        0.078 / 9.80328, 0.042 / 9.80328, 0.020 / 9.80328, 0.009 / 9.80328,
        0.004 / 9.80328, 0.002 / 9.80328, 0.002 / 9.80328, 0.001 / 9.80328,
        0.001 / 9.80328, 0.000 / 9.80328, 0.000 / 9.80328, 0.000 / 9.80328,
        0.000 / 9.80328, 0.000 / 9.80328, 0.000 / 9.80328, 0.000 / 9.80328,
        0.000 / 9.80328, 0.000 / 9.80328, 0.000 / 9.80328),
};

/// <summary>
/// CIE RGB color matching functions -  Samuel J. Williamson and Herman Z.
/// Cummins, `Light and Color in Nature and Art', 1983.
/// </summary>
basic_spectrum const RGB_CMF[3] = {
    basic_spectrum(0.000, 0.001, 0.002, 0.002, -0.003, -0.012, -0.026, -0.039,
                   -0.049, -0.058, -0.072, -0.089, -0.093, -0.071, -0.032,
                   0.023, 0.091, 0.168, 0.245, 0.309, 0.344, 0.340, 0.297,
                   0.227, 0.160, 0.102, 0.059, 0.031, 0.017, 0.008, 0.004),
    basic_spectrum(0.000, 0.000, -0.001, -0.001, 0.001, 0.007, 0.015, 0.025,
                   0.039, 0.057, 0.085, 0.129, 0.175, 0.203, 0.215, 0.212,
                   0.197, 0.171, 0.136, 0.098, 0.062, 0.036, 0.018, 0.008,
                   0.003, 0.001, 0.000, 0.000, 0.000, 0.000, 0.000),
    basic_spectrum(0.012, 0.037, 0.115, 0.248, 0.312, 0.317, 0.298, 0.230,
                   0.145, 0.083, 0.048, 0.027, 0.012, 0.005, 0.001, -0.001,
                   -0.001, -0.001, -0.001, -0.001, 0.000, 0.000, 0.000, 0.000,
                   0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000),
};

void spd2xyz(basic_spectrum const &spd,
             hi::basic_vector3<basic_spectrum::value_type> &xyz) {
  xyz[0] = hi::dot(spd, XYZ_CMF[0]);
  xyz[1] = hi::dot(spd, XYZ_CMF[1]);
  xyz[2] = hi::dot(spd, XYZ_CMF[2]);
}

void spd2xyz(std::size_t const i, basic_spectrum::value_type const spd,
             hi::basic_vector3<basic_spectrum::value_type> &xyz) {
  basic_spectrum::value_type const value = spd * basic_spectrum::Size;
  xyz[0] = value * XYZ_CMF[0][i];
  xyz[1] = value * XYZ_CMF[1][i];
  xyz[2] = value * XYZ_CMF[2][i];
}

basic_spectrum::value_type XyzCmf(std::size_t const i) {
  // return (XYZ_CMF[0][i]+XYZ_CMF[1][i]+XYZ_CMF[2][i])/3;
  return XYZ_CMF[1][i] * basic_spectrum::value_type(9.80328);
}

void rgb2spd(spectrum_vector const &rgb, basic_spectrum &spd) {
  static hi::basic_vector3<basic_spectrum::value_type> const lambda_c(
      641, 508, 426);  // 基底関数の中心波長

  static basic_spectrum::value_type const sigma_min = 80;
  static basic_spectrum::value_type const sigma_max = 150;

  basic_spectrum::value_type const int12 =
      (0 == (rgb[0] + rgb[1])) ? 0
                               : std::abs(rgb[0] - rgb[1]) / (rgb[0] + rgb[1]);
  basic_spectrum::value_type const int23 =
      (0 == (rgb[2] + rgb[1])) ? 0
                               : std::abs(rgb[2] - rgb[1]) / (rgb[2] + rgb[1]);

  hi::basic_vector3<basic_spectrum::value_type> sigma;
  sigma[0] = hi::lerp(sigma_max, sigma_min, int12);
  sigma[1] = hi::lerp(sigma_max, sigma_min, int23);
  sigma[2] = std::min(sigma[0], sigma[1]);

  hi::basic_vector3<basic_spectrum::value_type> const inv_sigma(
      hi::rsqrt(basic_spectrum::value_type(2 * M_PI) * sigma[0]),
      hi::rsqrt(basic_spectrum::value_type(2 * M_PI) * sigma[1]),
      hi::rsqrt(basic_spectrum::value_type(2 * M_PI) * sigma[2]));

  // 基底ベクトル
  basic_spectrum f[3];
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t k = 0; k < basic_spectrum::Size; ++k) {
      basic_spectrum::value_type const lambda =
          static_cast<basic_spectrum::value_type>(400 + 10 * k);
      basic_spectrum::value_type const x = (lambda - lambda_c[i]) / sigma[i];
      f[i][k] = std::exp(hi::square_of(x) / -2) * inv_sigma[i];
    }
  }

  // 線形方程式を解く $ax=b$
  hi::basic_vector3<basic_spectrum::value_type> a[3];
  for (std::size_t i = 0; i < 3; ++i) {
    a[i][0] = hi::dot(f[i], XYZ_CMF[0]);
    a[i][1] = hi::dot(f[i], XYZ_CMF[1]);
    a[i][2] = hi::dot(f[i], XYZ_CMF[2]);
  }

  hi::basic_vector3<basic_spectrum::value_type> x;
  hi::CCIR601_1_RGB2XYZ(rgb, x);
  hi::SolveLinearSystem(a[0], a[1], a[2], x);

  spd = 0;
  for (std::size_t i = 0; i < 3; ++i) {
    f[i] *= x[i];
    spd += f[i];
  }
}
}
// end of namespace hi
