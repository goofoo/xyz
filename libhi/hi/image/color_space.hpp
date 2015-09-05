#pragma once

namespace hi {
template <typename T>
struct value_type_of {
  typedef typename T::value_type type;
};

template <typename T>
struct value_type_of<T *> {
  typedef T type;
};

template <typename T>
struct value_type_of<T const *> {
  typedef T type;
};

//------------------------------------------------------------------------------
//
// XYZ -> Yxy
//
//------------------------------------------------------------------------------
template <typename S, typename T>
inline void XYZ2Yxy(S const &xyz, T &yxy) {
  typedef typename hi::value_type_of<T>::type R;
  typedef typename hi::value_type_of<S>::type Q;
  Q const sum = xyz[0] + xyz[1] + xyz[2];
  if (sum > Q(0)) {
    yxy[0] = static_cast<R>(xyz[1]);
    yxy[1] = static_cast<R>(xyz[0] / sum);
    yxy[2] = static_cast<R>(xyz[1] / sum);
  } else {
    yxy[0] = yxy[1] = yxy[2] = R(0);
  }
}

template <typename S, typename T>
inline void Yxy2XYZ(S const &yxy, T &xyz) {
  typedef typename hi::value_type_of<T>::type R;
  typedef typename hi::value_type_of<S>::type Q;
  if (yxy[2] > Q(0)) {
    xyz[0] = static_cast<R>(yxy[0] * yxy[1] / yxy[2]);
    xyz[1] = static_cast<R>(yxy[0]);
    xyz[2] = static_cast<R>(yxy[0] * (Q(1) - yxy[1] - yxy[2]) / yxy[2]);
  } else {
    xyz[0] = xyz[1] = xyz[2] = R(0);
  }
}

//------------------------------------------------------------------------------
//
// CCIR601-1 (D65)
//
//------------------------------------------------------------------------------

template <typename S, typename T>
inline void CCIR601_1_RGB2XYZ(S const &rgb, T &xyz) {
  typedef typename hi::value_type_of<T>::type R;
  typedef typename hi::value_type_of<S>::type Q;
  xyz[0] = static_cast<R>(Q(0.6069) * rgb[0] + Q(0.1735) * rgb[1] +
                          Q(0.2003) * rgb[2]);
  xyz[1] = static_cast<R>(Q(0.2989) * rgb[0] + Q(0.5866) * rgb[1] +
                          Q(0.1145) * rgb[2]);
  xyz[2] = static_cast<R>(Q(0.0000) * rgb[0] + Q(0.0661) * rgb[1] +
                          Q(1.1162) * rgb[2]);
}

template <typename S, typename T>
inline void CCIR601_1_XYZ2RGB(S const &xyz, T &rgb) {
  typedef typename hi::value_type_of<T>::type R;
  typedef typename hi::value_type_of<S>::type Q;
  rgb[0] = static_cast<R>(Q(1.9104) * xyz[0] - Q(0.5338) * xyz[1] -
                          Q(0.2891) * xyz[2]);
  rgb[1] = static_cast<R>(Q(-0.9844) * xyz[0] + Q(1.9985) * xyz[1] -
                          Q(0.0279) * xyz[2]);
  rgb[2] = static_cast<R>(Q(0.0585) * xyz[0] - Q(0.1187) * xyz[1] +
                          Q(0.9017) * xyz[2]);
}

//------------------------------------------------------------------------------
//
// CCIR709
//
//------------------------------------------------------------------------------

template <typename S, typename T>
inline void CCIR709_RGB2XYZ(S const &rgb, T &xyz) {
  typedef typename hi::value_type_of<T>::type R;
  typedef typename hi::value_type_of<S>::type Q;
  xyz[0] =
      static_cast<R>(Q(0.412) * rgb[0] + Q(0.358) * rgb[1] + Q(0.180) * rgb[2]);
  xyz[1] =
      static_cast<R>(Q(0.213) * rgb[0] + Q(0.715) * rgb[1] + Q(0.072) * rgb[2]);
  xyz[2] =
      static_cast<R>(Q(0.019) * rgb[0] + Q(0.119) * rgb[1] + Q(0.950) * rgb[2]);
}

template <typename S, typename T>
inline void CCIR709_XYZ2RGB(S const &xyz, T &rgb) {
  typedef typename hi::value_type_of<T>::type R;
  typedef typename hi::value_type_of<S>::type Q;
  rgb[0] =
      static_cast<R>(Q(3.241) * xyz[0] - Q(1.537) * xyz[1] - Q(0.499) * xyz[2]);
  rgb[1] = static_cast<R>(Q(-0.969) * xyz[0] + Q(1.876) * xyz[1] +
                          Q(0.042) * xyz[2]);
  rgb[2] =
      static_cast<R>(Q(0.056) * xyz[0] - Q(0.204) * xyz[1] + Q(1.057) * xyz[2]);
}

//------------------------------------------------------------------------------
//
// ITU
//
//------------------------------------------------------------------------------

template <typename S, typename T>
inline void ITU_RGB2XYZ(S const &rgb, T &xyz) {
  typedef typename hi::value_type_of<T>::type R;
  typedef typename hi::value_type_of<S>::type Q;
  xyz[0] = static_cast<R>(Q(0.4305) * rgb[0] + Q(0.3415) * rgb[1] +
                          Q(0.1784) * rgb[2]);
  xyz[1] = static_cast<R>(Q(0.2220) * rgb[0] + Q(0.7067) * rgb[1] +
                          Q(0.0713) * rgb[2]);
  xyz[2] = static_cast<R>(Q(0.0202) * rgb[0] + Q(0.1295) * rgb[1] +
                          Q(0.9394) * rgb[2]);
}

template <typename S, typename T>
inline void ITU_XYZ2RGB(S const &xyz, T &rgb) {
  typedef typename hi::value_type_of<T>::type R;
  typedef typename hi::value_type_of<S>::type Q;
  rgb[0] = static_cast<R>(Q(3.0527) * xyz[0] - Q(1.3928) * xyz[1] -
                          Q(0.4759) * xyz[2]);
  rgb[1] = static_cast<R>(Q(-0.9689) * xyz[0] + Q(1.8756) * xyz[1] +
                          Q(0.0417) * xyz[2]);
  rgb[2] = static_cast<R>(Q(0.0585) * xyz[0] - Q(0.2286) * xyz[1] +
                          Q(1.0690) * xyz[2]);
}

}  // end of namespace tgir
