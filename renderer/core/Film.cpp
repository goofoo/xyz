#include "Film.hpp"
#include "geom/Path.hpp"

namespace tgir {
void Film::Clear() { x_.assign(x_.size(), Film::value_type(0)); }

void Film::Resize(std::size_t const w, std::size_t const h) {
  width_ = w;
  height_ = h;
  x_.assign(w * h, Film::value_type(0));
}

void Film::Deposite(__in tgir::PixelDescriptor const &value,
                    __in tgir::Real const &w) {
  if (w <= tgir::Real(0)) {
    return;
  }

  union local_type {
    double d;
    hi::qword q;
  };

#ifdef CONFIG_CONCURRENT_LOCK
  if (tgir::Path::INVALID_CONSTRIBUTION != value.first) {
    tgir::Film::value_type const xyz(value.second * w);
    tgir::Film::value_type &_val = x_[value.first];
    local_type _old;
    local_type _new;
    for (std::size_t c = 0; c < 3; ++c) {
      _old.d = _val[c];
      do {
        _new.d = _old.d + xyz[c];
      } while (
          !hi::qword_cmpxchg(reinterpret_cast<hi::qword volatile &>(_val[c]),
                             _old.q, _new.q));
    }
  }
#else
  if (tgir::Path::INVALID_CONSTRIBUTION != value.first) {
    hi::synchronized lock(monitor_);
    tgir::Film::value_type const xyz(value.second * w);
    x_[value.first] += xyz;
  }
#endif
}

void Film::Deposite(__in std::vector<tgir::PixelDescriptor> const &values,
                    __in tgir::Real const &w) {
  if (w <= tgir::Real(0)) {
    return;
  }

  union local_type {
    double d;
    hi::qword q;
  };

#ifdef CONFIG_CONCURRENT_LOCK
  for (std::size_t i = 0, size = values.size(); i < size; ++i) {
    if (tgir::Path::INVALID_CONSTRIBUTION != values[i].first) {
      tgir::Film::value_type const xyz(values[i].second * w);
      tgir::Film::value_type &_val = x_[values[i].first];
      local_type _old;
      local_type _new;
      for (std::size_t c = 0; c < 3; ++c) {
        _old.d = _val[c];
        do {
          _new.d = _old.d + xyz[c];
        } while (
            !hi::qword_cmpxchg(reinterpret_cast<hi::qword volatile &>(_val[c]),
                               _old.q, _new.q));
      }
    }
  }
#else
  hi::synchronized lock(monitor_);
  for (std::size_t i = 0, size = values.size(); i < size; ++i) {
    if (tgir::Path::INVALID_CONSTRIBUTION != values[i].first) {
      tgir::Film::value_type const xyz(values[i].second * w);
      x_[values[i].first] += xyz;
    }
  }
#endif
}

bool Film::Save(std::tchar_t const *filename, tgir::Real const mpp) const {
  if (hi::ends_with(filename, _TEXT(".float"))) {
    return SaveAsFloat(filename, mpp);
  }

  if (hi::ends_with(filename, _TEXT(".pfm"))) {
    return SaveAsPfm(filename, mpp);
  }
  /*
      if (hi::ends_with(filename, _TEXT(".hdr")))
      {
        return SaveAsHdr(filename, mpp);
      }
  */
  return false;
}

// mpp: mutation per. pixels
bool Film::SaveAsFloat(std::tchar_t const *filename,
                       tgir::Real const mpp) const {
  std::ofstream out(filename, std::ios::out | std::ios::binary);

  if (!out) {
    return false;
  }

  tgir::Real const s = hi::rcp(mpp);

  hi::basic_vector3<float> xyz;
  hi::basic_vector3<float> rgb;

  std::size_t ignore_pixel = 0;

  for (std::size_t p = 0, size = x_.size(); p < size; ++p) {
    xyz = x_[p];
    xyz *= s;

    hi::CCIR601_1_XYZ2RGB(xyz, rgb);

    /*
          if (!::_finite(rgb[0])
            || !::_finite(rgb[1])
            || !::_finite(rgb[2]))
          {
            rgb = 0; // 無限大は0に!!
            ++ignore_pixel;
          }
    */

    for (int c = 0; c < 3; ++c) {
      if (rgb[c] < 0) {
        rgb[c] = 0;
      }
    }

    out.write(reinterpret_cast<char *>(rgb.data()), 3 * sizeof(float));
  }
  std::tcerr << _TEXT(" ignored pixels: ") << ignore_pixel << std::endl;

  return true;
}

bool Film::SaveAsPfm(std::tchar_t const *filename, tgir::Real const mpp) const {
  std::ofstream out(filename, std::ios::out | std::ios::binary);
  if (!out) {
    return false;
  }
  out.imbue(std::locale("C"));

  out << "PF\n" << GetWidth() << " " << GetHeight()
      << "\n"
         "-1\n";  // ピクセルのアスペクト比(=1), マイナスなのはデータが Little
                  // Endian だから

  tgir::Real const s = hi::rcp(mpp);

  hi::basic_vector3<float> xyz;
  hi::basic_vector3<float> rgb;

  std::size_t ignore_pixel = 0;

  for (std::size_t p = 0, size = x_.size(); p < size; ++p) {
    xyz = x_[p];
    xyz *= s;

    hi::CCIR601_1_XYZ2RGB(xyz, rgb);

    if (!::_finite(rgb[0]) || !::_finite(rgb[1]) || !::_finite(rgb[2])) {
      rgb = 0;  // 無限大は0に!!
      ++ignore_pixel;
    }

    /*
          for (int c = 0; c < 3; ++c)
          {
            if (rgb[c] < 0)
            {
              rgb[c] = 0;
            }
          }
    */

    out.write(reinterpret_cast<char *>(rgb.data()), 3 * sizeof(float));
  }

  if (ignore_pixel > 0) {
    std::tcerr << _TEXT(" 無効画素数: ") << ignore_pixel << std::endl;
  }

  return true;
}

bool Film::SaveAsPfmRgb(std::tchar_t const *filename,
                        tgir::Real const mpp) const {
  std::ofstream out(filename, std::ios::out | std::ios::binary);

  if (!out) {
    return false;
  }

  out << "PF\n" << GetWidth() << " " << GetHeight()
      << "\n"
         "-1\n";  // ピクセルのアスペクト比(=1), マイナスなのはデータが Little
                  // Endian だから

  tgir::Real const s = hi::rcp(mpp);

  hi::basic_vector3<float> rgb;

  for (std::size_t p = 0, size = x_.size(); p < size; ++p) {
    rgb = x_[p];
    rgb *= s;

    out.write(reinterpret_cast<char *>(rgb.data()), 3 * sizeof(float));
  }

  return true;
}

}  // end of namespace tgir

/*
    // tentfilter
    color_type const rgb(hi::spd2rgb(rad));

    u *= GetWidth()(); u -= 0.5; // [-0.5,w-0.5)
    v *= GetHeight()(); v -= 0.5; // [0-0.5,h-0.5)

    int x = static_cast<int>(std::floor(u));
    int y = static_cast<int>(std::floor(v));

    u -= x; // [0,1)
    v -= y; // [0,1)

    double s[2*2]; // sup

    s[0] = (1 - u) * (1 - v);
    s[1] = (    u) * (1 - v);
    s[2] = (1 - u) * (    v);
    s[3] = (    u) * (    v);

    if (x < 0) { s[0] = s[2] = 0; }
    else if ((x+1) >= int(GetWidth()())) { s[1] = s[3] = 0; }
    if (y < 0) { s[0] = s[1] = 0; }
    else if ((y+1) >= int(GetHeight()())) { s[2] = s[3] = 0; }

    double k = 0;
    for (std::size_t q = 0; q < 4; ++q)
    {
      k += s[q];
    }
    k = 1 / k;

    for (std::size_t q = 0, j = 0; j < 2; ++j)
    {
      for (std::size_t i = 0; i < 2; ++i, ++q)
      {
        if (s[q] > 0)
        {
          std::size_t p = (y + j) * GetWidth()() + (x + i);
          x_[p] += rgb * (w * k * s[q]);
        }
      }
    }
*/
