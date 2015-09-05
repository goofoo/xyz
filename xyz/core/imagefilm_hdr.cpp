#include "imagefilm.hpp"

//
// Radiance (.hdr) 読み書き用ユーティリティ
//
namespace {
std::size_t const E_OFFSET = 128;
std::size_t const MIN_SCANLINE_LENGTH = 8;
std::size_t const MAX_SCANLINE_LENGTH = 0x7fff;  // 32767
std::size_t const MIN_RUNLENGTH = 4;

typedef hi::basic_vector4<hi::ubyte> xyze_t;

template <typename T>
inline void hdr_encode(hi::basic_vector3<T> const &src,
                       hi::basic_vector4<hi::ubyte> &dst) {
  T const d = hi::max(hi::max(src[0], src[1]), src[2]);

  if (d <= 1e-32) {
    dst = 0;
  } else {
    int e;

    T const f = std::frexp(d, &e) * 255.9999 / d;
    dst[0] = static_cast<hi::ubyte>(src[0] * f);
    dst[1] = static_cast<hi::ubyte>(src[1] * f);
    dst[2] = static_cast<hi::ubyte>(src[2] * f);
    dst[3] = static_cast<hi::ubyte>(e + E_OFFSET);
  }
}

template <typename T>
inline void hdr_decode(hi::basic_vector4<hi::ubyte> const &src,
                       hi::basic_vector3<T> &dst) {
  if (src[3]) {
    T const f = std::ldexp(T(1), src[3] - (E_OFFSET + 8));
    dst[0] = (src[0] + T(0.5)) * f;
    dst[1] = (src[1] + T(0.5)) * f;
    dst[2] = (src[2] + T(0.5)) * f;
  } else {
    dst = 0;
  }
}

void hdr_write_line(std::vector<xyze_t> const &scanline, std::ostream &out) {
  int size = static_cast<int>(scanline.size());

  if ((size < MIN_SCANLINE_LENGTH) || (MAX_SCANLINE_LENGTH < size)) {
    out.write(reinterpret_cast<char const *>(&scanline[0]),
              sizeof(xyze_t) * size);
    return;
  }

  // write magic header
  out.put(2);
  out.put(2);

  out.put((size >> 8) & 0xff);
  out.put((size >> 0) & 0xff);

  for (int c = 0; c < 4; ++c) {
    for (int j = 0, count = 0; j < size; j += count) {
      // find next run
      int begin;  // run の開始位置
      for (begin = j; begin < size; begin += count) {
        for (count = 1; (count < 127) && (begin + count < size) &&
                        (scanline[begin + count][c] == scanline[begin][c]);
             ++count) {
        }

        // long enough
        if (count >= MIN_RUNLENGTH) {
          break;
        }
      }

      if (((begin - j) > 1) && ((begin - j) < MIN_RUNLENGTH)) {
        int code = j + 1;
        while (scanline[code++][c] == scanline[j][c]) {
          // short run
          if (code == begin) {
            out.put(begin - j + 128);
            out.put(scanline[j][c]);
            j = begin;
            break;
          }
        }
      }

      // write out non-run
      while (j < begin) {
        int code = begin - j;
        if (code > 128) {
          code = 128;
        }
        out.put(code);
        while (code--) {
          out.put(scanline[j++][c]);
        }
      }

      // write out run
      if (count >= MIN_RUNLENGTH) {
        out.put(count + 128);
        out.put(scanline[begin][c]);
      } else {
        count = 0;
      }
    }
  }
}
}

namespace xyz {
bool ImageFilm::SaveAsHdr(std::tchar_t const *filename,
                          float_t const mpp) const {
  std::ofstream out(filename, std::ios::out | std::ios::binary);

  if (!out) {
    return false;
  }
  out.imbue(std::locale("C"));

  out << "#?RADIANCE\n"
         "# Made with Mikado Renderer\n"
#ifndef CONFIG_HDR_STRICT
         "FORMAT=32-bit_rle_rgbe\n"
#else
         "FORMAT=32-bit_rle_xyze\n"
#endif
         "EXPOSURE=1.0000000000000\n"
         "\n"
#ifndef CONFIG_HDR_STRICT
         "-Y " << height() << " +X " << width() << "\n";
#else
         "+Y " << height() << " +X " << width() << "\n";
#endif

  float_t const s = hi::rcp(mpp);

  std::size_t const width = this->width();
  std::size_t const height = this->height();

  ImageFilm::color_t xyz;
#ifndef CONFIG_HDR_STRICT
  ImageFilm::color_t rgb;
#endif
  std::vector<xyze_t> scanline(width);

  std::size_t ignore_pixel = 0;
#ifndef CONFIG_HDR_STRICT
  for (std::size_t p = (height - 1) * width, y = 0; y < height;
       ++y, p -= width * 2)
#else
  for (std::size_t p = 0, y = 0; y < height; ++y)
#endif
  {
    for (std::size_t x = 0; x < width; ++x, ++p) {
      xyz = m_[p];
      xyz *= s;

#ifndef CONFIG_HDR_STRICT
      hi::CCIR601_1_XYZ2RGB(xyz, rgb);

      if (!::_finite(rgb[0]) || !::_finite(rgb[1]) || !::_finite(rgb[2])) {
        rgb = 0;  // 無限大は0に!!
        ++ignore_pixel;
      }
#if 0
        for (int c = 0; c < 3; ++c)
        {
          if (rgb[c] < 0)
          {
            rgb[c] = 0;
          }
        }
#endif
      ::hdr_encode(rgb, scanline[x]);
#else
      ::hdr_encode(xyz, scanline[x]);
#endif
    }
    ::hdr_write_line(scanline, out);
  }

  return true;
}

}  // end of namespace xyz
