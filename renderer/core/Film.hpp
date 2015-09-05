#ifndef __TGIR_CORE_FILM_HPP__
#define __TGIR_CORE_FILM_HPP__

namespace tgir {
class Path;

typedef std::pair<std::size_t, tgir::SpectrumVector> PixelDescriptor;

class Film {
 public:
  typedef tgir::PixelDescriptor::second_type value_type;

 public:
  inline Film() : width_(0), height_(0), x_() {}

  inline Film(std::size_t const width, std::size_t const height)
      : width_(width), height_(height), x_(width * height, value_type(0)) {}

 public:
  //__declspec(property(get=GetWidth)) std::size_t Width;
  //__declspec(property(get=GetHeight)) std::size_t Height;
  //__declspec(property(get=GetCount)) std::size_t Count;

  inline std::size_t GetWidth() const { return width_; }
  inline std::size_t GetHeight() const { return height_; }
  inline std::size_t GetCount() const { return x_.size(); }

  inline Film::value_type &operator[](__in std::size_t const p) {
    return x_[p];
  }

  void Clear();
  void Resize(__in std::size_t const, __in std::size_t const);
  void Deposite(__in tgir::PixelDescriptor const &value,
                __in tgir::Real const &w = 1);
  void Deposite(__in std::vector<tgir::PixelDescriptor> const &values,
                __in tgir::Real const &w = 1);

  bool Save(std::tchar_t const *filename, tgir::Real const) const;
  bool SaveAsPfmRgb(std::tchar_t const *filename, tgir::Real const) const;

 private:
  bool SaveAsFloat(std::tchar_t const *, tgir::Real const) const;
  bool SaveAsPfm(std::tchar_t const *, tgir::Real const) const;
  bool SaveAsHdr(std::tchar_t const *, tgir::Real const) const;

 private:
  std::size_t width_;
  std::size_t height_;
  std::vector<value_type> x_;

  Film(Film const &);
  Film &operator=(Film const &);

#ifdef CONFIG_MULTI_THREADING
#ifndef CONCURRENT_LOCK
  hi::monitor_object monitor_;
#endif
#endif
};

}  // end of namespace tgir

#endif __TGIR_CORE_FILM_HPP__
