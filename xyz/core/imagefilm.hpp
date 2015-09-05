#ifndef XYZ_FILM_HPP_
#define XYZ_FILM_HPP_

namespace xyz {
class ImageFilm {
 public:
  typedef pixel_descriptor_t::second_type color_t;

 public:
  inline ImageFilm() : width_(0), height_(0), m_() {}

  inline ImageFilm(__in std::size_t const width, __in std::size_t const height)
      : width_(width), height_(height), m_(width * height, color_t(0)) {}

 public:
  inline std::size_t width() const { return width_; }
  inline std::size_t height() const { return height_; }
  inline std::size_t count() const { return m_.size(); }

  inline ImageFilm::color_t &operator[](__in std::size_t const p) {
    return m_[p];
  }

  void Clear();
  void Resize(__in std::size_t const, __in std::size_t const);
  void Deposite(__in pixel_descriptor_t const &value,
                __in float_t const fWeightingFactor = 1);
  void Deposite(__in std::vector<pixel_descriptor_t> const &values,
                __in float_t const fWeightingFactor = 1);

  bool Save(std::tchar_t const *filename, float_t const) const;
  bool SaveAsPfmRgb(std::tchar_t const *filename, float_t const) const;

 private:
  bool SaveAsFloat(std::tchar_t const *, float_t const) const;
  bool SaveAsPfm(std::tchar_t const *, float_t const) const;
  bool SaveAsHdr(std::tchar_t const *, float_t const) const;

 private:
  std::size_t width_;
  std::size_t height_;
  std::vector<color_t> m_;

  ImageFilm(ImageFilm const &);
  ImageFilm &operator=(ImageFilm const &);

#ifdef CONFIG_MULTI_THREADING
#ifndef CONCURRENT_LOCK
  hi::monitor_object monitor_;
#endif
#endif
};

}  // end of namespace xyz

#endif XYZ_FILM_HPP_
