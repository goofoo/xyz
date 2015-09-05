#pragma once

#include <hi/math/promote_traits.hpp>

namespace hi {
template <typename T>
class basic_vector4 sealed {
 public:
  typedef T value_type;
  static std::size_t const N = 4;

 public:
  // interface for STL
  static std::size_t size() { return N; }
  static std::size_t max_size() { return size(); }

 public:
  inline basic_vector4() {}

  template <typename S>
  inline basic_vector4(basic_vector4<S> const &v) {
    (*this) = v;
  }

  inline explicit basic_vector4(T const &s) { (*this) = s; }

  inline basic_vector4(T const &x, T const &y, T const &z, T const &w) {
    (*this)[0] = x;
    (*this)[1] = y;
    (*this)[2] = z;
    (*this)[3] = w;
  }

  inline T const *data() const { return data_; }
  inline T *data() { return data_; }

  inline T const &operator[](std::size_t const &i) const {
    assert(i < size());
    return data_[i];
  }
  inline T &operator[](std::size_t const &i) {
    assert(i < size());
    return data_[i];
  }

#define VECTOR_IMPLEMENTS(Op)                                       \
  template <typename S>                                             \
  inline basic_vector4<T> &operator Op(basic_vector4<S> const &a) { \
    (*this)[0] Op static_cast<T>(a[0]);                             \
    (*this)[1] Op static_cast<T>(a[1]);                             \
    (*this)[2] Op static_cast<T>(a[2]);                             \
    (*this)[3] Op static_cast<T>(a[3]);                             \
    return *this;                                                   \
  }                                                                 \
  template <typename S>                                             \
  inline basic_vector4<T> &operator Op(S const a) {                 \
    (*this)[0] Op static_cast<T>(a);                                \
    (*this)[1] Op static_cast<T>(a);                                \
    (*this)[2] Op static_cast<T>(a);                                \
    (*this)[3] Op static_cast<T>(a);                                \
    return *this;                                                   \
  }
  VECTOR_IMPLEMENTS(= )
  VECTOR_IMPLEMENTS(+= )
  VECTOR_IMPLEMENTS(-= )
  VECTOR_IMPLEMENTS(*= )
  VECTOR_IMPLEMENTS(/= )
#undef VECTOR_IMPLEMENTS

 private:
  T data_[N];
};

#define VECTOR_IMPLEMENTS(Op)                                            \
  template <typename T, typename S>                                      \
  inline basic_vector4<typename promote_traits<T, S>::type> operator Op( \
      T const a, basic_vector4<S> const &b) {                            \
    typedef typename promote_traits<T, S>::type R;                       \
    return basic_vector4<R>(a Op b[0], a Op b[1], a Op b[2], a Op b[3]); \
  }                                                                      \
  template <typename T, typename S>                                      \
  inline basic_vector4<typename promote_traits<T, S>::type> operator Op( \
      basic_vector4<T> const &a, S const b) {                            \
    typedef typename promote_traits<T, S>::type R;                       \
    return basic_vector4<R>(a[0] Op b, a[1] Op b, a[2] Op b, a[3] Op b); \
  }                                                                      \
  template <typename T, typename S>                                      \
  inline basic_vector4<typename promote_traits<T, S>::type> operator Op( \
      basic_vector4<T> const &a, basic_vector4<S> const &b) {            \
    typedef typename promote_traits<T, S>::type R;                       \
    return basic_vector4<R>(a[0] Op b[0], a[1] Op b[1], a[2] Op b[2],    \
                            a[3] Op b[3]);                               \
  }
VECTOR_IMPLEMENTS(+)
VECTOR_IMPLEMENTS(-)
VECTOR_IMPLEMENTS(*)
VECTOR_IMPLEMENTS(/ )
#undef VECTOR_IMPLEMENTS

#define VECTOR_IMPLEMENTS(Op)                                      \
  template <typename T>                                            \
  inline basic_vector4<T> operator Op(basic_vector4<T> const &a) { \
    return basic_vector4<T>(Op a[0], Op a[1], Op a[2], Op a[3]);   \
  }
VECTOR_IMPLEMENTS(+)
VECTOR_IMPLEMENTS(-)
#undef VECTOR_IMPLEMENTS

#define VECTOR_IMPLEMENTS(Op)                                    \
  template <typename T>                                          \
  inline bool operator Op(basic_vector4<T> const &a,             \
                          basic_vector4<T> const &b) {           \
    return (a[0] Op b[0]) && (a[1] Op b[1]) && (a[2] Op b[2]) && \
           (a[3] Op b[3]);                                       \
  }
VECTOR_IMPLEMENTS(< )
VECTOR_IMPLEMENTS(> )
VECTOR_IMPLEMENTS(<= )
VECTOR_IMPLEMENTS(>= )
VECTOR_IMPLEMENTS(== )
#undef VECTOR_IMPLEMENTS

template <typename T, typename S>
inline typename promote_traits<T, S>::type dot(basic_vector4<T> const &a,
                                               basic_vector4<S> const &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

template <typename T, typename S, typename R>
inline typename promote_traits<typename promote_traits<T, S>, R>::type cross(
    basic_vector4<T> const &a, basic_vector4<S> const &b,
    basic_vector4<T> const &c) {
  typedef typename promote_traits<typename promote_traits<T, S>, R>::type Q;
  return basic_vector4<Q>(a[1] * b[2] * c[3] - a[3] * b[2] * c[1],
                          a[2] * b[3] * c[0] - a[0] * b[3] * c[2],
                          a[3] * b[0] * c[1] - a[1] * b[0] * c[3],
                          a[0] * b[1] * c[2] - a[2] * b[1] * c[0]);
}

template <typename T>
inline T l1_norm(basic_vector4<T> const &a) {
  return std::abs(a[0]) + std::abs(a[1]) + std::abs(a[2]) + std::abs(a[3]);
}

template <typename T>
inline T max_norm(basic_vector4<T> const &a) {
  return hi::max(hi::max(std::abs(a[0]), std::abs(a[1])),
                 hi::max(std::abs(a[2]), std::abs(a[3])));
}

template <typename T>
inline T sum(basic_vector4<T> const &a) {
  return a[0] + a[1] + a[2] + a[3];
}

template <typename T>
inline T prod(basic_vector4<T> const &a) {
  return a[0] * a[1] * a[2] * a[3];
}

template <typename T>
inline std::size_t min_index_of(basic_vector4<T> const &a) {
  return ((a[0] <= a[1]) && (a[0] <= a[2])) ? 0 : (a[1] <= a[2]) ? 1 : 2;
}

template <typename T>
inline std::size_t max_index_of(basic_vector4<T> const &a) {
  return ((a[0] >= a[1]) && (a[0] >= a[2])) ? 0 : (a[1] >= a[2]) ? 1 : 2;
}
}

// ƒ‰ƒCƒuƒ‰ƒŠ‚Ì min,max ‚Ì“ÁŽê‰»
namespace hi {
template <typename T, typename S>
inline basic_vector4<typename promote_traits<T, S>::type> min(
    basic_vector4<T> const &a, basic_vector4<S> const &b) {
  typedef typename promote_traits<T, S>::type R;
  return basic_vector4<R>(hi::min(a[0], b[0]), hi::min(a[1], b[1]),
                          hi::min(a[2], b[2]), hi::min(a[3], b[3]));
}

template <typename T>
inline basic_vector4<T> max(basic_vector4<T> const &a,
                            basic_vector4<T> const &b) {
  typedef typename promote_traits<T, S>::type R;
  return basic_vector4<R>(hi::max(a[0], b[0]), hi::max(a[1], b[1]),
                          hi::max(a[2], b[2]), hi::max(a[3], b[3]));
}

template <typename T, typename S>
inline void min(basic_vector4<T> &a, basic_vector4<S> const &b) {
  hi::min(a[0], b[0]);
  hi::min(a[1], b[1]);
  hi::min(a[2], b[2]);
  hi::min(a[3], b[3]);
}

template <typename T, typename S>
inline void max(basic_vector4<T> &a, basic_vector4<S> const &b) {
  hi::max(a[0], b[0]);
  hi::max(a[1], b[1]);
  hi::max(a[2], b[2]);
  hi::max(a[3], b[3]);
}

}  // end of namespace hi

namespace hi {
//[TODO]ŽÀ‘•‚·‚é
#if 0
  template <typename T>
  inline T Det(
    basic_vector4<T> const & a, basic_vector4<T> const & b,
    basic_vector4<T> const & c, basic_vector4<T> const & d)
  {
    return hi::dot(hi::cross(a, b), c);
  }

  template <typename T>
  T SolveLinearSystem(
    basic_vector4<T> const & a, basic_vector4<T> const & b,
    basic_vector4<T> const & c, basic_vector4<T> const & d, basic_vector4<T> & y)
  {
    basic_vector4<T> const bxc = hi::cross<T,T>(b, c);

    T const det = hi::dot(bxc, a);

    if (T(0) != det)
    {
      basic_vector4<T> const axy = hi::cross(a, y);

      T const inv_det = hi::rcp(det);
      y[0] =  inv_det * hi::dot(bxc, y);
      y[1] =  inv_det * hi::dot(axy, c);
      y[2] = -inv_det * hi::dot(axy, b);
    }

    return det;
  }
#endif

template <typename T>
std::tistream &operator>>(std::tistream &in, basic_vector4<T> &v) {
  return in >> v[0] >> v[1] >> v[2] >> v[3];
}

}  // end of namespace hi
