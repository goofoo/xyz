#pragma once

#include <hi/lang.hpp>
#include <hi/math/promote_traits.hpp>

namespace hi {
template <typename T>
class basic_vector3 sealed {
 public:
  typedef T value_type;
  static std::size_t const N = 3;

 public:
  // interface for STL
  static std::size_t size() { return N; }
  static std::size_t max_size() { return size(); }

 public:
  inline basic_vector3() {}

  template <typename S>
  inline basic_vector3(basic_vector3<S> const &v) {
    (*this) = v;
  }

  inline explicit basic_vector3(T const &s) { (*this) = s; }

  inline basic_vector3(T const &x, T const &y, T const &z) {
    (*this)[0] = x;
    (*this)[1] = y;
    (*this)[2] = z;
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

#define VECTOR_IMPLEMENTS(Op)                                         \
  template <typename S>                                               \
  inline basic_vector3<T> &operator Op(basic_vector3<S> const &rhs) { \
    (*this)[0] Op static_cast<T>(rhs[0]);                             \
    (*this)[1] Op static_cast<T>(rhs[1]);                             \
    (*this)[2] Op static_cast<T>(rhs[2]);                             \
    return *this;                                                     \
  }                                                                   \
  template <typename S>                                               \
  inline basic_vector3<T> &operator Op(S const rhs) {                 \
    (*this)[0] Op static_cast<T>(rhs);                                \
    (*this)[1] Op static_cast<T>(rhs);                                \
    (*this)[2] Op static_cast<T>(rhs);                                \
    return *this;                                                     \
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
  inline basic_vector3<typename promote_traits<T, S>::type> operator Op( \
      T const a, basic_vector3<S> const &b) {                            \
    typedef typename promote_traits<T, S>::type R;                       \
    return basic_vector3<R>(a Op b[0], a Op b[1], a Op b[2]);            \
  }                                                                      \
  template <typename T, typename S>                                      \
  inline basic_vector3<typename promote_traits<T, S>::type> operator Op( \
      basic_vector3<T> const &a, S const b) {                            \
    typedef typename promote_traits<T, S>::type R;                       \
    return basic_vector3<R>(a[0] Op b, a[1] Op b, a[2] Op b);            \
  }                                                                      \
  template <typename T, typename S>                                      \
  inline basic_vector3<typename promote_traits<T, S>::type> operator Op( \
      basic_vector3<T> const &a, basic_vector3<S> const &b) {            \
    typedef typename promote_traits<T, S>::type R;                       \
    return basic_vector3<R>(a[0] Op b[0], a[1] Op b[1], a[2] Op b[2]);   \
  }
VECTOR_IMPLEMENTS(+)
VECTOR_IMPLEMENTS(-)
VECTOR_IMPLEMENTS(*)
VECTOR_IMPLEMENTS(/ )
#undef VECTOR_IMPLEMENTS

#define VECTOR_IMPLEMENTS(Op)                                      \
  template <typename T>                                            \
  inline basic_vector3<T> operator Op(basic_vector3<T> const &a) { \
    return basic_vector3<T>(Op a[0], Op a[1], Op a[2]);            \
  }
VECTOR_IMPLEMENTS(+)
VECTOR_IMPLEMENTS(-)
#undef VECTOR_IMPLEMENTS

#define VECTOR_IMPLEMENTS(Op)                                  \
  template <typename T, typename S>                            \
  inline bool operator Op(basic_vector3<T> const &a,           \
                          basic_vector3<S> const &b) {         \
    return (a[0] Op b[0]) && (a[1] Op b[1]) && (a[2] Op b[2]); \
  }
VECTOR_IMPLEMENTS(< )
VECTOR_IMPLEMENTS(> )
VECTOR_IMPLEMENTS(<= )
VECTOR_IMPLEMENTS(>= )
VECTOR_IMPLEMENTS(== )
#undef VECTOR_IMPLEMENTS

template <typename T, typename S>
inline typename promote_traits<T, S>::type dot(basic_vector3<T> const &a,
                                               basic_vector3<S> const &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <typename T, typename S>
inline basic_vector3<typename promote_traits<T, S>::type> cross(
    basic_vector3<T> const &a, basic_vector3<S> const &b) {
  typedef typename promote_traits<T, S>::type R;
  return basic_vector3<R>(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                          a[0] * b[1] - a[1] * b[0]);
}

template <typename T>
inline T l1_norm(basic_vector3<T> const &a) {
  return std::abs(a[0]) + std::abs(a[1]) + std::abs(a[2]);
}

template <typename T>
inline T max_norm(basic_vector3<T> const &a) {
  return hi::max(hi::max(std::abs(a[0]), std::abs(a[1])), std::abs(a[2]));
}

template <typename T>
inline T sum(basic_vector3<T> const &a) {
  return a[0] + a[1] + a[2];
}

template <typename T>
inline T prod(basic_vector3<T> const &a) {
  return a[0] * a[1] * a[2];
}

template <typename T>
inline std::size_t min_index_of(basic_vector3<T> const &a) {
  return ((a[0] <= a[1]) && (a[0] <= a[2])) ? 0 : (a[1] <= a[2]) ? 1 : 2;
}

template <typename T>
inline std::size_t max_index_of(basic_vector3<T> const &a) {
  return ((a[0] >= a[1]) && (a[0] >= a[2])) ? 0 : (a[1] >= a[2]) ? 1 : 2;
}

template <typename T>
inline T absolute_error(basic_vector3<T> const &a, basic_vector3<T> const &b) {
  return hi::length(a - b);
}

template <typename T>
inline T relative_error(basic_vector3<T> const &a, basic_vector3<T> const &b) {
  T const d = hi::length(b);
  return (d > T(0)) ? absolute_error(a, b) / d : T(0);
}
}

// ƒ‰ƒCƒuƒ‰ƒŠ‚Ì min,max ‚Ì“ÁŽê‰»
namespace hi {
template <typename T, typename S>
inline basic_vector3<typename promote_traits<T, S>::type> min(
    basic_vector3<T> const &a, basic_vector3<S> const &b) {
  typedef typename promote_traits<T, S>::type R;
  return basic_vector3<R>(hi::min(a[0], b[0]), hi::min(a[1], b[1]),
                          hi::min(a[2], b[2]));
}

template <typename T, typename S>
inline basic_vector3<T> max(basic_vector3<T> const &a,
                            basic_vector3<S> const &b) {
  typedef typename promote_traits<T, S>::type R;
  return basic_vector3<R>(hi::max(a[0], b[0]), hi::max(a[1], b[1]),
                          hi::max(a[2], b[2]));
}

template <typename T, typename S>
inline void min(basic_vector3<T> &a, basic_vector3<S> const &b) {
  hi::min(a[0], b[0]);
  hi::min(a[1], b[1]);
  hi::min(a[2], b[2]);
}

template <typename T, typename S>
inline void max(basic_vector3<T> &a, basic_vector3<S> const &b) {
  hi::max(a[0], b[0]);
  hi::max(a[1], b[1]);
  hi::max(a[2], b[2]);
}

}  // end of namespace hi

namespace hi {
template <typename T>
inline T Luminance(basic_vector3<T> const &a) {
  return T(0.299) * a[0] + T(0.587) * a[1] + T(0.114) * a[2];
}

template <typename T>
inline T Det(basic_vector3<T> const &a, basic_vector3<T> const &b,
             basic_vector3<T> const &c) {
  return hi::dot(hi::cross(a, b), c);
}

template <typename T>
T SolveLinearSystem(basic_vector3<T> const &a, basic_vector3<T> const &b,
                    basic_vector3<T> const &c, basic_vector3<T> &y) {
  basic_vector3<T> const bxc = hi::cross(b, c);

  T const det = hi::dot(bxc, a);

  if (T(0) != det) {
    basic_vector3<T> const axy = hi::cross(a, y);

    T const inv_det = hi::rcp(det);
    y[0] = inv_det * hi::dot(bxc, y);
    y[1] = inv_det * hi::dot(axy, c);
    y[2] = -inv_det * hi::dot(axy, b);
  }

  return det;
}

template <typename T>
std::tistream &operator>>(std::tistream &in, basic_vector3<T> &v) {
  return in >> v[0] >> v[1] >> v[2];
}

}  // end of namespace hi
