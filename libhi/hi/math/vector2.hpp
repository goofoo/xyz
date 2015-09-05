#pragma once

#include <hi/lang.hpp>
#include <hi/math/promote_traits.hpp>

namespace hi {
template <typename T>
class basic_vector2 sealed {
 public:
  typedef T value_type;
  static std::size_t const N = 2;

 public:
  // interface for STL
  static std::size_t size() { return N; }
  static std::size_t max_size() { return size(); }

 public:
  inline basic_vector2() {}

  template <typename S>
  inline basic_vector2(basic_vector2<S> const &v) {
    (*this) = v;
  }

  inline explicit basic_vector2(T const &s) { (*this) = s; }

  inline basic_vector2(T const &x, T const &y) {
    (*this)[0] = x;
    (*this)[1] = y;
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
  inline basic_vector2<T> &operator Op(basic_vector2<S> const &a) { \
    (*this)[0] Op static_cast<T>(a[0]);                             \
    (*this)[1] Op static_cast<T>(a[1]);                             \
    return *this;                                                   \
  }                                                                 \
  template <typename S>                                             \
  inline basic_vector2<T> &operator Op(S const a) {                 \
    (*this)[0] Op static_cast<T>(a);                                \
    (*this)[1] Op static_cast<T>(a);                                \
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
  inline basic_vector2<typename promote_traits<T, S>::type> operator Op( \
      T const a, basic_vector2<S> const &b) {                            \
    typedef typename promote_traits<T, S>::type R;                       \
    return basic_vector2<R>(a Op b[0], a Op b[1]);                       \
  }                                                                      \
  template <typename T, typename S>                                      \
  inline basic_vector2<typename promote_traits<T, S>::type> operator Op( \
      basic_vector2<T> const &a, S const b) {                            \
    typedef typename promote_traits<T, S>::type R;                       \
    return basic_vector2<R>(a[0] Op b, a[1] Op b);                       \
  }                                                                      \
  template <typename T, typename S>                                      \
  inline basic_vector2<typename promote_traits<T, S>::type> operator Op( \
      basic_vector2<T> const &a, basic_vector2<S> const &b) {            \
    typedef typename promote_traits<T, S>::type R;                       \
    return basic_vector2<R>(a[0] Op b[0], a[1] Op b[1]);                 \
  }
VECTOR_IMPLEMENTS(+)
VECTOR_IMPLEMENTS(-)
VECTOR_IMPLEMENTS(*)
VECTOR_IMPLEMENTS(/ )
#undef VECTOR_IMPLEMENTS

#define VECTOR_IMPLEMENTS(Op)                                      \
  template <typename T>                                            \
  inline basic_vector2<T> operator Op(basic_vector2<T> const &a) { \
    return basic_vector2<T>(Op a[0], Op a[1]);                     \
  }
VECTOR_IMPLEMENTS(+)
VECTOR_IMPLEMENTS(-)
#undef VECTOR_IMPLEMENTS

#define VECTOR_IMPLEMENTS(Op)                          \
  template <typename T>                                \
  inline bool operator Op(basic_vector2<T> const &a,   \
                          basic_vector2<T> const &b) { \
    return (a[0] Op b[0]) && (a[1] Op b[1]);           \
  }
VECTOR_IMPLEMENTS(< )
VECTOR_IMPLEMENTS(> )
VECTOR_IMPLEMENTS(<= )
VECTOR_IMPLEMENTS(>= )
VECTOR_IMPLEMENTS(== )
#undef VECTOR_IMPLEMENTS

template <typename T, typename S>
inline typename promote_traits<T, S>::type dot(basic_vector2<T> const &a,
                                               basic_vector2<S> const &b) {
  return a[0] * b[0] + a[1] * b[1];
}

template <typename T, typename S>
inline typename promote_traits<T, S>::type cross(basic_vector2<T> const &a,
                                                 basic_vector2<S> const &b) {
  return a[0] * b[1] - a[1] * b[0];
}

template <typename T>
inline T l1_norm(basic_vector2<T> const &a) {
  return std::abs(a[0]) + std::abs(a[1]);
}

template <typename T>
inline T max_norm(basic_vector2<T> const &a) {
  return hi::max(std::abs(a[0]), std::abs(a[1]));
}

template <typename T>
inline T sum(basic_vector2<T> const &a) {
  return a[0] + a[1];
}

template <typename T>
inline T prod(basic_vector2<T> const &a) {
  return a[0] * a[1];
}

template <typename T>
inline std::size_t min_index_of(basic_vector2<T> const &a) {
  return (a[0] <= a[1]) ? 0 : 1;
}

template <typename T>
inline std::size_t max_index_of(basic_vector2<T> const &a) {
  return (a[0] >= a[1]) ? 0 : 1;
}
}

// ƒ‰ƒCƒuƒ‰ƒŠ‚Ì min,max ‚Ì“ÁŽê‰»
namespace hi {
template <typename T, typename S>
inline basic_vector2<typename promote_traits<T, S>::type> min(
    basic_vector2<T> const &a, basic_vector2<T> const &b) {
  typedef typename promote_traits<T, S>::type R;
  return basic_vector2<R>(hi::min(a[0], b[0]), hi::min(a[1], b[1]));
}

template <typename T, typename S>
inline basic_vector2<typename promote_traits<T, S>::type> max(
    basic_vector2<T> const &a, basic_vector2<T> const &b) {
  typedef typename promote_traits<T, S>::type R;
  return basic_vector2<R>(hi::max(a[0], b[0]), hi::max(a[1], b[1]));
}

template <typename T, typename S>
inline void min(basic_vector2<T> &a, basic_vector2<S> const &b) {
  hi::min(a[0], b[0]);
  hi::min(a[1], b[1]);
}

template <typename T, typename S>
inline void max(basic_vector2<T> &a, basic_vector2<S> const &b) {
  hi::max(a[0], b[0]);
  hi::max(a[1], b[1]);
}

}  // end of namespace hi

namespace hi {
template <typename T, typename S>
inline typename promote_traits<T, S>::type Det(basic_vector2<T> const &a,
                                               basic_vector2<S> const &b) {
  return hi::cross(a, b);
}

template <typename T>
T SolveLinearSystem(basic_vector2<T> const &a, basic_vector2<T> const &b,
                    basic_vector2<T> &y) {
  T const det = hi::Det(a, b);

  if (T(0) != det) {
    T const inv_det = hi::rcp(det);
    T const tmp = inv_det * hi::Det(a, y);
    y[0] = inv_det * hi::Det(y, b);
    y[1] = tmp;
  }

  return det;
}

template <typename T>
std::tistream &operator>>(std::tistream &in, basic_vector2<T> &v) {
  return in >> v[0] >> v[1];
}
}
