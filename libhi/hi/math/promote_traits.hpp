#pragma once

namespace hi {
/// <summary>
/// type = C ? S : T;
/// </summary>
template <bool C, typename S, typename T>
struct if_ {
  typedef S type;
};

template <typename S, typename T>
struct if_<false, S, T> {
  typedef T type;
};

/// <summary>
/// T ‚Æ S ‚Ì‰‰ŽZŒã‚ÌŒ^‚ª type
/// </summary>
template <typename S, typename T>
struct promote_traits;

#define PROMOTE_TRAITS_IMPLEMENTS(T) \
  template <>                        \
  struct promote_traits<T, T> {      \
    typedef T type;                  \
  };

PROMOTE_TRAITS_IMPLEMENTS(int);
PROMOTE_TRAITS_IMPLEMENTS(float);
PROMOTE_TRAITS_IMPLEMENTS(double);
PROMOTE_TRAITS_IMPLEMENTS(long double);

#undef PROMOTE_TRAITS_IMPLEMENTS

#define PROMOTE_TRAITS_IMPLEMENTS(S, T) \
  template <>                           \
  struct promote_traits<S, T> {         \
    typedef S type;                     \
  };                                    \
  template <>                           \
  struct promote_traits<T, S> {         \
    typedef S type;                     \
  };

PROMOTE_TRAITS_IMPLEMENTS(float, int);

PROMOTE_TRAITS_IMPLEMENTS(double, int);
PROMOTE_TRAITS_IMPLEMENTS(double, float);

PROMOTE_TRAITS_IMPLEMENTS(long double, int);
PROMOTE_TRAITS_IMPLEMENTS(long double, float);
PROMOTE_TRAITS_IMPLEMENTS(long double, double);

#undef PROMOTE_TRAITS_IMPLEMENTS

}  // end of namespace hi
