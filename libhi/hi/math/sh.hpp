#include <hi/lang.hpp>

namespace hi {
template <typename T>
T SphericalHamonics(int const l, int const m, T const theta, T const phi);

#define MATERIALIZE(T)                                            \
  template <>                                                     \
  T SphericalHamonics<T>(int const, int const, T const, T const); \
  template <>                                                     \
  T SphericalHamonics<T>(int const, int const, T const, T const); \
  template <>                                                     \
  T SphericalHamonics<T>(int const, int const, T const, T const); \
  template <>                                                     \
  T SphericalHamonics<T>(int const, int const, T const, T const)
MATERIALIZE(float);
MATERIALIZE(double);
MATERIALIZE(long double);
#undef MATERIALIZE
}
// end of namespace hi
