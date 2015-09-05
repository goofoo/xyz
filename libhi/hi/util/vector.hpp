#pragma once

#include <hi/lang/lang.hpp>

namespace hi {
template <typename T, typename A>
inline void trim(std::vector<T, A> &v) {
  std::vector<T, A>(v).swap(v);  // the swap trick
}
}
