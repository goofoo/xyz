#pragma once

#include <hi/lang/lang.hpp>

namespace hi {
inline bool starts_with(std::tstring const &s, std::tstring const &h) {
  if (s.length() < h.length()) {
    return false;
  }

  return s.substr(0, h.length()) == h;
}

inline bool ends_with(std::tstring const &s, std::tstring const &h) {
  if (s.length() < h.length()) {
    return false;
  }

  return s.substr(s.length() - h.length(), h.length()) == h;
}

inline void trim(std::tstring &s) {
  if (s.empty()) {
    return;
  }

  std::size_t const begin = s.find_first_not_of(_TEXT(" \t\r\n"));
  if (std::tstring::npos == begin) {
    std::tstring().swap(s);
  } else {
    std::size_t const end = s.find_last_not_of(_TEXT(" \t\r\n"));
    s.substr(begin, end - begin + 1).swap(s);
  }
}
}  // end of namespace hi
