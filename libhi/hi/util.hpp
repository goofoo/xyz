#pragma once

#include <hi/lang.hpp>
#include <hi/util/mt19937.hpp>
#include <hi/util/vector.hpp>
#include <hi/util/string.hpp>

namespace hi {
template <typename T>
struct delete_ptr {
  void operator()(T *p) const { delete p; }

  void operator()(T const *p) const { delete p; }
};

template <typename T>
struct delete_ary {
  void operator()(T *p) const { delete[] p; }

  void operator()(T const *p) const { delete[] p; }
};

template <typename T>
struct single_stack {
  single_stack(T &ref, T val) : ref_(ref), val_(val) { std::swap(ref, val); }

  ~single_stack() { ref_ = val_; }

 private:
  T &ref_;
  T val_;

 private:
  single_stack(single_stack const &);
  single_stack &operator=(single_stack const &);
};
}
