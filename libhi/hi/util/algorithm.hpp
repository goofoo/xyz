#include <hi/lang.hpp>

namespace hi {
// 単純選択ソート
template <typename Iterator>
void selection_sort(Iterator first, Iterator last) {
  for (; first != last; ++first) {
    std::iter_swap(first, std::min_element(first, last));
  }
}

// マージソート
template <typename Iterator>
void marge_sort(Iterator first, Iterator last) {
  typename std::iterator_traits<Iterator>::difference_type const threshold = 36;
  std::iterator_traits<Iterator>::difference_type const size =
      std::distance(first, last);

  if (size < threshold) {
    selection_sort(first, last);
  } else {
    Iterator middle = first;
    std::advance(middle, size / 2);

    marge_sort(first, middle);
    marge_sort(middle, last);

    std::inplace_merge(first, middle, last);
  }
}
}
