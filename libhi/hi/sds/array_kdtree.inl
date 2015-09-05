#pragma once

namespace hi
{
  template <typename T>
  array_kdtree<T>::array_kdtree()
  {
  }

  template <typename T>
  void array_kdtree<T>::clear(std::size_t const max_size)
  {
#if 0
    data_.swap(std::vector<T>());
    if (max_size > 0)
    {
      data_.reserve(max_size + 1);
      data_.push_back(T());
    }
#else
    data_.clear();
    data_.reserve(max_size + 1);
    data_.push_back(T());
#endif

    min_ =  std::numeric_limits<typename T::vector_type::value_type>::infinity();
    max_ = -std::numeric_limits<typename T::vector_type::value_type>::infinity();
  }

  template <typename T>
  void array_kdtree<T>::push_back(T const & value)
  {
    data_.push_back(value);
    hi::min(min_, value.val());
    hi::max(max_, value.val());
  }

  template <typename T>
  void array_kdtree<T>::build()
  {
    if (data_.size() <= 2)
    {
      return;
    }

    // ary を通りがけ順(in-order traversal)で整列させると同時に、
    // これを利用して i が親ノードであれば左右のノードのインデックスに
    // それぞれ 2*i, 2*i+1 でアクセスできるように ary を整列させる。

    std::vector<std::size_t> ary(data_.size()); // 
    std::vector<std::size_t> iot(data_.size()); // in-order traversal

    for (std::size_t i = 0, size = data_.size(); i < size; ++i)
    {
      ary[i] = iot[i] = i;
    }

    build(ary, iot, 1, 1, data_.size() - 1);

    // メモリコンパクションとスワップ
    std::vector<T> temp(data_.size());
    for (std::size_t i = 0, size = data_.size(); i < size; ++i)
    {
      temp[i] = data_[ary[i]];
    }
    temp.swap(data_);

  //std::fprintf(stderr, "add photons: %d\n", data_.size() - 1);
  }

  template <typename T>
  void array_kdtree<T>::build(
    std::vector<std::size_t> & ary, std::vector<std::size_t> & iot,
    std::size_t i, std::size_t start, std::size_t end)
  {
    std::size_t median = find_median_index(start, end); // 要素中でちょうど中央に来るインデックス
    std::size_t axis = hi::max_index_of(max_ - min_); // 分割軸を決定する

  //std::cerr << "begin split" << std::endl;
    split(iot, start, end, median, axis); // 中央値で左右に振り分ける
  //std::cerr << "end split" << std::endl;

    ary[i] = iot[median];
    data_[ ary[i] ].key() = axis;

    if (median > start) // 左
    {
      if (start < median - 1)
      {
        float max = max_[axis];
        max_[axis] = data_[ ary[i] ].val()[axis];
        build(ary, iot, i + i, start, median - 1);
        max_[axis] = max;
      }
      else
      {
        // 残り1つ
        ary[i+i] = iot[start];
      }
    }

    if (median < end) // 右
    {
      if (median + 1 < end)
      {
        float min = min_[axis];
        min_[axis] = data_[ ary[i] ].val()[axis];
        build(ary, iot, i + i + 1, median + 1, end);
        min_[axis] = min;
      }
      else
      {
        // 残り1つ
        ary[i+i+1] = iot[end];
      }
    }
  }

  template <typename T>
  void array_kdtree<T>::split(
    std::vector<std::size_t> & iot,
    std::size_t start, std::size_t end,
    std::size_t median, std::size_t axis)
  {
    std::size_t left = start;
    std::size_t right = end;

  //std::cerr << "split [" << iot.size() << "] (" << start << ", " << end << ")" << std::endl;

    while (right > left)
    {
      //std::cerr << "key [" << data_.size() << "] (" << iot[right] << ")" << std::endl;

      // 右端をキーに設定
      T::value_type key = data_[ iot[right] ].val()[axis];

      std::size_t i = left - 1;
      std::size_t j = right;

      // キーより小さいものは左へ、大きいものは右へ移動
      for (; ;)
      {
        while (data_[ iot[++i] ].val()[axis] < key)
        {
          if (i >= right)
          {
            break;
          }
        }
        while (data_[ iot[--j] ].val()[axis] > key)
        {
          if (j <= left)
          {
            break;
          }
        }
        if (i >= j)
        {
          break;
        }

        std::swap(iot[i], iot[j]);
      }
 
      // 左から検索した位置(i)とキー(right)を入れ替える
      // (iを中心に大きいものと小さいものが左右に分かれた)
      std::swap(iot[i], iot[right]);

      // iが分割点より右にあったら右の位置を更新
      if (i >= median)
      {
        right = i - 1;
      }

      // 同上
      if (i <= median)
      {
        left = i + 1;
      }
    }
  }

  // [start,end] の要素から左寄り完全平衡 kd-tree のルートとなるインデックスを検索
  template <typename T>
  std::size_t array_kdtree<T>::find_median_index(std::size_t start, std::size_t end)
  {
    std::size_t size = end - start + 1;
    std::size_t median = 1;

    // 最左ノードインデックスを見つける
    while ((median << 2) <= size)
    {
      median += median;
    }

    if ((3 * median) <= size)
    {
      // 右部分木に左と同じ深さの葉があるので
      return (median + median) + start - 1;
    }
    else
    {
      // ここで median は右部分木の要素数である
      // 中間点 = 要素数 - 右部分木の要素数
      return end - median + 1;
    }
  }

  template <typename T>
  template <typename Fn>
  void array_kdtree<T>::kNN_query(vector_type const & p, std::size_t const i, Fn & fn) const
  {
    vector_type const d = data_[i].val() - p;

    if (i + i + 1 < data_.size() - 1)
    {
      typename vector_type::value_type const dist = d[ data_[i].key() ]; //[TODO]`hi::real_t'に依存しないように書き直す

      if (dist < 0)
      {
        if ((dist * dist) < fn.max_search_range())
        {
          kNN_query<Fn>(p, i + i, fn); // 左の子
        }

        kNN_query<Fn>(p, i + i + 1, fn); // 右の子
      }
      else
      {
        kNN_query<Fn>(p, i + i, fn); // 左の子

        if ((dist * dist) < fn.max_search_range())
        {
          kNN_query<Fn>(p, i + i + 1, fn); // 右の子
        }
      }
    }

    if (hi::length_squared(d) < fn.max_search_range())
    {
      fn.push_back(data_[i], d);
    }
  }

} // end of namespace hi
