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

    // ary ��ʂ肪����(in-order traversal)�Ő��񂳂���Ɠ����ɁA
    // ����𗘗p���� i ���e�m�[�h�ł���΍��E�̃m�[�h�̃C���f�b�N�X��
    // ���ꂼ�� 2*i, 2*i+1 �ŃA�N�Z�X�ł���悤�� ary �𐮗񂳂���B

    std::vector<std::size_t> ary(data_.size()); // 
    std::vector<std::size_t> iot(data_.size()); // in-order traversal

    for (std::size_t i = 0, size = data_.size(); i < size; ++i)
    {
      ary[i] = iot[i] = i;
    }

    build(ary, iot, 1, 1, data_.size() - 1);

    // �������R���p�N�V�����ƃX���b�v
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
    std::size_t median = find_median_index(start, end); // �v�f���ł��傤�ǒ����ɗ���C���f�b�N�X
    std::size_t axis = hi::max_index_of(max_ - min_); // �����������肷��

  //std::cerr << "begin split" << std::endl;
    split(iot, start, end, median, axis); // �����l�ō��E�ɐU�蕪����
  //std::cerr << "end split" << std::endl;

    ary[i] = iot[median];
    data_[ ary[i] ].key() = axis;

    if (median > start) // ��
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
        // �c��1��
        ary[i+i] = iot[start];
      }
    }

    if (median < end) // �E
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
        // �c��1��
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

      // �E�[���L�[�ɐݒ�
      T::value_type key = data_[ iot[right] ].val()[axis];

      std::size_t i = left - 1;
      std::size_t j = right;

      // �L�[��菬�������͍̂��ցA�傫�����͉̂E�ֈړ�
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
 
      // �����猟�������ʒu(i)�ƃL�[(right)�����ւ���
      // (i�𒆐S�ɑ傫�����̂Ə��������̂����E�ɕ����ꂽ)
      std::swap(iot[i], iot[right]);

      // i�������_���E�ɂ�������E�̈ʒu���X�V
      if (i >= median)
      {
        right = i - 1;
      }

      // ����
      if (i <= median)
      {
        left = i + 1;
      }
    }
  }

  // [start,end] �̗v�f���獶��芮�S���t kd-tree �̃��[�g�ƂȂ�C���f�b�N�X������
  template <typename T>
  std::size_t array_kdtree<T>::find_median_index(std::size_t start, std::size_t end)
  {
    std::size_t size = end - start + 1;
    std::size_t median = 1;

    // �ō��m�[�h�C���f�b�N�X��������
    while ((median << 2) <= size)
    {
      median += median;
    }

    if ((3 * median) <= size)
    {
      // �E�����؂ɍ��Ɠ����[���̗t������̂�
      return (median + median) + start - 1;
    }
    else
    {
      // ������ median �͉E�����؂̗v�f���ł���
      // ���ԓ_ = �v�f�� - �E�����؂̗v�f��
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
      typename vector_type::value_type const dist = d[ data_[i].key() ]; //[TODO]`hi::real_t'�Ɉˑ����Ȃ��悤�ɏ�������

      if (dist < 0)
      {
        if ((dist * dist) < fn.max_search_range())
        {
          kNN_query<Fn>(p, i + i, fn); // ���̎q
        }

        kNN_query<Fn>(p, i + i + 1, fn); // �E�̎q
      }
      else
      {
        kNN_query<Fn>(p, i + i, fn); // ���̎q

        if ((dist * dist) < fn.max_search_range())
        {
          kNN_query<Fn>(p, i + i + 1, fn); // �E�̎q
        }
      }
    }

    if (hi::length_squared(d) < fn.max_search_range())
    {
      fn.push_back(data_[i], d);
    }
  }

} // end of namespace hi
