#pragma once

namespace hi {
/*
  格納データに要求されるポリシー

  class Data
  {
  public:
    typedef hi::Vector3 vector_type;
    typedef vector_type::value_type value_type;

    vector_type const & val() const { return val_; }
    std::size_t       & key()       { return key_; }
    std::size_t const & key() const { return key_; }

  private:
    vector_type val_;
    std::size_t key_;
  };


  クエリ関数に要求されるポリシー

  class Query
  {
  public:
    Data::value_type max_search_range() const;

    void push_back(Data const & data, Data::vector_type const & d);

  private:
  };
*/

/// <summary>
/// Array KDTree (静的な多次元点データへのkNNクエリを実現するクラス)
/// </summary>
template <typename T>
class array_kdtree {
 public:
  typedef T value_type;
  typedef typename T::vector_type vector_type;

 public:  // some interfaces for STL
  typedef typename std::vector<T>::const_iterator const_iterator;

  const_iterator begin() const { return data_.begin() + 1; }
  const_iterator end() const { return data_.end(); }

 public:
  array_kdtree();

  void clear(std::size_t const max_size = 0);
  void push_back(T const &value);

  void build();

  // kNN query
  template <typename Fn>
  inline void kNN_query(vector_type const &p, Fn &fn) const {
    fn.begin_query();
    kNN_query(p, 1, fn);
    fn.end_query();
  }

 private:
  void build(std::vector<std::size_t> &ary, std::vector<std::size_t> &iot,
             std::size_t i, std::size_t start, std::size_t end);

  void split(std::vector<std::size_t> &indices, std::size_t start,
             std::size_t end, std::size_t median, std::size_t axis);

  static std::size_t find_median_index(std::size_t start, std::size_t end);

  template <typename Fn>
  void kNN_query(vector_type const &p, std::size_t const i, Fn &fn) const;

 private:
  array_kdtree(array_kdtree const &);
  array_kdtree &operator=(array_kdtree const &);

 private:
  // 格納データ
  std::vector<T> data_;
  vector_type min_;  // Bounding Box
  vector_type max_;  //
};

}  // end of namespace hi

#include "array_kdtree.inl"

/*
NOTE: Photon mapping tricks
  A practical guide to global illumination using Photon mapping,
  SIGGRAPH 2001 Cource 38, Chapter 5, August 14, 2001.

  (1) Maximum search radius for surfaces
    Lr = \sum_{i=1}^{n} f_r(x,wi,wo) P[i] / A
      <= n * P_max / ((pi)^2 * r_max^2)
           f_r = 1 / (pi)
           \sum_{i=1}^{n} P[i] = n * P_max
           A = (pi) * r_max^2
       = Lt

      Lt = (threshold)
      n = (number of nearest photons)
      P_max = (maximum power of Photon in Photon-map)
      r_max = (maximum serch radius)

    r_max^2 = n * P_max / (Lt*(pi)^2)

  (2) Maximum search radius for volumes
    Lr = \sum_{i=1}^{n} f(x,wi,w) P[i] / (rho(x) * V)
      <= (3/16) * n * P_max / ((pi)^2 * rho(x) * r_max^3)
           f = 1 / (4 * (pi))
           \sum_{i=1}^{n} P[i] = P_max
           V = (4/3) * (pi) * r^3
       = Lt

      Lt = (threshold)
      n = (number of nearest photons)
      P_max = (maximum power of Photon in Photon-map)
      r_max = (maximum serch radius)

    r_max^3 = ((3/16) * n * P_max) / ((pi)^2 * rho(x) * Lt)

  Lt = 0.05の場合
    (1) r_max^2 (approx) 2 * n * P_max
    (2) r_max^3 (approx) 2 * (3/16) * n * P_max / rho
*/
