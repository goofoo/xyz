#ifndef XYZ_GEOM_INTERSECTION_HPP_
#define XYZ_GEOM_INTERSECTION_HPP_

namespace xyz {
//
// Intersection Parameter
//
class Intersection {
 public:
  // Init
  inline void Init() {
    t_min_ = float_t(1e-5);  // 10μm
    t_max_ = std::numeric_limits<float_t>::infinity();
  }

  inline float_t t_min() const { return t_min_; }
  inline float_t t_max() const { return t_max_; }
  inline float_t f() const { return f_; }
  inline float_t g() const { return g_; }
  inline bool is_back_side() const { return is_back_side_; }

  inline void set_t_min(__in float_t const &t_min) { t_min_ = t_min; }
  inline void set_t_max(__in float_t const &t_max) { t_max_ = t_max; }
  inline void set_f(__in float_t const &f) { f_ = f; }
  inline void set_g(__in float_t const &g) { g_ = g; }
  inline void set_is_back_side(__in bool const &is_back_side) {
    is_back_side_ = is_back_side;
  }

 private:
  // 光線の有効範囲
  float_t t_min_;
  float_t t_max_;

  // 交差パラメータ
  float_t f_;
  float_t g_;

  // 背面に当たったか
  bool is_back_side_;
};

}  // end of namespace xyz

#endif XYZ_GEOM_INTERSECTION_HPP_
