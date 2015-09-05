#ifndef __TGIR_GEOM_INTERSECTION_HPP__
#define __TGIR_GEOM_INTERSECTION_HPP__

namespace tgir {
//
// Intersection Parameter
//
class Intersection {
 public:
  // Init
  inline void Init() {
    t_min_ = tgir::Real(1e-5);  // 10μm
    t_max_ = std::numeric_limits<tgir::Real>::infinity();
  }

  inline tgir::Real t_min() const { return t_min_; }
  inline tgir::Real t_max() const { return t_max_; }
  inline tgir::Real f() const { return f_; }
  inline tgir::Real g() const { return g_; }
  inline bool is_back_side() const { return is_back_side_; }

  inline void set_t_min(__in tgir::Real const &t_min) { t_min_ = t_min; }
  inline void set_t_max(__in tgir::Real const &t_max) { t_max_ = t_max; }
  inline void set_f(__in tgir::Real const &f) { f_ = f; }
  inline void set_g(__in tgir::Real const &g) { g_ = g; }
  inline void set_is_back_side(__in bool const &is_back_side) {
    is_back_side_ = is_back_side;
  }

 private:
  // 光線の有効範囲
  tgir::Real t_min_;
  tgir::Real t_max_;

  // 交差パラメータ
  tgir::Real f_;
  tgir::Real g_;

  // 背面に当たったか
  bool is_back_side_;
};

}  // end of namespace tgir

#endif __TGIR_GEOM_INTERSECTION_HPP__
