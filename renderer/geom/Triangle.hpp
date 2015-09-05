#ifndef __TGIR_GEOM_TRIANGLE_HPP__
#define __TGIR_GEOM_TRIANGLE_HPP__

//#define HI_TRIANGLE_HAS_MINMAX

namespace tgir {
class Intersection;

class Triangle {
  HI_DISALLOW_COPY_AND_ASSIGN(Triangle);

 public:
  typedef double value_type;
  typedef hi::basic_vector2<value_type> Vector2;
  typedef hi::basic_vector3<value_type> Vector3;

  Triangle(std::size_t const bsdf, Vector3 const &vp0, Vector3 const &vp1,
           Vector3 const &vp2, value_type const s0, value_type const t0,
           value_type const s1, value_type const t1, value_type const s2,
           value_type const t2);

  Triangle(std::size_t const bsdf, Vector3 const &vp0, Vector3 const &vp1,
           Vector3 const &vp2, Vector3 const &ns0, Vector3 const &ns1,
           Vector3 const &ns2, value_type const s0, value_type const t0,
           value_type const s1, value_type const t1, value_type const s2,
           value_type const t2);

  //
  // gui
  //
  void Render() const;

  //
  // kd-tree
  //
  inline void min(tgir::Vector3 &v) const;
  inline void max(tgir::Vector3 &v) const;
  inline value_type min(std::size_t const k) const;
  inline value_type max(std::size_t const k) const;

  int Overlap(std::size_t const axis, value_type const value) const;

  //
  // shading
  //
  inline tgir::Triangle const *FindIntersection(
      __in tgir::Vector3 const &vOrigin, __in tgir::Vector3 const &vDirection,
      __out tgir::Intersection *const pParam) const;

  inline std::size_t bsdf() const { return bsdf_; }

  inline value_type area() const {
    return hi::length(hi::cross(e2_, e1_)) * tgir::Real(0.5);
  }

  void GeometricBasis(__in tgir::Intersection const &param,
                      __out tgir::Vector3 *const pvShadingNormal,
                      __out tgir::Vector3 *const pvGeometricNormal,
                      __out tgir::Vector3 *const pvTangent,
                      __out tgir::Vector3 *const pvBinormal) const;

  void ShadingBasis(__in tgir::Intersection const &param,
                    __out tgir::Vector3 *const pvShadingNormal,
                    __out tgir::Vector3 *const pvGeometricNormal,
                    __out tgir::Vector3 *const pvTangent,
                    __out tgir::Vector3 *const pvBinormal) const;

  void Triangle::Sample(
      __in tgir::Real const &fXi1,            ///< [in]
      __in tgir::Real const &fXi2,            ///< [in]
      __out tgir::Vector3 *const pvPosition,  ///< [out] 三角形上のサンプル位置
      __out tgir::Vector3 *const pvShadingNormal,  ///< [out] シェーディング法線
      __out tgir::Vector3 *const pvGeometricNormal,  ///< [out] 幾何的法線
      __out tgir::Vector3 *const pvTangentVector,    ///< [out] 接線ベクトル
      __out tgir::Vector3 *const pvBinormalVector    ///< [out] 従法線ベクトル
      ) const;

 private:
  tgir::Triangle const *FindIntersectionCW(
      tgir::Vector3 const &vOrigin, tgir::Vector3 const &vDirection,
      tgir::Intersection *const pParam) const;

  tgir::Triangle const *FindIntersectionDF(
      tgir::Vector3 const &vOrigin, tgir::Vector3 const &vDirection,
      tgir::Intersection *const pParam) const;

 private:
  // material
  std::size_t const bsdf_;
  std::size_t const isDoubleFace_;

  // vertex positions
  tgir::Vector3 const v0_;
  tgir::Vector3 const e1_;
  tgir::Vector3 const e2_;

  // geometric normal vector (z-axis of local basis)
  tgir::Vector3 const ng_;

  // shading normal vectors
  tgir::Vector3 const n0_;
  tgir::Vector3 const s1_;
  tgir::Vector3 const s2_;

  // texture coordinates
  tgir::Vector2 const uv0_;
  tgir::Vector2 const uv1_;
  tgir::Vector2 const uv2_;

#ifdef HI_TRIANGLE_HAS_MINMAX
  // aabb (axis-aligned bounding box)
  tgir::Vector3 const min_;
  tgir::Vector3 const max_;
#endif HI_TRIANGLE_HAS_MINMAX
};

inline void Triangle::min(tgir::Vector3 &v) const {
#ifdef HI_TRIANGLE_HAS_MINMAX
  hi::min(v, min_);
#else
  v[0] = hi::min(v[0], v0_[0], v0_[0] + e1_[0], v0_[0] + e2_[0]);
  v[1] = hi::min(v[1], v0_[1], v0_[1] + e1_[1], v0_[1] + e2_[1]);
  v[2] = hi::min(v[2], v0_[2], v0_[2] + e1_[2], v0_[2] + e2_[2]);
#endif HI_TRIANGLE_HAS_MINMAX
}

inline void Triangle::max(tgir::Vector3 &v) const {
#ifdef HI_TRIANGLE_HAS_MINMAX
  hi::max(v, max_);
#else
  v[0] = hi::max(v[0], v0_[0], v0_[0] + e1_[0], v0_[0] + e2_[0]);
  v[1] = hi::max(v[1], v0_[1], v0_[1] + e1_[1], v0_[1] + e2_[1]);
  v[2] = hi::max(v[2], v0_[2], v0_[2] + e1_[2], v0_[2] + e2_[2]);
#endif HI_TRIANGLE_HAS_MINMAX
}

inline Triangle::value_type Triangle::min(std::size_t const k) const {
#ifdef HI_TRIANGLE_HAS_MINMAX
  return min_[k];
#else
  return hi::min(v0_[k], v0_[k] + e1_[k], v0_[k] + e2_[k]);
#endif HI_TRIANGLE_HAS_MINMAX
}

inline Triangle::value_type Triangle::max(std::size_t const k) const {
#ifdef HI_TRIANGLE_HAS_MINMAX
  return max_[k];
#else
  return hi::max(v0_[k], v0_[k] + e1_[k], v0_[k] + e2_[k]);
#endif HI_TRIANGLE_HAS_MINMAX
}

inline tgir::Triangle const *Triangle::FindIntersection(
    tgir::Vector3 const &vOrigin, tgir::Vector3 const &vDirection,
    tgir::Intersection *const pParam) const {
  return isDoubleFace_ ? FindIntersectionDF(vOrigin, vDirection, pParam)
                       : FindIntersectionCW(vOrigin, vDirection, pParam);
}

// Barycentric Coordinates
inline tgir::Vector3 Barycenter(tgir::Vector3 const &v0,
                                tgir::Vector3 const &e1,
                                tgir::Vector3 const &e2, tgir::Real const &f,
                                tgir::Real const &g) {
  return v0 + f * e1 + g * e2;
}

// {T, N, B}
inline tgir::Vector3 NormalVector(tgir::Vector3 const &e1,
                                  tgir::Vector3 const &e2) {
  return hi::cross(e2, e1);
}

// {T, N, B}
inline tgir::Vector3 TangentVector(tgir::Vector3 const &N,
                                   tgir::Vector3 const &B) {
  return hi::cross(N, B);
}

// {T, N, B}
inline tgir::Vector3 BinormalVector(tgir::Vector3 const &T,
                                    tgir::Vector3 const &N) {
  return hi::cross(T, N);
}

}  // end of namespace tgir

#endif __TGIR_GEOM_TRIANGLE_HPP__
