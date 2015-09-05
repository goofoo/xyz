#ifndef XYZ_GEOM_TRIANGLE_HPP_
#define XYZ_GEOM_TRIANGLE_HPP_

//#define HI_TRIANGLE_HAS_MINMAX

namespace xyz {
class Intersection;

class Triangle {
  HI_DISALLOW_COPY_AND_ASSIGN(Triangle);

 public:
  typedef float value_type;
  typedef hi::basic_vector2<value_type> float2_t;
  typedef hi::basic_vector3<value_type> float3_t;

  Triangle(std::size_t const bsdf, float3_t const &vp0, float3_t const &vp1,
           float3_t const &vp2, value_type const s0, value_type const t0,
           value_type const s1, value_type const t1, value_type const s2,
           value_type const t2);

  Triangle(std::size_t const bsdf, float3_t const &vp0, float3_t const &vp1,
           float3_t const &vp2, float3_t const &ns0, float3_t const &ns1,
           float3_t const &ns2, value_type const s0, value_type const t0,
           value_type const s1, value_type const t1, value_type const s2,
           value_type const t2);

  //
  // gui
  //
  void Render() const;

  //
  // kd-tree
  //
  inline void min(xyz::float3_t &v) const;
  inline void max(xyz::float3_t &v) const;
  inline value_type min(std::size_t const k) const;
  inline value_type max(std::size_t const k) const;
  inline float3_t min() const;
  inline float3_t max() const;

  inline int morton_code(xyz::float3_t const &min,
                         xyz::float3_t const &range) const;

  int Overlap(std::size_t const axis, value_type const value) const;

  //
  // shading
  //
  inline Triangle const *FindIntersection(
      __in xyz::float3_t const &vOrigin, __in xyz::float3_t const &vDirection,
      __out Intersection *const pParam) const;

  inline std::size_t bsdf() const { return bsdf_; }

  inline value_type area() const {
    return hi::length(hi::cross(e2_, e1_)) * float_t(0.5);
  }

  void GeometricBasis(__in Intersection const &param,
                      __out xyz::float3_t *const pvShadingNormal,
                      __out xyz::float3_t *const pvGeometricNormal,
                      __out xyz::float3_t *const pvTangent,
                      __out xyz::float3_t *const pvBinormal) const;

  void ShadingBasis(__in Intersection const &param,
                    __out xyz::float3_t *const pvShadingNormal,
                    __out xyz::float3_t *const pvGeometricNormal,
                    __out xyz::float3_t *const pvTangent,
                    __out xyz::float3_t *const pvBinormal) const;

  void Sample(
      __in xyz::float_t const &fXi1,          ///< [in]
      __in xyz::float_t const &fXi2,          ///< [in]
      __out xyz::float3_t *const pvPosition,  ///< [out] 三角形上のサンプル位置
      __out xyz::float3_t *const pvShadingNormal,  ///< [out] シェーディング法線
      __out xyz::float3_t *const pvGeometricNormal,  ///< [out] 幾何学的法線
      __out xyz::float3_t *const pvTangentVector,    ///< [out] 接線ベクトル
      __out xyz::float3_t *const pvBinormalVector    ///< [out] 従法線ベクトル
      ) const;

 private:
  Triangle const *FindIntersectionCW(xyz::float3_t const &vOrigin,
                                     xyz::float3_t const &vDirection,
                                     Intersection *const pParam) const;

  Triangle const *FindIntersectionDF(xyz::float3_t const &vOrigin,
                                     xyz::float3_t const &vDirection,
                                     Intersection *const pParam) const;

 private:
  // material
  std::size_t const bsdf_;
  std::size_t const isDoubleFace_;

  // vertex positions
  float3_t const v0_;
  float3_t const e1_;
  float3_t const e2_;

  // geometric normal vector (z-axis of local basis)
  float3_t const ng_;

  // shading normal vectors
  float3_t const n0_;
  float3_t const s1_;
  float3_t const s2_;

  // texture coordinates
  float2_t const uv0_;
  float2_t const uv1_;
  float2_t const uv2_;

#ifdef HI_TRIANGLE_HAS_MINMAX
  // aabb (axis-aligned bounding box)
  float3_t const min_;
  float3_t const max_;
#endif HI_TRIANGLE_HAS_MINMAX
};

inline void Triangle::min(xyz::float3_t &v) const {
#ifdef HI_TRIANGLE_HAS_MINMAX
  hi::min(v, min_);
#else
  v[0] = hi::min<xyz::float_t>(v[0], v0_[0], v0_[0] + e1_[0], v0_[0] + e2_[0]);
  v[1] = hi::min<xyz::float_t>(v[1], v0_[1], v0_[1] + e1_[1], v0_[1] + e2_[1]);
  v[2] = hi::min<xyz::float_t>(v[2], v0_[2], v0_[2] + e1_[2], v0_[2] + e2_[2]);
#endif HI_TRIANGLE_HAS_MINMAX
}

inline void Triangle::max(xyz::float3_t &v) const {
#ifdef HI_TRIANGLE_HAS_MINMAX
  hi::max(v, max_);
#else
  v[0] = hi::max<xyz::float_t>(v[0], v0_[0], v0_[0] + e1_[0], v0_[0] + e2_[0]);
  v[1] = hi::max<xyz::float_t>(v[1], v0_[1], v0_[1] + e1_[1], v0_[1] + e2_[1]);
  v[2] = hi::max<xyz::float_t>(v[2], v0_[2], v0_[2] + e1_[2], v0_[2] + e2_[2]);
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

inline Triangle::float3_t Triangle::min() const {
#ifdef HI_TRIANGLE_HAS_MINMAX
  return min_;
#else
  return Triangle::float3_t(this->min(0), this->min(1), this->min(2));
#endif
}

inline Triangle::float3_t Triangle::max() const {
#ifdef HI_TRIANGLE_HAS_MINMAX
  return max_;
#else
  return Triangle::float3_t(this->max(0), this->max(1), this->max(2));
#endif
}

inline int Triangle::morton_code(xyz::float3_t const &min,
                                 xyz::float3_t const &range) const {
  float3_t const center((this->min() + this->max()) * value_type(0.5));
  float3_t const center_normalized((center - min) / range);

  int const x = static_cast<int>(center_normalized[0] * 1024);
  int const y = static_cast<int>(center_normalized[1] * 1024);
  int const z = static_cast<int>(center_normalized[2] * 1024);

  return hi::morton_code(x, y, z);
}

inline Triangle const *Triangle::FindIntersection(
    xyz::float3_t const &vOrigin, xyz::float3_t const &vDirection,
    Intersection *const pParam) const {
  return isDoubleFace_ ? FindIntersectionDF(vOrigin, vDirection, pParam)
                       : FindIntersectionCW(vOrigin, vDirection, pParam);
}

// Barycentric Coordinates
inline xyz::float3_t Barycenter(xyz::float3_t const &v0,
                                xyz::float3_t const &e1,
                                xyz::float3_t const &e2, xyz::float_t const &f,
                                xyz::float_t const &g) {
  return v0 + f * e1 + g * e2;
}

// {T, N, B}
inline float3_t NormalVector(xyz::float3_t const &e1, xyz::float3_t const &e2) {
  return hi::cross(e2, e1);
}

// {T, N, B}
inline float3_t TangentVector(xyz::float3_t const &N, xyz::float3_t const &B) {
  return hi::cross(N, B);
}

// {T, N, B}
inline float3_t BinormalVector(xyz::float3_t const &T, xyz::float3_t const &N) {
  return hi::cross(T, N);
}

}  // end of namespace xyz

#endif XYZ_GEOM_TRIANGLE_HPP_
