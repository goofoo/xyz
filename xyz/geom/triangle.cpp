#include "triangle.hpp"
#include "intersection.hpp"
#include "../core/config.hpp"

namespace xyz {
/// basis vectors of right-hand coordinate system: [T, N, B]
///   tangent vector: T, normal vector: N, bitangent vector: B = TxN
/// 定義: B = cross(T, N)

Triangle::Triangle(std::size_t const mat, float3_t const &vp0,
                   float3_t const &vp1, float3_t const &vp2,
                   value_type const s0, value_type const t0,
                   value_type const s1, value_type const t1,
                   value_type const s2, value_type const t2)
    : bsdf_(mat),
      isDoubleFace_(true),
      v0_(vp0),
      e1_(vp1 - vp0),
      e2_(vp2 - vp0),
      ng_(hi::normalize(NormalVector(e1_, e2_))),
      n0_(ng_),
      s1_(0),
      s2_(0),
      uv0_(s0, t0),
      uv1_(s1, t1),
      uv2_(s2, t2)
#ifdef HI_TRIANGLE_HAS_MINMAX
      ,
      min_(hi::min(hi::min(vp0, vp1), vp2)),
      max_(hi::max(hi::max(vp0, vp1), vp2))
#endif HI_TRIANGLE_HAS_MINMAX
{
}

Triangle::Triangle(std::size_t const mat, float3_t const &vp0,
                   float3_t const &vp1, float3_t const &vp2,
                   float3_t const &ns0, float3_t const &ns1,
                   float3_t const &ns2, value_type const s0,
                   value_type const t0, value_type const s1,
                   value_type const t1, value_type const s2,
                   value_type const t2)
    : bsdf_(mat),
      isDoubleFace_(true),
      v0_(vp0),
      e1_(vp1 - vp0),
      e2_(vp2 - vp0),
      ng_(hi::normalize(NormalVector(e1_, e2_))),
      n0_(ns0),
      s1_(ns1 - ns0),
      s2_(ns2 - ns0),
      uv0_(s0, t0),
      uv1_(s1, t1),
      uv2_(s2, t2)
#ifdef HI_TRIANGLE_HAS_MINMAX
      ,
      min_(hi::min(hi::min(vp0, vp1), vp2)),
      max_(hi::max(hi::max(vp0, vp1), vp2))
#endif HI_TRIANGLE_HAS_MINMAX
{
}

void Triangle::Render() const {
  ::glBegin(GL_TRIANGLES);
  {
    hi::glNormal(n0_.data());
    hi::glVertex<3>(v0_.data());
    hi::glNormal(n0_[0] + s1_[0], n0_[1] + s1_[1], n0_[2] + s1_[2]);
    hi::glVertex(v0_[0] + e1_[0], v0_[1] + e1_[1], v0_[2] + e1_[2]);
    hi::glNormal(n0_[0] + s2_[0], n0_[1] + s2_[1], n0_[2] + s2_[2]);
    hi::glVertex(v0_[0] + e2_[0], v0_[1] + e2_[1], v0_[2] + e2_[2]);
  }
  ::glEnd();
}

int Triangle::Overlap(std::size_t const axis, value_type const value) const {
  if (max(axis) < value) {
    return 1;
  }

  if (value <= min(axis)) {
    return 2;
  }

  return 3;
}

Triangle const *Triangle::FindIntersectionCW(
    __in xyz::float3_t const &vOrigin, __in xyz::float3_t const &vDirection,
    __out Intersection *const pParam) const {
  typedef xyz::float3_t vector;
  typedef vector::value_type real;

  vector const e2xe1 = hi::cross(e2_, e1_);

  real const det = hi::dot(vDirection, e2xe1);
  if (det > EPSILON) {
    return nullptr;  // 裏側
  }
  real const inv_det = hi::rcp(det);

  vector const vmo = v0_ - vOrigin;
  real const t = inv_det * hi::dot(vmo, e2xe1);
  if ((t < pParam->t_min()) || (pParam->t_max() < t)) {
    return nullptr;
  }

  vector const dxvmo = hi::cross(vDirection, vmo);

  real const f = inv_det * hi::dot(e2_, dxvmo);
  if ((f < 0) || (f > 1)) {
    return nullptr;
  }

  real const g = inv_det * -hi::dot(e1_, dxvmo);
  if ((g < 0) || ((f + g) > 1)) {
    return nullptr;
  }

  pParam->set_t_max(t);
  pParam->set_f(f);
  pParam->set_g(g);
  pParam->set_is_back_side(false);

  return this;
}

Triangle const *Triangle::FindIntersectionDF(
    __in xyz::float3_t const &vOrigin, __in xyz::float3_t const &vDirection,
    __out Intersection *const pParam) const {
  typedef xyz::float3_t vector;
  typedef vector::value_type real;

  vector const e2xe1 = hi::cross(e2_, e1_);

  real const det = hi::dot(vDirection, e2xe1);
  if (std::fabs(det) < EPSILON) {
    return nullptr;
  }
  real const inv_det = hi::rcp(det);

  vector const vmo = v0_ - vOrigin;
  real const t = inv_det * hi::dot(vmo, e2xe1);
  if ((t < pParam->t_min()) || (pParam->t_max() < t)) {
    return nullptr;
  }

  vector const dxvmo = hi::cross(vDirection, vmo);

  real const f = inv_det * hi::dot(e2_, dxvmo);
  if ((f < 0) || (f > 1)) {
    return nullptr;
  }

  real const g = inv_det * -hi::dot(e1_, dxvmo);
  if ((g < 0) || ((f + g) > 1)) {
    return nullptr;
  }

  pParam->set_t_max(t);
  pParam->set_f(f);
  pParam->set_g(g);
  pParam->set_is_back_side(det > 0);  // 裏側

  return this;
}

void Triangle::GeometricBasis(__in Intersection const &param,
                              __out xyz::float3_t *const pvShadingNormal,
                              __out xyz::float3_t *const pvGeometricNormal,
                              __out xyz::float3_t *const pvTangent,
                              __out xyz::float3_t *const pvBinormal) const {
  *pvShadingNormal =
      hi::normalize(Barycenter(n0_, s1_, s2_, param.f(), param.g()));
  *pvGeometricNormal = ng_;

  if (param.is_back_side()) {
    *pvGeometricNormal = -*pvGeometricNormal;
    *pvShadingNormal = -*pvShadingNormal;
  }

  *pvTangent = float3_t(-(*pvGeometricNormal)[2], 0, (*pvGeometricNormal)[0]);
  value_type const fLengthSquared = hi::length_squared(*pvTangent);
  if (fLengthSquared > 0) {
    // 接ベクトルから求める
    *pvTangent *= hi::rsqrt(fLengthSquared);
    *pvBinormal = BinormalVector(*pvTangent, *pvGeometricNormal);
  } else {
    // 従接ベクトルから求める
    *pvBinormal = float3_t(0, 0, ((*pvGeometricNormal)[1] > 0) ? 1 : -1);
    *pvTangent = TangentVector(*pvGeometricNormal, *pvBinormal);
  }
}

void Triangle::ShadingBasis(__in Intersection const &param,
                            __out xyz::float3_t *const pvShadingNormal,
                            __out xyz::float3_t *const pvGeometricNormal,
                            __out xyz::float3_t *const pvTangent,
                            __out xyz::float3_t *const pvBinormal) const {
  *pvShadingNormal =
      hi::normalize(Barycenter(n0_, s1_, s2_, param.f(), param.g()));
  *pvGeometricNormal = ng_;

  if (param.is_back_side()) {
    *pvGeometricNormal = -*pvGeometricNormal;
    *pvShadingNormal = -*pvShadingNormal;
  }

  *pvTangent = float3_t(-(*pvShadingNormal)[2], 0, (*pvShadingNormal)[0]);
  value_type const fLengthSquared = hi::length_squared(*pvTangent);
  if (fLengthSquared > 0) {
    // 接ベクトルから求める
    *pvTangent *= hi::rsqrt(fLengthSquared);
    *pvBinormal = BinormalVector(*pvTangent, *pvShadingNormal);
  } else {
    // 従接ベクトルから求める
    *pvBinormal = float3_t(0, 0, ((*pvShadingNormal)[1] > 0) ? 1 : -1);
    *pvTangent = TangentVector(*pvShadingNormal, *pvBinormal);
  }
}

/*!
 * SamplePlanarTriangle(real xi_1,  real xi_2)
 * {
 *   Compute the warping function (xi_1, xi_2) -> (s, t).
 *     s <- sqrt(xi_1);
 *     t <-      xi_2 ;
 *   Plug the warped coords into the original parametrization.
 *     P <- (1-s)*A + s*(1-t)*B + s*t*C;
 *     return P;
 * }
 */
void Triangle::Sample(
    __in xyz::float_t const &fXi1,          ///< [in]
    __in xyz::float_t const &fXi2,          ///< [in]
    __out xyz::float3_t *const pvPosition,  ///< [out] 三角形上のサンプル位置
    __out xyz::float3_t *const pvShadingNormal,  ///< [out] シェーディング法線
    __out xyz::float3_t *const pvGeometricNormal,  ///< [out] 幾何学的法線
    __out xyz::float3_t *const pvTangentVector,    ///< [out] 接線ベクトル
    __out xyz::float3_t *const pvBinormalVector  ///< [out] 従法線ベクトル
    ) const {
  float_t const s = std::sqrt(fXi1);
  float_t const t = (fXi2);
  float_t const u = s * (1 - t);
  float_t const v = s * (t);

  *pvPosition = Barycenter(v0_, e1_, e2_, u, v);
  *pvShadingNormal = hi::normalize(Barycenter(n0_, s1_, s2_, u, v));
  *pvGeometricNormal = ng_;

  // 球面座標系をベースに基底を求める
  *pvTangentVector =
      float3_t(-(*pvGeometricNormal)[2], 0, (*pvGeometricNormal)[0]);
  value_type const fLengthSquared = hi::length_squared(*pvTangentVector);
  if (fLengthSquared > 0) {
    // 接ベクトルから求める
    *pvTangentVector *= hi::rsqrt(fLengthSquared);
    *pvBinormalVector = BinormalVector(*pvTangentVector, *pvGeometricNormal);
  } else {
    // 従接ベクトルから求める
    *pvBinormalVector = float3_t(0, 0, ((*pvGeometricNormal)[1] > 0) ? 1 : -1);
    *pvTangentVector = TangentVector(*pvGeometricNormal, *pvBinormalVector);
  }
}

}  // end of namespace xyz
