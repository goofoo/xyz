#ifndef XYZ_GEOM_PATHVERTEX_HPP_
#define XYZ_GEOM_PATHVERTEX_HPP_

#include "../core/config.hpp"
#include "../geom/triangle.hpp"

namespace xyz {
class Intersection;

class PathVertex {
 public:
  PathVertex() : pGeometry(reinterpret_cast<Triangle const *>(0xDEADBEAF)) {}

  Triangle const *pGeometry;  ///< 三角形

  bool bBackSide;  ///< 交点は面の裏か
  bool bSpecular;  ///< デルタ関数による反射か
  // bool bModified;

  float3_t vPosition;           ///< 交点の座標
  float3_t vIncomingDirection;  ///< 入射方向
  float3_t vOutgoingDirection;  ///< 射出方向
  float3_t vGeometricNormal;    ///< 幾何学的法線(geometric normal)
  float3_t vShadingNormal;      ///< シェーディング法線(shading noraml)
  float3_t vTangent;            ///< 接ベクトル
  float3_t vBinormal;           ///< 従法線ベクトル

  float3_t power;  ///< importance (potential) or power

  float_t fBSDFxIPDF;  ///< BSDF/PDF の値

  float_t fSamplingPrev;  ///<
                          ///この経路頂点から一つ前の経路頂点をサンプリングする射影立体角を測度とした確率(視点から光源の方向)
  float_t fSamplingNext;  ///<
                          ///この経路頂点から一つ次の経路頂点をサンプリングする射影立体角を測度とした確率(光源から視点の方向)

  float_t fIncomingCosThetaShading;
  float_t fOutgoingCosThetaGeometric;

  float_t fGeometricFactor;  ///< この経路頂点と一つ前の経路頂点間の幾何学項
                             ///\frac{\cos\theta \cdot \cos\theta}{r^{2}}

  inline void SetGeometricBasis(Intersection const &param) {
    pGeometry->GeometricBasis(param, &vShadingNormal, &vGeometricNormal,
                              &vTangent, &vBinormal);
  }

  inline void SetShadingBasis(Intersection const &param) {
    pGeometry->ShadingBasis(param, &vShadingNormal, &vGeometricNormal,
                            &vTangent, &vBinormal);
  }

  bool SetGeometricBasis(Intersection const &param,
                         PathVertex const &oldLightPathVertex);
  bool SetShadingBasis(Intersection const &param,
                       PathVertex const &oldEyePathVertex);
};
}
// end of namespace xyz

#endif XYZ_GEOM_PATHVERTEX_HPP_
