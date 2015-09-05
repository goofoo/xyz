#ifndef __TGIR_GEOM_PATHVERTEX_HPP__
#define __TGIR_GEOM_PATHVERTEX_HPP__

#include "core/config.hpp"
#include "geom/Triangle.hpp"

namespace tgir {
class Intersection;

class PathVertex {
 public:
  tgir::Triangle const *pGeometry;  ///< 三角形

  tgir::Vector3 vPosition;           ///< 交点の座標
  tgir::Vector3 vIncomingDirection;  ///< 入射方向
  tgir::Vector3 vOutgoingDirection;  ///< 射出方向
  tgir::Vector3 vGeometricNormal;    ///< 幾何学的法線(geometric normal)
  tgir::Vector3 vShadingNormal;      ///< シェーディング法線(shading noraml)
  tgir::Vector3 vTangent;            ///< 接ベクトル
  tgir::Vector3 vBinormal;           ///< 従法線ベクトル

  tgir::Spectrum sQuantum;  ///< importance (potential) or power

  bool back_side;  ///< 交点は面の裏か
  bool bSpecular;  ///< デルタ関数による反射か
  // bool bModified;

  tgir::Real rSamplingPrev;  ///<
  ///この経路頂点から一つ前の経路頂点をサンプリングする射影立体角を測度とした確率
  tgir::Real rSamplingNext;  ///<
  ///この経路頂点から一つ次の経路頂点をサンプリングする射影立体角を測度とした確率
  tgir::Real rGeometricFactor;  ///< この経路頂点と一つ前の経路頂点間の幾何学項
                                ///\frac{\cos\theta \cdot \cos\theta}{r^{2}}

  inline void SetGeometricBasis(tgir::Intersection const &param) {
    pGeometry->GeometricBasis(param, &vShadingNormal, &vGeometricNormal,
                              &vTangent, &vBinormal);
  }

  inline void SetShadingBasis(tgir::Intersection const &param) {
    pGeometry->ShadingBasis(param, &vShadingNormal, &vGeometricNormal,
                            &vTangent, &vBinormal);
  }

  bool SetGeometricBasis(tgir::Intersection const &param,
                         tgir::PathVertex const &oldLightPathVertex);
  bool SetShadingBasis(tgir::Intersection const &param,
                       tgir::PathVertex const &oldEyePathVertex);
};
}
// end of namespace tgir

#endif __TGIR_GEOM_PATHVERTEX_HPP__
