#include "PathVertex.hpp"
#include "Intersection.hpp"

namespace tgir {
bool PathVertex::SetGeometricBasis(tgir::Intersection const &param,
                                   tgir::PathVertex const &oldLightPathVertex) {
  pGeometry->GeometricBasis(param, &vShadingNormal, &vGeometricNormal,
                            &vTangent, &vBinormal);

  if (-hi::dot(oldLightPathVertex.vOutgoingDirection, vShadingNormal) <= 0) {
    pGeometry = nullptr;
    return false;
  }

  sQuantum = oldLightPathVertex.sQuantum *
             hi::dot(oldLightPathVertex.vOutgoingDirection, vShadingNormal) /
             hi::dot(oldLightPathVertex.vOutgoingDirection,
                     vGeometricNormal);  // modify BSDF
  vIncomingDirection = -oldLightPathVertex.vOutgoingDirection;
  vPosition = oldLightPathVertex.vOutgoingDirection * param.t_max() +
              oldLightPathVertex.vPosition;
  rGeometricFactor =
      hi::dot(oldLightPathVertex.vOutgoingDirection,
              oldLightPathVertex.vGeometricNormal) *
      hi::dot(oldLightPathVertex.vOutgoingDirection, vShadingNormal) /
      hi::square_of(param.t_max());
  rGeometricFactor = std::abs(rGeometricFactor);  // 非負値化
  back_side = param.is_back_side();
  bSpecular = false;

  return true;
}

bool PathVertex::SetShadingBasis(tgir::Intersection const &param,
                                 tgir::PathVertex const &oldEyePathVertex) {
  pGeometry->ShadingBasis(param, &vShadingNormal, &vGeometricNormal, &vTangent,
                          &vBinormal);

  if (-hi::dot(oldEyePathVertex.vOutgoingDirection, vShadingNormal) <= 0) {
    return false;
  }

  sQuantum = oldEyePathVertex.sQuantum;
  vOutgoingDirection = oldEyePathVertex.vOutgoingDirection;
  // vIncomingDirection  =-oldEyePathVertex.vOutgoingDirection;
  vPosition = oldEyePathVertex.vOutgoingDirection * param.t_max() +
              oldEyePathVertex.vPosition;
  rGeometricFactor =
      hi::dot(oldEyePathVertex.vOutgoingDirection,
              oldEyePathVertex.vShadingNormal) *
      hi::dot(oldEyePathVertex.vOutgoingDirection, vGeometricNormal) /
      hi::square_of(param.t_max());
  rGeometricFactor = std::abs(rGeometricFactor);  // 非負値化
  back_side = param.is_back_side();
  bSpecular = false;

  return true;
}
}
// end of namespace tgir
