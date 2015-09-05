#include "pathvertex.hpp"
#include "intersection.hpp"

namespace xyz {
bool PathVertex::SetGeometricBasis(Intersection const &param,
                                   PathVertex const &oldLightPathVertex) {
  pGeometry->GeometricBasis(param, &vShadingNormal, &vGeometricNormal,
                            &vTangent, &vBinormal);

  if (-hi::dot(oldLightPathVertex.vOutgoingDirection, vShadingNormal) <= 0) {
    pGeometry = nullptr;
    return false;
  }

  power = oldLightPathVertex.power *
          hi::dot(oldLightPathVertex.vOutgoingDirection, vShadingNormal) /
          hi::dot(oldLightPathVertex.vOutgoingDirection,
                  vGeometricNormal);  // modify BSDF
  vIncomingDirection = -oldLightPathVertex.vOutgoingDirection;
  vPosition = oldLightPathVertex.vOutgoingDirection * param.t_max() +
              oldLightPathVertex.vPosition;
  fGeometricFactor =
      hi::dot(oldLightPathVertex.vOutgoingDirection,
              oldLightPathVertex.vGeometricNormal) *
      hi::dot(oldLightPathVertex.vOutgoingDirection, vShadingNormal) /
      hi::square_of(param.t_max());
  fGeometricFactor = std::abs(fGeometricFactor);  // 非負値化
  bBackSide = param.is_back_side();
  bSpecular = false;

  return true;
}

bool PathVertex::SetShadingBasis(Intersection const &param,
                                 PathVertex const &oldEyePathVertex) {
  pGeometry->ShadingBasis(param, &vShadingNormal, &vGeometricNormal, &vTangent,
                          &vBinormal);

  if (-hi::dot(oldEyePathVertex.vOutgoingDirection, vShadingNormal) <= 0) {
    return false;
  }

  power = oldEyePathVertex.power;
  vOutgoingDirection = oldEyePathVertex.vOutgoingDirection;
  // vIncomingDirection  =-oldEyePathVertex.vOutgoingDirection;
  vPosition = oldEyePathVertex.vOutgoingDirection * param.t_max() +
              oldEyePathVertex.vPosition;
  fGeometricFactor =
      hi::dot(oldEyePathVertex.vOutgoingDirection,
              oldEyePathVertex.vShadingNormal) *
      hi::dot(oldEyePathVertex.vOutgoingDirection, vGeometricNormal) /
      hi::square_of(param.t_max());
  fGeometricFactor = std::abs(fGeometricFactor);  // 非負値化
  bBackSide = param.is_back_side();
  bSpecular = false;

  return true;
}
}
// end of namespace xyz
