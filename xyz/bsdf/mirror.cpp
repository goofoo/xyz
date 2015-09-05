#include "mirror.hpp"
#include "../core/scene.hpp"
#include "../geom/pathvertex.hpp"
#include "../geom/intersection.hpp"

namespace xyz {
void Mirror::Render() const {
  hi::basic_vector4<GLfloat> const color(0, 0.5483f, 0.6549f,
                                         1);  // Patch 18, Cyan
  GLfloat const shininess = 20;
  ::glEnable(GL_LIGHTING);
  ::glDisable(GL_BLEND);
  hi::glMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color.data());
  hi::glMaterial(GL_FRONT_AND_BACK, GL_SPECULAR, color.data());
  hi::glMaterial(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
}

bool Mirror::BoundsImportance(PathVertex *const pLightPathVertex,
                              std::mt19937_64 &random) const {
  // vOutgoingDirection には 入射してきた方向が入っている
  float_t const fCosTheta_s = -hi::dot(pLightPathVertex->vOutgoingDirection,
                                       pLightPathVertex->vShadingNormal);
  if (fCosTheta_s <= 0) {
    return false;
  }
  float_t const fCosTheta_g = -hi::dot(pLightPathVertex->vOutgoingDirection,
                                       pLightPathVertex->vGeometricNormal);

  // 反射ベクトルを求める
  ComputeReflectedVector(pLightPathVertex->vShadingNormal, fCosTheta_s,
                         &pLightPathVertex->vOutgoingDirection);

  float_t const fCosTheta_o = hi::dot(pLightPathVertex->vOutgoingDirection,
                                      pLightPathVertex->vGeometricNormal);
  if (fCosTheta_o <= 0) {
    return false;
  }

  // NOTE: ここでエネルギーが増えている可能性がある
  float_t const fShadingFactor = fCosTheta_o / fCosTheta_g;
  if (random.next<float_t>() >= fShadingFactor) {
    return false;
  }

  pLightPathVertex->bSpecular = true;

  return true;
}

float3_t Mirror::CalculateWeightedDirectLighting(__inout PathVertex *const,
                                                 __in PathVertex const &,
                                                 __in Scene const &) const {
  return float3_t(0);
}

bool Mirror::NextScatteringDirection(__inout PathVertex *const pPathVertex,
                                     __inout IPrimarySample &sample) const {
  UNREFERENCED_PARAMETER(sample);

  // 入力時のpathVertex.vOutgoingDirectionには，IncomingDirectionが格納されている．
  float_t const cosTheta = -hi::dot(pPathVertex->vOutgoingDirection,
                                    pPathVertex->vShadingNormal);  // \cos\theta

  // 反射ベクトルを求める
  ComputeReflectedVector(pPathVertex->vShadingNormal, cosTheta,
                         &pPathVertex->vOutgoingDirection);
  if (hi::dot(pPathVertex->vOutgoingDirection, pPathVertex->vGeometricNormal) <=
      0) {
    return false;  // 打ち止め
  }

  /*
      // 薄膜の厚さ
      static float_t const thinness = 1090;

      // 波長
      float_t const lambda = 400 + (700 - 400) * wavelength;

      // 薄膜の屈折率
      float_t const ior = 1.400 + 50.0 / (lambda - 230);
      float_t const iorSquared = hi::square_of(ior);

      float_t const sinThetaSquared = float_t(1) - hi::square_of(cosTheta); //
     \sin^2\theta

      float_t const delta = (2*M_PI) * thinness / lambda * std::sqrt(iorSquared
     - sinThetaSquared); // \mu_{1} = 1 としている(真空なので)

      // 金属反射なので \mu_{1} < \mu_{f} < \mu_{2} よって deletaDash = 0
      // \mu_{1} = 1;
      // \mu_{f} = 1.4 以上
      // \mu_{2} = 金属なので100ぐらい?
      static float_t const deletaDash = 0;

      float_t const cosDeltaSquared = hi::square_of(std::cos(delta+deletaDash));

      pathVertex.power *= cosDeltaSquared;
  */

  pPathVertex->bSpecular = true;
  pPathVertex->fSamplingPrev = 1;
  pPathVertex->fSamplingNext = 1;

  return true;
}

float_t Mirror::GetDensityVariance() const { return 1; }

bool Mirror::NextScatteringDirection(PathVertex &pathVertex, float3_t const &,
                                     bool const &) const {
  // 入力時のpathVertex.vOutgoingDirectionには，IncomingDirectionが格納されている．
  float_t const cosTheta = -hi::dot(pathVertex.vOutgoingDirection,
                                    pathVertex.vShadingNormal);  // \cos\theta

  // 反射ベクトルを求める
  ComputeReflectedVector(pathVertex.vShadingNormal, cosTheta,
                         &pathVertex.vOutgoingDirection);

  if (hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal) <=
      0) {
    return false;  // 打ち止め
  }

#if 0
    // 薄膜の厚さ
    static float_t const thinness = 1090;

    // 波長
    float_t const lambda = 400 + (700 - 400) * sample.WavelengthMu();

    // 薄膜の屈折率
    float_t const ior = 1.400 + 50.0 / (lambda - 230);
    float_t const iorSquared = hi::square_of(ior);

    float_t const sinThetaSquared = float_t(1) - hi::square_of(cosTheta); // \sin^2\theta

    float_t const delta = (2*M_PI) * thinness / lambda * std::sqrt(iorSquared - sinThetaSquared); // \mu_{1} = 1 としている(真空なので)

    // 金属反射なので \mu_{1} < \mu_{f} < \mu_{2} よって deletaDash = 0
    // \mu_{1} = 1;
    // \mu_{f} = 1.4 以上
    // \mu_{2} = 金属なので100ぐらい?
    static float_t const deletaDash = 0;

    float_t const cosDeltaSquared = hi::square_of(std::cos(delta+deletaDash));

    pathVertex.power *= cosDeltaSquared;
#endif

  pathVertex.bSpecular = true;
  pathVertex.fSamplingPrev = 1;
  pathVertex.fSamplingNext = 1;

  return true;
}

}  // end of namespace xyz

namespace xyz {
// pPathVertex には，交点とその基底ベクトル，vIncomingDirection，
// および前の頂点と同じ power の値が設定されているものとする．
//   pPathVertex->vIncomingDirection = -pPreviousPathVertex->vOutgoingDirection;
//   pPathVertex->power              =  pPreviousPathVertex->power;
// また，条件として[hi::dot(pPathVertex->vIncomingDirection,
// pPathVertex->vShadingNormal) > 0]を課す．
//
// pPathVertex->fIncomingCosThetaShading は設定済み
// pPathVertex->fOutgoingCosThetaGeometric を設定する．
//
// NOTE: Incoming や Outgoing は光エネルギーの流れる方向を示しており
//       実装上の光線追跡の方向を示すものではない
bool Mirror::LightsToCameraScatteringDirection(
    __inout PathVertex *const pPathVertex,
    __inout IPrimarySample &sample) const {
  UNREFERENCED_PARAMETER(sample);

  pPathVertex->bSpecular = true;

  pPathVertex->vOutgoingDirection = -pPathVertex->vIncomingDirection;
  ComputeReflectedVector(pPathVertex->vShadingNormal,
                         pPathVertex->fIncomingCosThetaShading,
                         &(pPathVertex->vOutgoingDirection));

  // 幾何学的法線と散乱方向ベクトルの内積を求める．
  pPathVertex->fOutgoingCosThetaGeometric =
      hi::dot(pPathVertex->vOutgoingDirection, pPathVertex->vGeometricNormal);
  if (pPathVertex->fOutgoingCosThetaGeometric <= 0) {
    return false;
  }

  pPathVertex->fSamplingPrev = hi::rcp(pPathVertex->fIncomingCosThetaShading);
  pPathVertex->fSamplingNext = hi::rcp(pPathVertex->fOutgoingCosThetaGeometric);

  // 幾何学的法線と入射方向ベクトルの内積を求める．
  float_t const fIncomingCosThetaGeometric =
      hi::dot(pPathVertex->vIncomingDirection, pPathVertex->vGeometricNormal);

  pPathVertex->fBSDFxIPDF =
      pPathVertex->fOutgoingCosThetaGeometric / fIncomingCosThetaGeometric;
  // pPathVertex->power *= albedo_;

  return true;
}

// pPathVertex には，交点とその基底ベクトル，vOutgoingDirection，
// および前の頂点と同じ power の値が設定されているものとする．
//   pPathVertex->vOutgoingDirection = -pPreviousPathVertex->vIncomingDirection;
//   pPathVertex->power              =  pPreviousPathVertex->power;
// また，条件として[hi::dot(pPathVertex->vOutgingDirection,
// pPathVertex->vShadingNormal) > 0]を課す．
//
// pPathVertex->fOutgoingCosThetaGeometric は設定済み
// pPathVertex->fIncomingCosThetaShading を設定する．
//
// NOTE: Incoming や Outgoing は光エネルギーの流れる方向を示しており
//       実装上の光線追跡の方向を示すものではない
bool Mirror::CameraToLightsScatteringDirection(
    __inout PathVertex *const pPathVertex,
    __inout IPrimarySample &sample) const {
  UNREFERENCED_PARAMETER(sample);

  pPathVertex->bSpecular = true;

  // 鏡面反射では(fIncomingCosThetaShading ==
  // fOutgoingCosThetaShading)だからここで設定する．
  pPathVertex->fIncomingCosThetaShading =
      hi::dot(pPathVertex->vOutgoingDirection, pPathVertex->vShadingNormal);

  // 鏡面反射方向をサンプリングする
  pPathVertex->vIncomingDirection = -pPathVertex->vOutgoingDirection;
  ComputeReflectedVector(pPathVertex->vShadingNormal,
                         pPathVertex->fIncomingCosThetaShading,
                         &(pPathVertex->vIncomingDirection));

  float_t const fIncomingCosThetaGeometric =
      hi::dot(pPathVertex->vIncomingDirection, pPathVertex->vGeometricNormal);
  if (fIncomingCosThetaGeometric > 0) {
    return false;
  }

  pPathVertex->fSamplingPrev = hi::rcp(pPathVertex->fIncomingCosThetaShading);
  pPathVertex->fSamplingNext = hi::rcp(pPathVertex->fOutgoingCosThetaGeometric);

  pPathVertex->fBSDFxIPDF = 1;
  // pPathVertex->power *= albedo_;

  return true;
}
}  // end of namespace xyz
