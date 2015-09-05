#include "metal.hpp"
#include "../core/scene.hpp"
#include "../geom/pathvertex.hpp"
#include "../geom/intersection.hpp"

namespace xyz {
//
// 金属
//
void Metal::Render() const {
  hi::basic_vector3<float> rgb;
  hi::CCIR601_1_XYZ2RGB(albedo_, rgb);

  ::glEnable(GL_LIGHTING);
  ::glDisable(GL_BLEND);
  hi::glMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, rgb.data());
  hi::glMaterial(GL_FRONT_AND_BACK, GL_SPECULAR, rgb.data());
  hi::glMaterial(GL_FRONT_AND_BACK, GL_SHININESS,
                 static_cast<GLfloat>(shininess_));
}

bool Metal::BoundsImportance(PathVertex *const pLightPathVertex,
                             std::mt19937_64 &random) const {
  float_t const nk1 = -hi::dot(pLightPathVertex->vOutgoingDirection,
                               pLightPathVertex->vShadingNormal);
  if (nk1 <= 0.0) {
    return false;
  }
  float_t const fCosTheta_g = -hi::dot(pLightPathVertex->vOutgoingDirection,
                                       pLightPathVertex->vGeometricNormal);

  // ハーフベクトルをサンプリング(シェーディング法線を使用)
  float3_t vHalfVector;
  ComputeGlossyHalfwayVector(random.next<float_t>(), random.next<float_t>(),
                             shininess_, pLightPathVertex->vTangent,
                             pLightPathVertex->vShadingNormal,
                             pLightPathVertex->vBinormal, &vHalfVector);

  // 反射光線を計算
  float_t const hk =
      -hi::dot(vHalfVector, pLightPathVertex->vOutgoingDirection);
  pLightPathVertex->vOutgoingDirection += (2 * hk) * vHalfVector;

  float_t const nk2 = hi::dot(pLightPathVertex->vOutgoingDirection,
                              pLightPathVertex->vShadingNormal);
  if (nk2 <= 0) {
    return false;  // 打ち止め
  }

  float_t const fCosTheta_o = hi::dot(pLightPathVertex->vGeometricNormal,
                                      pLightPathVertex->vOutgoingDirection);
  if (fCosTheta_o <= 0) {
    return false;  // 打ち止め
  }

  // NOTE: ここで，エネルギーが増大することに留意する．
  float_t const fShadingFactor = FakeFresnelTerm(hk) / hi::max(nk1, nk2) * nk2 *
                                 (fCosTheta_o / fCosTheta_g);
  if (random.next<float_t>() >= albedo_[1] * fShadingFactor) {
    return false;  // 反射率で刈る
  }

  pLightPathVertex->bSpecular = false;

  return true;
}

float3_t Metal::CalculateWeightedDirectLighting(
    __inout PathVertex *const pEyePathVertex,
    __in PathVertex const &lightPathVertex, __in Scene const &scene) const {
  pEyePathVertex->power *= albedo_;

#ifdef CONFIG_EXPLICIT_GLOSSY
  if (!lightPathVertex.pGeometry) {
    return float3_t(0);
  }

  float3_t const vIncomingDirection =
      hi::normalize(pEyePathVertex->vPosition - lightPathVertex.vPosition);

  float_t const nk3 =
      hi::dot(vIncomingDirection, lightPathVertex.vGeometricNormal);
  if ((nk3 <= 0) ||
      (hi::dot(vIncomingDirection, lightPathVertex.vShadingNormal) <= 0)) {
    return float3_t(0);
  }

  float_t const nk2 =
      -hi::dot(vIncomingDirection, pEyePathVertex->vShadingNormal);
  if ((nk2 <= 0) ||
      (-hi::dot(vIncomingDirection, pEyePathVertex->vGeometricNormal) <= 0)) {
    return float3_t(0);
  }

  Intersection param;
  if (scene.FindIntersection(lightPathVertex.vPosition, vIncomingDirection,
                             &param) != pEyePathVertex->pGeometry) {
    return float3_t(0);
  }

  float_t const w = nk3 / hi::square_of(param.t_max());
  float_t const G = nk2 * w;
  float3_t const vHalfVector =
      -hi::normalize(vIncomingDirection + pEyePathVertex->vOutgoingDirection);
  float_t const hk = -hi::dot(vHalfVector, vIncomingDirection);
  float_t const nh = hi::dot(pEyePathVertex->vShadingNormal, vHalfVector);
  float_t const nk1 = -hi::dot(pEyePathVertex->vShadingNormal,
                               pEyePathVertex->vOutgoingDirection);
  float_t const fSamplingNext =
      (shininess_ + 1) * (M_1_PI / 8) * std::pow(nh, shininess_) / hk;
  float_t const fDensity = fSamplingNext * w * scene.GetLightArea();

  // 重みとして，光源上の点のサンプリング確率を除し，光源方向のサンプリング確率を掛ける
  return lightPathVertex.power * pEyePathVertex->power * (G * M_1_PI) *
         fSamplingNext / (std::max(nk1, nk2)) * FakeFresnelTerm(hk) *
         hi::rcp(1 + hi::square_of(fDensity));
#else
  return float3_t(0);
#endif
}

bool Metal::NextScatteringDirection(__inout PathVertex *const pEyePathVertex,
                                    __inout IPrimarySample &sample) const {
  float3_t vHalfVector;
  float_t const nh = ComputeGlossyHalfwayVector(
      sample.next(), sample.next(), shininess_, pEyePathVertex->vTangent,
      pEyePathVertex->vShadingNormal, pEyePathVertex->vBinormal, &vHalfVector);

  float_t const nk1 = -hi::dot(pEyePathVertex->vShadingNormal,
                               pEyePathVertex->vOutgoingDirection);
  float_t const hk = -hi::dot(vHalfVector, pEyePathVertex->vOutgoingDirection);
  pEyePathVertex->vOutgoingDirection += (2 * hk) * vHalfVector;  // 反射光線
  float_t const nk2 = hi::dot(pEyePathVertex->vShadingNormal,
                              pEyePathVertex->vOutgoingDirection);
  if ((nk2 <= 0) || (hi::dot(pEyePathVertex->vGeometricNormal,
                             pEyePathVertex->vOutgoingDirection) <= 0)) {
    return false;  // 打ち止め
  }

  pEyePathVertex->power *= FakeFresnelTerm(hk) / hi::max(nk1, nk2) * nk2;
#ifdef CONFIG_EXPLICIT_GLOSSY
  pEyePathVertex->bSpecular = false;
  pEyePathVertex->fSamplingPrev = 1;
  pEyePathVertex->fSamplingNext =
      (shininess_ + 1) * (M_1_PI / 8) * std::pow(nh, shininess_) / hk;
#else
  pEyePathVertex->bSpecular = true;
  pEyePathVertex->fSamplingPrev = 1;
  pEyePathVertex->fSamplingNext = 1;
#endif

  return true;
}

float_t Metal::GetDensityVariance() const {
  return std::sqrt(1 +
                   (hi::square_of(XYZ_CONFIG_kMaxBranchCount) - 1) /
                       hi::square_of(shininess_ + 1));
}

bool Metal::NextScatteringDirection(PathVertex &pathVertex, float3_t const &mu,
                                    bool const &bAdjoint) const {
  if (bAdjoint) {
    // ハーフベクトルを求める
    float3_t vHalfVector;
    ComputeGlossyHalfwayVector(mu[0], mu[1], shininess_, pathVertex.vTangent,
                               pathVertex.vGeometricNormal,
                               pathVertex.vBinormal, &vHalfVector);

    float_t const nk1 =
        -hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal);
    float_t const hk = -hi::dot(pathVertex.vOutgoingDirection, vHalfVector);
    pathVertex.vOutgoingDirection += (2 * hk) * vHalfVector;
    float_t const nk2 =
        hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal);
    if ((nk2 <= 0) || (hi::dot(pathVertex.vOutgoingDirection,
                               pathVertex.vShadingNormal) <= 0)) {
      return false;  // 打ち止め
    }

    pathVertex.power *=
        albedo_ * (FakeFresnelTerm(hk) / hi::max(nk1, nk2) * nk2);
  } else {
    // ハーフベクトルを求める
    float3_t vHalfVector;
    ComputeGlossyHalfwayVector(mu[0], mu[1], shininess_, pathVertex.vTangent,
                               pathVertex.vShadingNormal, pathVertex.vBinormal,
                               &vHalfVector);

    float_t const nk1 =
        -hi::dot(pathVertex.vOutgoingDirection, pathVertex.vShadingNormal);
    float_t const hk = -hi::dot(pathVertex.vOutgoingDirection, vHalfVector);
    pathVertex.vOutgoingDirection += (2 * hk) * vHalfVector;
    float_t const nk2 =
        hi::dot(pathVertex.vOutgoingDirection, pathVertex.vShadingNormal);
    if ((nk2 <= 0) || (hi::dot(pathVertex.vOutgoingDirection,
                               pathVertex.vGeometricNormal) <= 0)) {
      return false;  // 打ち止め
    }

    pathVertex.power *=
        albedo_ * (FakeFresnelTerm(hk) / hi::max(nk1, nk2) * nk2);
  }

  pathVertex.bSpecular = true;
  pathVertex.fSamplingPrev = 1;
  pathVertex.fSamplingNext = 1;

  return true;
}
}
// end of namespace xyz

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
bool Metal::LightsToCameraScatteringDirection(
    __inout PathVertex *const pPathVertex,
    __inout IPrimarySample &sample) const {
  // デルタ反射(or屈折)を含まないのでfalse
  pPathVertex->bSpecular = false;

  // ハーフベクトルをサンプリングする．
  float3_t vHalfVector;
  float_t const nh = ComputeGlossyHalfwayVector(
      sample.next(), sample.next(), shininess_, pPathVertex->vTangent,
      pPathVertex->vShadingNormal, pPathVertex->vBinormal, &vHalfVector);

  // ハーフベクトルと(入射＆散乱)方向ベクトルとの内積を求める．
  float_t const hk = hi::dot(vHalfVector, pPathVertex->vIncomingDirection);

  // ハーフベクトルと散乱方向ベクトルから入射方向ベクトルを求める．
  pPathVertex->vOutgoingDirection =
      (2 * hk) * vHalfVector - pPathVertex->vIncomingDirection;

  // 幾何学的法線と散乱方向ベクトルの内積を求める．
  pPathVertex->fOutgoingCosThetaGeometric =
      hi::dot(pPathVertex->vGeometricNormal, pPathVertex->vOutgoingDirection);

  if (pPathVertex->fOutgoingCosThetaGeometric <= 0) {
    return false;  // 幾何学的法線の裏方向に入射方向ベクトルが向いたら打ち止める．
  }

  // シェーディング法線と散乱方向ベクトルの内積を求める．
  float_t const fOutgoingCosThetaShading =
      hi::dot(pPathVertex->vShadingNormal, pPathVertex->vOutgoingDirection);
  if (fOutgoingCosThetaShading <= 0) {
    return false;
  }

  // 幾何学的法線と入射方向ベクトルの内積を求める．
  float_t const fIncomingCosThetaGeometric =
      hi::dot(pPathVertex->vGeometricNormal, pPathVertex->vIncomingDirection);

  pPathVertex->fBSDFxIPDF =
      FakeFresnelTerm(hk) * pPathVertex->fIncomingCosThetaShading /
      hi::max(fOutgoingCosThetaShading, pPathVertex->fIncomingCosThetaShading) *
      pPathVertex->fOutgoingCosThetaGeometric / fIncomingCosThetaGeometric;

  pPathVertex->fSamplingPrev = pPathVertex->fSamplingNext =
      (shininess_ + 1) * (M_1_PI / 8) * std::pow(nh, shininess_) / hk;

  pPathVertex->fSamplingPrev *= hi::rcp(pPathVertex->fIncomingCosThetaShading);
  pPathVertex->fSamplingNext *=
      hi::rcp(pPathVertex->fOutgoingCosThetaGeometric);

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
bool Metal::CameraToLightsScatteringDirection(
    __inout PathVertex *const pPathVertex,
    __inout IPrimarySample &sample) const {
  // デルタ反射(or屈折)を含まないのでfalse
  pPathVertex->bSpecular = false;

  // ハーフベクトルをサンプリングする．
  float3_t vHalfVector;
  float_t const nh = ComputeGlossyHalfwayVector(
      sample.next(), sample.next(), shininess_, pPathVertex->vTangent,
      pPathVertex->vShadingNormal, pPathVertex->vBinormal, &vHalfVector);

  // ハーフベクトルと(入射＆散乱)方向ベクトルとの内積を求める．
  float_t const hk = hi::dot(vHalfVector, pPathVertex->vOutgoingDirection);

  // ハーフベクトルと散乱方向ベクトルから入射方向ベクトルを求める．
  pPathVertex->vIncomingDirection =
      (2 * hk) * vHalfVector - pPathVertex->vOutgoingDirection;

  float_t const fIncomingCosThetaGeometric =
      hi::dot(pPathVertex->vGeometricNormal, pPathVertex->vIncomingDirection);
  if (fIncomingCosThetaGeometric <= 0) {
    return false;  // 幾何学的法線の裏方向に入射方向ベクトルが向いたら打ち止める．
  }

  // シェーディング法線と入射方向ベクトルの内積を求める．
  pPathVertex->fIncomingCosThetaShading =
      hi::dot(pPathVertex->vShadingNormal, pPathVertex->vIncomingDirection);
  if (pPathVertex->fIncomingCosThetaShading <= 0) {
    return false;
  }

  // シェーディング法線と散乱方向ベクトルの内積を求める．
  float_t const fOutgoingCosThetaShading =
      hi::dot(pPathVertex->vShadingNormal, pPathVertex->vOutgoingDirection);

  pPathVertex->fBSDFxIPDF =
      FakeFresnelTerm(hk) * pPathVertex->fIncomingCosThetaShading /
      hi::max(fOutgoingCosThetaShading, pPathVertex->fIncomingCosThetaShading);

  pPathVertex->fSamplingPrev = pPathVertex->fSamplingNext =
      (shininess_ + 1) * (M_1_PI / 8) * std::pow(nh, shininess_) / hk;

  pPathVertex->fSamplingPrev *= hi::rcp(pPathVertex->fIncomingCosThetaShading);
  pPathVertex->fSamplingNext *=
      hi::rcp(pPathVertex->fOutgoingCosThetaGeometric);

  return true;
}
}  // end of namespace xyz

/*
視点から光源方向へのサンプリング(fSamplingPrev)では，
シェーディング法線と入射方向ベクトルの内積(fIncomingCosThetaShading)で，
サンプリング分布を除す

光源から視点方向へのサンプリング(fSamplingNext)では，
幾何学的法線と散乱方向ベクトルの内積(fOutgoingCosThetaGeometric)で，
サンプリング分布を除す
*/

/*
TODO:
fOutgoingCosThetaShading
fIncomingCosThetaShading

fOutgoingCosThetaGeometric
fIncomingCosThetaGeometric

をそれぞれ経路頂点に持っておくようにする．
*/

namespace xyz {
/*
  virtual float_t Metal::ProjectedPDF(
    float3_t const & wi, ///< [in] incoming direction
    float3_t const & wo, ///< [in] outgoing direction
    float3_t const & ng, ///< [in] geometric normal
    float3_t const & ns  ///< [in] shading normal
 ) const
  {
    // p(w) = p'(w) / |cos\theta|
    //   p'(w) = ph(h) / dot(w,h) / 4
    //   ph(h) = (s+1) / 2pi * dot(n,h)^{s}
    float3_t const h = hi::normalize(wi + wo);
    float_t const ph = std::pow(hi::dot(n,h), shininess_);
    return (M_1_PI/8) * (shininess_+1) * ph / hi::dot(n, wi);
  }
*/
float_t Metal::BSDF(std::size_t const &is,  ///< [in] index of spectrum
                    float3_t const &wi,     ///< [in] incoming direction
                    float3_t const &wo,     ///< [in] outgoing direction
                    float3_t const &ns      ///< [in] shading normal
                    ) const {
  // f_{s}(w_{i} -> w_{o})
  float3_t const h = hi::normalize(wi + wo);
  float_t const hw = hi::dot(h, wi);
  float_t const hn = hi::dot(h, ns);
  return albedo_[is] * (M_1_PI / 8) * (shininess_ + 1) *
         std::pow(hn, shininess_) /
         (hw * std::max(hi::dot(ns, wi), hi::dot(ns, wo))) *
         FakeFresnelTerm(hw);
}
}
// end of namespace xyz
