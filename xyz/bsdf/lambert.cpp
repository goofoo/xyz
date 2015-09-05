#include "lambert.hpp"
#include "../core/scene.hpp"
#include "../geom/pathvertex.hpp"
#include "../geom/intersection.hpp"

namespace xyz {
void Lambert::Render() const {
  hi::basic_vector3<float> rgb;
  hi::CCIR601_1_XYZ2RGB(albedo_, rgb);

  ::glEnable(GL_LIGHTING);
  ::glDisable(GL_BLEND);
  hi::glMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, rgb.data());
}

bool Lambert::BoundsImportance(PathVertex *const pLightPathVertex,
                               std::mt19937_64 &random) const {
  // 光源からの追跡を想定している(Adjoint)
  // vOutgoingDirection には 入射してきた方向が入っている
  float_t const fCosTheta_s = -hi::dot(pLightPathVertex->vOutgoingDirection,
                                       pLightPathVertex->vShadingNormal);
  if (fCosTheta_s <= 0) {
    return false;  // シェーディング法線の裏に当たったら，素直に終了する
  }
  float_t const fCosTheta_g = -hi::dot(pLightPathVertex->vOutgoingDirection,
                                       pLightPathVertex->vShadingNormal);

  /// NOTE: ここで，エネルギーが増大する可能性に留意する
  float_t const fShadingFactor = fCosTheta_s / fCosTheta_g;
  if (random.next<float_t>() >= albedo_[1] * fShadingFactor) {
    return false;  // 反射率で刈る
  }

  pLightPathVertex->bSpecular = false;

  // 幾何学的法線に基づいて，反射方向を決める
  ComputeDiffusedVector(
      random.next<float_t>(), random.next<float_t>(),
      pLightPathVertex->vTangent, pLightPathVertex->vGeometricNormal,
      pLightPathVertex->vBinormal, &pLightPathVertex->vOutgoingDirection);

  return true;
  /*
      if (random.next<float_t>() >= albedo_[1])
      {
        return false; // 反射率で刈る
      }

      eyePathVertex.bSpecular = false;

  #if 0
      // こっちは視点からの重要度を決めるとき
      ComputeDiffusedVector(random.next<float_t>(), random.next<float_t>(),
        eyePathVertex.vTangent, eyePathVertex.vShadingNormal,
  eyePathVertex.vBinormal,
        eyePathVertex.vOutgoingDirection);

      return hi::dot(eyePathVertex.vOutgoingDirection,
  eyePathVertex.vGeometricNormal) > 0;
  #else
      // 光源からの重要度を計算するので幾何法線で反射方向を決める
      ComputeDiffusedVector(random.next<float_t>(), random.next<float_t>(),
        eyePathVertex.vTangent, eyePathVertex.vGeometricNormal,
  eyePathVertex.vBinormal,
        &eyePathVertex.vOutgoingDirection);

      return true;
  #endif
  */
}

bool Lambert::BoundsPhoton(PathVertex &lightPathVertex,
                           std::mt19937_64 &random) const {
  if (random.next<float_t>() >= CONFIG_ALBEDO) {
    return false;
  }

  lightPathVertex.bSpecular = false;

#if 0
    ComputeDiffusedVector(random.next<float_t>(), random.next<float_t>(),
      lightPathVertex.vTangent, lightPathVertex.vGeometricNormal, lightPathVertex.vBinormal,
      lightPathVertex.vOutgoingDirection);

    return hi::dot(lightPathVertex.vOutgoingDirection, lightPathVertex.vShadingNormal) > 0;
#else
  ComputeDiffusedVector(
      random.next<float_t>(), random.next<float_t>(), lightPathVertex.vTangent,
      lightPathVertex.vGeometricNormal, lightPathVertex.vBinormal,
      &lightPathVertex.vOutgoingDirection);

  return true;
#endif
}

bool Lambert::BoundsPhotonWithImportonMap(PathVertex &lightPathVertex,
                                          ImportonMap const &importonMap,
                                          ImportonQuery &importons,
                                          std::mt19937_64 &random) const {
  if (random.next<float_t>() >= CONFIG_ALBEDO) {
    return false;
  }

  lightPathVertex.bSpecular = false;
  importonMap.kNN_query(lightPathVertex.vPosition, importons);
  importons.SetDiffusedVector(random);

  return true;
}

bool Lambert::BoundsPhotonWithParticleFilter(PathVertex &lightPathVertex,
                                             ImportonMap const &importonMap,
                                             ParticleQuery &importons,
                                             std::mt19937_64 &random) const {
  if (random.next<float_t>() >= CONFIG_ALBEDO) {
    return false;
  }

  lightPathVertex.bSpecular = false;
  importonMap.kNN_query(lightPathVertex.vPosition, importons);

  return true;
}

float3_t Lambert::CalculateWeightedDirectLighting(
    __inout PathVertex *const pEyePathVertex,
    __in PathVertex const &lightPathVertex, __in Scene const &scene) const {
  pEyePathVertex->power *= albedo_;

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

  // 幾何学項
  float_t const G = nk2 * nk3 / hi::square_of(param.t_max());

  // 暗黙的サンプリングの経路密度
  float_t const fDensityOfImplicit = M_1_PI * G;

  // 明示的サンプリングの経路密度
  float_t const fDensityOfExplicit = hi::rcp(scene.GetLightArea());

  // 経路密度の比
  float_t const fDensity = fDensityOfImplicit / fDensityOfExplicit;

  // 重みとして，光源上の点のサンプリング確率を除し，光源方向のサンプリング確率を掛ける
  return lightPathVertex.power * pEyePathVertex->power *
         ((G * (M_1_PI * M_1_PI)) *
          hi::rcp(1 + hi::square_of(fDensity)));  // 重み
}

// BRDFのPDF(シェーディング法線)に基づいて反射方向をサンプリングする．
bool Lambert::NextScatteringDirection(__inout PathVertex *const pEyePathVertex,
                                      __inout IPrimarySample &sample) const {
  pEyePathVertex->bSpecular = false;
  pEyePathVertex->fSamplingPrev = 1;
  pEyePathVertex->fSamplingNext =
      M_1_PI * ComputeDiffusedVector(
                   sample.next(), sample.next(), pEyePathVertex->vTangent,
                   pEyePathVertex->vShadingNormal, pEyePathVertex->vBinormal,
                   &pEyePathVertex->vOutgoingDirection);
  return hi::dot(pEyePathVertex->vOutgoingDirection,
                 pEyePathVertex->vGeometricNormal) > 0;
}

float3_t Lambert::CalculateWeightedDirectLighting(
    __inout PathVertex *const pEyePathVertex,
    __in PathVertex const &lightPathVertex,
    __in ImportanceQuery const &importons, __in Scene const &scene) const {
  pEyePathVertex->power *= albedo_;

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

  // 幾何学項
  float_t const G = nk2 * nk3 / hi::square_of(param.t_max());

  // 暗黙的サンプリングの経路密度
  float_t const fDensityOfImplicit =
      M_1_PI * G * importons.GetDensity(vIncomingDirection);

  // 明示的サンプリングの経路密度
  float_t const fDensityOfExplicit = hi::rcp(scene.GetLightArea());

  // 経路密度の比
  float_t const fDensity = fDensityOfImplicit / fDensityOfExplicit;

  // 重みとして，光源上の点のサンプリング確率を除し，光源方向のサンプリング確率を掛ける
  return lightPathVertex.power * pEyePathVertex->power *
         ((G * (M_1_PI * M_1_PI)) *
          hi::rcp(1 + hi::square_of(fDensity)));  // 重み
}

// BRDFのPDF(シェーディング法線)に基づいて反射方向をサンプリングする．
// ただし Photon Map に基づいて光源方向についての情報を同時に利用する．
bool Lambert::NextScatteringDirection(__inout PathVertex *const pEyePathVertex,
                                      __in ImportanceQuery const &importons,
                                      __inout IPrimarySample &sample) const {
  pEyePathVertex->bSpecular = false;

  float_t fIncomingCosThetaShading;
  float_t fAlpha;  // 1/PDF
  importons.SetDiffusedVector(sample, &(pEyePathVertex->vOutgoingDirection),
                              &fIncomingCosThetaShading,
                              &(pEyePathVertex->fSamplingNext), &fAlpha);
  pEyePathVertex->power *= fAlpha;  // パワーをスケーリング．

  return hi::dot(pEyePathVertex->vOutgoingDirection,
                 pEyePathVertex->vGeometricNormal) > 0;
}

float_t Lambert::GetDensityVariance() const {
  return XYZ_CONFIG_kMaxBranchCount;
}

/// 双方向経路追跡法用の反射光線生成器
bool Lambert::NextScatteringDirection(PathVertex &pathVertex,  ///< [in/out]
                                      float3_t const &mu,      ///< [in]
                                      bool const &bAdjoint     ///< [in]
                                      ) const {
#ifdef CONFIG_RUSSIAN_RULLET
  float_t const t = pathVertex.power;
  if (t < CONFIG_TERMINATE_THRESHOLD) {
    float_t const q = t / CONFIG_TERMINATE_THRESHOLD;
    if (mu[2] >= q) {
      return false;  // 打ち止め
    }
    rMu *= hi::rcp(q);
  }
#endif

  pathVertex.bSpecular = false;
  pathVertex.fSamplingPrev = M_1_PI;
  pathVertex.fSamplingNext = M_1_PI;

  // next direction and importance (potential) or power
  if (bAdjoint) {
    ComputeDiffusedVector(mu[0], mu[1], pathVertex.vTangent,
                          pathVertex.vGeometricNormal, pathVertex.vBinormal,
                          &pathVertex.vOutgoingDirection);
    return hi::dot(pathVertex.vOutgoingDirection, pathVertex.vShadingNormal) >
           0;
  } else {
    ComputeDiffusedVector(mu[0], mu[1], pathVertex.vTangent,
                          pathVertex.vShadingNormal, pathVertex.vBinormal,
                          &pathVertex.vOutgoingDirection);
    return hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal) >
           0;
  }
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
bool Lambert::LightsToCameraScatteringDirection(
    __inout PathVertex *const pPathVertex,
    __inout IPrimarySample &sample) const {
  // 幾何学的法線と入射方向ベクトルの内積を求める．
  float_t const fIncomingCosThetaGeometric =
      hi::dot(pPathVertex->vIncomingDirection, pPathVertex->vGeometricNormal);

  // BRDF/PDF を計算する．
  float_t const fAlpha =
      pPathVertex->fIncomingCosThetaShading / fIncomingCosThetaGeometric;

  // 延長確率を計算する．
  float_t const fProbPrev = std::min(float_t(1), albedo_[1]);
  float_t const fProbNext = std::min(float_t(1), albedo_[1] * fAlpha);

  // 完全拡散反射(ランバート)面なので鏡面反射はしない．
  pPathVertex->bSpecular = false;

  pPathVertex->fSamplingPrev =
      M_1_PI * fProbPrev;  // 視点から光源方向のサンプリング密度
  pPathVertex->fSamplingNext =
      M_1_PI * fProbNext;  // 光源から視点方向のサンプリング密度

  // BRDF/PDF で重み付けする．
  pPathVertex->power *= albedo_;
  pPathVertex->power *=
      fAlpha;  // TODO: ランバート面以外に対応する場合は修正が必要

  // TODO: ランバート面しか扱わないので現状はこんな感じ．
  pPathVertex->fBSDFxIPDF = hi::rcp(fProbNext);

  // 延長確率で刈る．
  if (sample.next() >= fProbNext) {
    return false;
  }

  // 幾何学的法線に基づいて，散乱方向ベクトルを決める．
  pPathVertex->fOutgoingCosThetaGeometric = ComputeDiffusedVector(
      sample.next(), sample.next(), pPathVertex->vTangent,
      pPathVertex->vGeometricNormal, pPathVertex->vBinormal,
      &(pPathVertex->vOutgoingDirection));
  float_t const fOutgoingCosThetaShaing =
      hi::dot(pPathVertex->vOutgoingDirection, pPathVertex->vShadingNormal);
  if (fOutgoingCosThetaShaing <= 0) {
    return false;
  }

  // シェーディング法線と入射方向ベクトルの内積を求める．
  // float_t const fIncomingCosThetaShading
  //  = hi::dot(pPathVertex->vIncomingDirection, pPathVertex->vShadingNormal);

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
bool Lambert::CameraToLightsScatteringDirection(
    __inout PathVertex *const pPathVertex,
    __inout IPrimarySample &sample) const {
  // 完全拡散反射(ランバート)面なので鏡面反射はしない．
  pPathVertex->bSpecular = false;

  // シェーディング法線をベースに入射方向ベクトルを求め
  // シェーディング法線と入射方向ベクトルの内積を得る．
  pPathVertex->fIncomingCosThetaShading =
      ComputeDiffusedVector(sample.next(), sample.next(), pPathVertex->vTangent,
                            pPathVertex->vShadingNormal, pPathVertex->vBinormal,
                            &(pPathVertex->vIncomingDirection));
  float_t const fIncomingCosThetaGeometric =
      hi::dot(pPathVertex->vIncomingDirection, pPathVertex->vGeometricNormal);
  if (fIncomingCosThetaGeometric <= 0) {
    return false;
  }

  // 幾何学的法線と散乱方向ベクトルの内積を求める．
  float_t const fOutgoingCosThetaGeometric =
      hi::dot(pPathVertex->vOutgoingDirection, pPathVertex->vGeometricNormal);

  // BRDF/PDF を計算する．
  float_t const fAlpha =
      pPathVertex->fIncomingCosThetaShading / fOutgoingCosThetaGeometric;
  float_t const pProbPrev = std::min(float_t(1), albedo_[1]);
  float_t const pProbNext = std::min(float_t(1), albedo_[1] * fAlpha);

  pPathVertex->fSamplingPrev =
      M_1_PI * pProbPrev;  // 視点から光源方向のサンプリング密度
  pPathVertex->fSamplingNext =
      M_1_PI * pProbNext;  // 光源から視点方向のサンプリング密度

  // BRDF/PDF で重み付けする．
  pPathVertex->power *= albedo_;

  // TODO: ランバート面しか扱わないので現状はこんな感じ．
  pPathVertex->fBSDFxIPDF = hi::rcp(pProbPrev);

  if (sample.next() >= pProbPrev) {
    return false;
  }

  return true;
}

bool Lambert::CameraToLightsScatteringDirection(
    __inout PathVertex *const pPathVertex,
    __in ImportanceQuery const &importons,
    __inout IPrimarySample &sample) const {
  // 完全拡散反射(ランバート)面なので鏡面反射はしない．
  pPathVertex->bSpecular = false;

  // シェーディング法線をベースに入射方向ベクトルを求め
  // シェーディング法線と入射方向ベクトルの内積を得る．
  float_t fEyeToLightAlpha;  // 1/PDF
  importons.SetDiffusedVector(
      sample, &(pPathVertex->vIncomingDirection),
      &(pPathVertex->fIncomingCosThetaShading),
      &(pPathVertex->fSamplingNext),  // NOTE: ここで設定した fSamplingNext
      // の内容は現在は使用していない．
      &fEyeToLightAlpha);

  float_t const fIncomingCosThetaGeometric =
      hi::dot(pPathVertex->vIncomingDirection, pPathVertex->vGeometricNormal);
  if (fIncomingCosThetaGeometric <= 0) {
    return false;
  }

  // 幾何学的法線と散乱方向ベクトルの内積を求める．
  float_t const fOutgoingCosThetaGeometric =
      hi::dot(pPathVertex->vOutgoingDirection, pPathVertex->vGeometricNormal);

  // BRDF/PDF を計算する．
  float_t const fLightToEyeAlpha =
      pPathVertex->fIncomingCosThetaShading / fOutgoingCosThetaGeometric;
  float_t const pProbPrev = std::min(float_t(1), albedo_[1]);
  float_t const pProbNext = std::min(float_t(1), albedo_[1] * fLightToEyeAlpha);

  pPathVertex->fSamplingPrev =
      M_1_PI *
      pProbPrev;  // 視点から光源方向のサンプリング密度(TODO: 重みは未修正)
  pPathVertex->fSamplingNext =
      M_1_PI * pProbNext;  // 光源から視点方向のサンプリング密度

  // TODO: ランバート面しか扱わないので現状はこんな感じ．
  pPathVertex->power *= albedo_;

  // BRDF/PDF で重み付けする．
  pPathVertex->fBSDFxIPDF =
      hi::rcp(pProbPrev) * fEyeToLightAlpha;  // fEyeToLightAlpha は追加要素

  if (sample.next() >= pProbPrev) {
    return false;
  }

  return true;
}
}

/*
視点から光源方向へのサンプリング(fSamplingPrev)では，
シェーディング法線と入射方向ベクトルの内積(fIncomingCosThetaShading)で，
サンプリング分布を除す

光源から視点方向へのサンプリング(fSamplingNext)では，
幾何学的法線と散乱方向ベクトルの内積(fOutgoingCosThetaGeometric)で，
サンプリング分布を除す
*/

namespace xyz {
/*
  virtual float_t Lambert::ProjectedPDF(
    float3_t const & wi, ///< [in] incoming direction
    float3_t const & wo, ///< [in] outgoing direction
    float3_t const & ng, ///< [in] geometric normal
    float3_t const & ns  ///< [in] shading normal
 ) const
  {
    // p(w) = p'(w) / |cos\theta|
    //   p'(w) = |cos\theta| / pi
    return M_1_PI;
  }
*/
float_t Lambert::BSDF(std::size_t const &is,  ///< [in] index of spectrum
                      float3_t const &wi,     ///< [in] incoming direction
                      float3_t const &wo,     ///< [in] outgoing direction
                      float3_t const &ns      ///< [in] shading normal
                      ) const {
  UNREFERENCED_PARAMETER(is);
  UNREFERENCED_PARAMETER(wi);
  UNREFERENCED_PARAMETER(wo);
  UNREFERENCED_PARAMETER(ns);

  // f_{s}(w_{i} -> w_{o})
  return albedo_[is] * M_1_PI;
}
}
// end of namespace xyz
