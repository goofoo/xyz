#include "Lambert.hpp"
#include "core/Scene.hpp"
#include "geom/PathVertex.hpp"
#include "geom/Intersection.hpp"
#include "sampler/PrimarySample.hpp"

namespace tgir {
void Lambert::Render() const {
  tgir::SpectrumVector xyz;
  hi::spd2xyz(albedo_, xyz);

  hi::basic_vector3<float> rgb;
  hi::CCIR601_1_XYZ2RGB(xyz, rgb);

  ::glEnable(GL_LIGHTING);
  ::glDisable(GL_BLEND);
  hi::glMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, rgb.data());
}

bool Lambert::BoundsImportance(tgir::PathVertex &eyePathVertex,
                               std::mt19937_64 &random) const {
  if (std::uniform_real_distribution<tgir::Real>()(random) >= CONFIG_ALBEDO) {
    return false;  // 反射率で刈る
  }

  eyePathVertex.bSpecular = false;

#if 0
    tgir::ComputeDiffusedVector(std::uniform_real_distribution<tgir::Real>()(random), std::uniform_real_distribution<tgir::Real>()(random),
      eyePathVertex.vTangent, eyePathVertex.vShadingNormal, eyePathVertex.vBinormal,
      eyePathVertex.vOutgoingDirection);

    return hi::dot(eyePathVertex.vOutgoingDirection, eyePathVertex.vGeometricNormal) > 0;
#else
  tgir::ComputeDiffusedVector(
      std::uniform_real_distribution<tgir::Real>()(random),
      std::uniform_real_distribution<tgir::Real>()(random),
      eyePathVertex.vTangent, eyePathVertex.vGeometricNormal,
      eyePathVertex.vBinormal, &eyePathVertex.vOutgoingDirection);

  return true;
#endif
}

bool Lambert::BoundsPhoton(tgir::PathVertex &lightPathVertex,
                           std::mt19937_64 &random) const {
  if (std::uniform_real_distribution<tgir::Real>()(random) >= CONFIG_ALBEDO) {
    return false;
  }

  lightPathVertex.bSpecular = false;

#if 0
    tgir::ComputeDiffusedVector(std::uniform_real_distribution<tgir::Real>()(random), std::uniform_real_distribution<tgir::Real>()(random),
      lightPathVertex.vTangent, lightPathVertex.vGeometricNormal, lightPathVertex.vBinormal,
      lightPathVertex.vOutgoingDirection);

    return hi::dot(lightPathVertex.vOutgoingDirection, lightPathVertex.vShadingNormal) > 0;
#else
  tgir::ComputeDiffusedVector(
      std::uniform_real_distribution<tgir::Real>()(random),
      std::uniform_real_distribution<tgir::Real>()(random),
      lightPathVertex.vTangent, lightPathVertex.vGeometricNormal,
      lightPathVertex.vBinormal, &lightPathVertex.vOutgoingDirection);

  return true;
#endif
}

bool Lambert::BoundsPhotonWithImportonMap(tgir::PathVertex &lightPathVertex,
                                          tgir::ImportonMap const &importonMap,
                                          tgir::ImportonQuery &importons,
                                          std::mt19937_64 &random) const {
  if (std::uniform_real_distribution<tgir::Real>()(random) >= CONFIG_ALBEDO) {
    return false;
  }

  lightPathVertex.bSpecular = false;
  importonMap.kNN_query(lightPathVertex.vPosition, importons);
  importons.SetDiffusedVector(random);

  return true;
}

bool Lambert::BoundsPhotonWithParticleFilter(
    tgir::PathVertex &lightPathVertex, tgir::ImportonMap const &importonMap,
    tgir::ParticleQuery &importons, std::mt19937_64 &random) const {
  if (std::uniform_real_distribution<tgir::Real>()(random) >= CONFIG_ALBEDO) {
    return false;
  }

  lightPathVertex.bSpecular = false;
  importonMap.kNN_query(lightPathVertex.vPosition, importons);

  return true;
}

tgir::Spectrum Lambert::CalculateWeightedDirectLighting(
    tgir::PathVertex &eyePathVertex, tgir::PathVertex const &lightPathVertex,
    std::size_t const uWavelengthIndex, tgir::Scene const &scene) const {
  eyePathVertex.sQuantum *= albedo_[uWavelengthIndex];

  if (!lightPathVertex.pGeometry) {
    return 0;
  }

  tgir::Vector3 const vIncomingDirection =
      hi::normalize(eyePathVertex.vPosition - lightPathVertex.vPosition);

  tgir::Real const nk3 =
      hi::dot(vIncomingDirection, lightPathVertex.vGeometricNormal);
  if ((nk3 <= 0) ||
      (hi::dot(vIncomingDirection, lightPathVertex.vShadingNormal) <= 0)) {
    return 0;
  }

  tgir::Real const nk2 =
      -hi::dot(vIncomingDirection, eyePathVertex.vShadingNormal);
  if ((nk2 <= 0) ||
      (-hi::dot(vIncomingDirection, eyePathVertex.vGeometricNormal) <= 0)) {
    return 0;
  }

  tgir::Intersection param;
  if (scene.FindIntersection(lightPathVertex.vPosition, vIncomingDirection,
                             &param) != eyePathVertex.pGeometry) {
    return 0;
  }

  tgir::Real const G = nk2 * nk3 / hi::square_of(param.t_max());  // 幾何学項
  tgir::Real const rDensity = M_1_PI * G * scene.GetLightArea();

  // 重みとして，光源上の点のサンプリング確率を除し，光源方向のサンプリング確率を掛ける
  return lightPathVertex.sQuantum * eyePathVertex.sQuantum *
         (G * (M_1_PI * M_1_PI)) *
         hi::rcp(1 + hi::square_of(rDensity));  // 重み
}

bool Lambert::GetScatterVector(tgir::PathVertex &eyePathVertex,
                               tgir::Real const,
                               std::mt19937_64 &random) const {
  eyePathVertex.bSpecular = false;
  eyePathVertex.rSamplingPrev = 1;
  eyePathVertex.rSamplingNext =
      M_1_PI * tgir::ComputeDiffusedVector(
                   std::uniform_real_distribution<tgir::Real>()(random),
                   std::uniform_real_distribution<tgir::Real>()(random),
                   eyePathVertex.vTangent, eyePathVertex.vShadingNormal,
                   eyePathVertex.vBinormal, &eyePathVertex.vOutgoingDirection);
  return hi::dot(eyePathVertex.vOutgoingDirection,
                 eyePathVertex.vGeometricNormal) > 0;
}

bool Lambert::GetScatterVector(tgir::PathVertex &eyePathVertex,
                               tgir::Real const,
                               tgir::ImportonMap const &importonMap,
                               tgir::ImportanceQuery &importons,
                               std::mt19937_64 &random) const {
  eyePathVertex.bSpecular = false;

  importonMap.kNN_query(eyePathVertex.vPosition, importons);
  importons.SetDiffusedVector(random);

  return hi::dot(eyePathVertex.vOutgoingDirection,
                 eyePathVertex.vGeometricNormal) > 0;
}

tgir::Real Lambert::GetDensityVariance() const {
  return TGIR_CONFIG_kMaxBranchCount;
}

/// 双方向経路追跡法用の反射光線生成器
bool Lambert::GetScatterVector(tgir::PathVertex &pathVertex,  ///< [in/out]
                               tgir::PrimarySample const &,   ///< [in]
                               tgir::Vector3 const &mu,       ///< [in]
                               bool const &bAdjoint           ///< [in]
                               ) const {
#ifdef CONFIG_RUSSIAN_RULLET
  tgir::Real const t = pathVertex.sQuantum;
  if (t < CONFIG_TERMINATE_THRESHOLD) {
    tgir::Real const q = t / CONFIG_TERMINATE_THRESHOLD;
    if (mu[2] >= q) {
      return false;  // 打ち止め
    }
    rMu *= hi::rcp(q);
  }
#endif

  pathVertex.bSpecular = false;
  pathVertex.rSamplingPrev = M_1_PI;
  pathVertex.rSamplingNext = M_1_PI;

  // next direction and importance (potential) or power
  if (bAdjoint) {
    tgir::ComputeDiffusedVector(
        mu[0], mu[1], pathVertex.vTangent, pathVertex.vGeometricNormal,
        pathVertex.vBinormal, &pathVertex.vOutgoingDirection);
    return hi::dot(pathVertex.vOutgoingDirection, pathVertex.vShadingNormal) >
           0;
  } else {
    tgir::ComputeDiffusedVector(mu[0], mu[1], pathVertex.vTangent,
                                pathVertex.vShadingNormal, pathVertex.vBinormal,
                                &pathVertex.vOutgoingDirection);
    return hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal) >
           0;
  }
}

}  // end of namespace tgir

namespace tgir {
/*
  virtual tgir::Real Lambert::ProjectedPDF(
    tgir::Vector3 const & wi, ///< [in] incoming direction
    tgir::Vector3 const & wo, ///< [in] outgoing direction
    tgir::Vector3 const & ng, ///< [in] geometric normal
    tgir::Vector3 const & ns  ///< [in] shading normal
 ) const
  {
    // p(w) = p'(w) / |cos\theta|
    //   p'(w) = |cos\theta| / pi
    return M_1_PI;
  }
*/
tgir::Real Lambert::BSDF(std::size_t const &is,    ///< [in] index of spectrum
                         tgir::Vector3 const &wi,  ///< [in] incoming direction
                         tgir::Vector3 const &wo,  ///< [in] outgoing direction
                         tgir::Vector3 const &ns   ///< [in] shading normal
                         ) const {
  UNREFERENCED_PARAMETER(is);
  UNREFERENCED_PARAMETER(wi);
  UNREFERENCED_PARAMETER(wo);
  UNREFERENCED_PARAMETER(ns);

  // f_{s}(w_{i} -> w_{o})
  return albedo_[is] * M_1_PI;
}
}
// end of namespace tgir
