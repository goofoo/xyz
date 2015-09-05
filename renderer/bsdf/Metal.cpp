#include "Metal.hpp"
#include "core/Scene.hpp"
#include "geom/PathVertex.hpp"
#include "geom/Intersection.hpp"
#include "sampler/PrimarySample.hpp"

namespace tgir {
//
// 金属
//
void Metal::Render() const {
  tgir::SpectrumVector xyz;
  hi::spd2xyz(albedo_, xyz);

  hi::basic_vector3<float> rgb;
  hi::CCIR601_1_XYZ2RGB(xyz, rgb);

  ::glEnable(GL_LIGHTING);
  ::glDisable(GL_BLEND);
  hi::glMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, rgb.data());
  hi::glMaterial(GL_FRONT_AND_BACK, GL_SPECULAR, rgb.data());
  hi::glMaterial(GL_FRONT_AND_BACK, GL_SHININESS,
                 static_cast<GLfloat>(shininess_));
}

tgir::Spectrum Metal::CalculateWeightedDirectLighting(
    tgir::PathVertex &eyePathVertex, tgir::PathVertex const &lightPathVertex,
    std::size_t const uWavelengthIndex, tgir::Scene const &scene) const {
  eyePathVertex.sQuantum *= albedo_[uWavelengthIndex];

#ifdef CONFIG_EXPLICIT_GLOSSY
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

  tgir::Real const w = nk3 / hi::square_of(param.t_max());
  tgir::Real const G = nk2 * w;
  tgir::Vector3 const vHalfVector =
      -hi::normalize(vIncomingDirection + eyePathVertex.vOutgoingDirection);
  tgir::Real const hk = -hi::dot(vHalfVector, vIncomingDirection);
  tgir::Real const nh = hi::dot(eyePathVertex.vShadingNormal, vHalfVector);
  tgir::Real const nk1 =
      -hi::dot(eyePathVertex.vShadingNormal, eyePathVertex.vOutgoingDirection);
  tgir::Real const rSamplingNext =
      (shininess_ + 1) * (M_1_PI / 8) * std::pow(nh, shininess_) / hk;
  tgir::Real const rDensity = rSamplingNext * w * scene.GetLightArea();

  // 重みとして，光源上の点のサンプリング確率を除し，光源方向のサンプリング確率を掛ける
  return lightPathVertex.sQuantum * eyePathVertex.sQuantum * (G * M_1_PI) *
         rSamplingNext / (std::max(nk1, nk2)) * tgir::FakeFresnelTerm(hk) *
         hi::rcp(1 + hi::square_of(rDensity));
#else
  return 0;
#endif
}

bool Metal::GetScatterVector(tgir::PathVertex &eyePathVertex, tgir::Real const,
                             std::mt19937_64 &random) const {
  tgir::Vector3 vHalfVector;
  tgir::Real const nh = tgir::ComputeGlossyHalfwayVector(
      std::uniform_real_distribution<tgir::Real>()(random),
      std::uniform_real_distribution<tgir::Real>()(random), shininess_,
      eyePathVertex.vTangent, eyePathVertex.vShadingNormal,
      eyePathVertex.vBinormal, &vHalfVector);

  tgir::Real const nk1 =
      -hi::dot(eyePathVertex.vShadingNormal, eyePathVertex.vOutgoingDirection);
  tgir::Real const hk = -hi::dot(vHalfVector, eyePathVertex.vOutgoingDirection);
  eyePathVertex.vOutgoingDirection += (2 * hk) * vHalfVector;  // 反射光線
  tgir::Real const nk2 =
      hi::dot(eyePathVertex.vShadingNormal, eyePathVertex.vOutgoingDirection);
  if ((nk2 <= 0) || (hi::dot(eyePathVertex.vGeometricNormal,
                             eyePathVertex.vOutgoingDirection) <= 0)) {
    return false;  // 打ち止め
  }

  eyePathVertex.sQuantum *= tgir::FakeFresnelTerm(hk) / hi::max(nk1, nk2) * nk2;
#ifdef CONFIG_EXPLICIT_GLOSSY
  eyePathVertex.bSpecular = false;
  eyePathVertex.rSamplingPrev = 1;
  eyePathVertex.rSamplingNext =
      (shininess_ + 1) * (M_1_PI / 8) * std::pow(nh, shininess_) / hk;
#else
  eyePathVertex.bSpecular = true;
  eyePathVertex.rSamplingPrev = 1;
  eyePathVertex.rSamplingNext = 1;
#endif

  return true;
}

tgir::Real Metal::GetDensityVariance() const {
  return std::sqrt(1 +
                   (hi::square_of(TGIR_CONFIG_kMaxBranchCount) - 1) /
                       hi::square_of(shininess_ + 1));
}

bool Metal::GetScatterVector(tgir::PathVertex &pathVertex,
                             tgir::PrimarySample const &sample,
                             tgir::Vector3 const &mu,
                             bool const &bAdjoint) const {
  if (bAdjoint) {
    // ハーフベクトルを求める
    tgir::Vector3 vHalfVector;
    tgir::ComputeGlossyHalfwayVector(
        mu[0], mu[1], shininess_, pathVertex.vTangent,
        pathVertex.vGeometricNormal, pathVertex.vBinormal, &vHalfVector);

    tgir::Real const nk1 =
        -hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal);
    tgir::Real const hk = -hi::dot(pathVertex.vOutgoingDirection, vHalfVector);
    pathVertex.vOutgoingDirection += (2 * hk) * vHalfVector;
    tgir::Real const nk2 =
        hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal);
    if ((nk2 <= 0) || (hi::dot(pathVertex.vOutgoingDirection,
                               pathVertex.vShadingNormal) <= 0)) {
      return false;  // 打ち止め
    }

    pathVertex.sQuantum *= albedo_[sample.WavelengthIndex()] *
                           tgir::FakeFresnelTerm(hk) / hi::max(nk1, nk2) * nk2;
  } else {
    // ハーフベクトルを求める
    tgir::Vector3 vHalfVector;
    tgir::ComputeGlossyHalfwayVector(
        mu[0], mu[1], shininess_, pathVertex.vTangent,
        pathVertex.vShadingNormal, pathVertex.vBinormal, &vHalfVector);

    tgir::Real const nk1 =
        -hi::dot(pathVertex.vOutgoingDirection, pathVertex.vShadingNormal);
    tgir::Real const hk = -hi::dot(pathVertex.vOutgoingDirection, vHalfVector);
    pathVertex.vOutgoingDirection += (2 * hk) * vHalfVector;
    tgir::Real const nk2 =
        hi::dot(pathVertex.vOutgoingDirection, pathVertex.vShadingNormal);
    if ((nk2 <= 0) || (hi::dot(pathVertex.vOutgoingDirection,
                               pathVertex.vGeometricNormal) <= 0)) {
      return false;  // 打ち止め
    }

    pathVertex.sQuantum *= albedo_[sample.WavelengthIndex()] *
                           tgir::FakeFresnelTerm(hk) / hi::max(nk1, nk2) * nk2;
  }

  pathVertex.bSpecular = true;
  pathVertex.rSamplingPrev = 1;
  pathVertex.rSamplingNext = 1;

  return true;
}
}
// end of namespace tgir

namespace tgir {
/*
  virtual tgir::Real Metal::ProjectedPDF(
    tgir::Vector3 const & wi, ///< [in] incoming direction
    tgir::Vector3 const & wo, ///< [in] outgoing direction
    tgir::Vector3 const & ng, ///< [in] geometric normal
    tgir::Vector3 const & ns  ///< [in] shading normal
 ) const
  {
    // p(w) = p'(w) / |cos\theta|
    //   p'(w) = ph(h) / dot(w,h) / 4
    //   ph(h) = (s+1) / 2pi * dot(n,h)^{s}
    tgir::Vector3 const h = hi::normalize(wi + wo);
    tgir::Real const ph = std::pow(hi::dot(n,h), shininess_);
    return (M_1_PI/8) * (shininess_+1) * ph / hi::dot(n, wi);
  }
*/
tgir::Real Metal::BSDF(std::size_t const &is,    ///< [in] index of spectrum
                       tgir::Vector3 const &wi,  ///< [in] incoming direction
                       tgir::Vector3 const &wo,  ///< [in] outgoing direction
                       tgir::Vector3 const &ns   ///< [in] shading normal
                       ) const {
  // f_{s}(w_{i} -> w_{o})
  tgir::Vector3 const h = hi::normalize(wi + wo);
  tgir::Real const hw = hi::dot(h, wi);
  tgir::Real const hn = hi::dot(h, ns);
  return albedo_[is] * (M_1_PI / 8) * (shininess_ + 1) *
         std::pow(hn, shininess_) /
         (hw * std::max(hi::dot(ns, wi), hi::dot(ns, wo))) *
         tgir::FakeFresnelTerm(hw);
}
}
// end of namespace tgir
