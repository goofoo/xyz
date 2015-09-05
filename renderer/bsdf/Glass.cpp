#include "Glass.hpp"
#include "core/Scene.hpp"
#include "geom/PathVertex.hpp"
#include "geom/Intersection.hpp"
#include "sampler/PrimarySample.hpp"

namespace tgir {
void Glass::Render() const {
  hi::basic_vector4<GLfloat> const color(0.3465f, 0.222f, 0.4402f,
                                         0.5f);  // Patch 10, Purple
  ::glEnable(GL_LIGHTING);
  ::glEnable(GL_BLEND);
  ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  hi::glMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color.data());
}

tgir::Spectrum Glass::CalculateWeightedDirectLighting(
    tgir::PathVertex &, tgir::PathVertex const &, std::size_t const,
    tgir::Scene const &) const {
  return 0;
}

bool Glass::GetScatterVector(tgir::PathVertex &pathVertex,
                             tgir::Real const wavelength,
                             std::mt19937_64 &random) const {
#if 0
    static tgir::Real const B1 = 1.738484030;
    static tgir::Real const B2 = 3.11168974e-1;
    static tgir::Real const B3 = 1.174908710;
    static tgir::Real const C1 = 1.36068604e-2;
    static tgir::Real const C2 = 6.15960463e-2;
    static tgir::Real const C3 = 1.21922711e+2;
    tgir::Real const L_sq = hi::square_of(0.4 + (0.7 - 0.4) * wavelength); // μm
    tgir::Real const mu_sq = (B1*L_sq/(L_sq-C1)) + (B2*L_sq/(L_sq-C2)) + (B3*L_sq/(L_sq-C3)) + 1;
    tgir::Real const ior_sq = pathVertex.back ? hi::rcp(mu_sq) : mu_sq;
#else
  tgir::Real const ior =
      1.400 + 50.0 / ((400 - 230) + (700 - 400) * wavelength);
  tgir::Real const ior_sq =
      hi::square_of(pathVertex.back_side ? hi::rcp(ior) : ior);
#endif
  tgir::Real const nk =
      -hi::dot(pathVertex.vOutgoingDirection, pathVertex.vShadingNormal);
  tgir::Real const b = ior_sq + nk * nk - 1;

  if (b > 0) {
    tgir::Real const g = std::sqrt(b);
    tgir::Real const gpc = g + nk;
    tgir::Real const gmc = g - nk;
    tgir::Real const alpha = gmc / gpc;
    tgir::Real const beta = (nk * gpc - 1) / (nk * gmc + 1);
    tgir::Real const fr = tgir::Real(0.5) * alpha * alpha * (beta * beta + 1);

    if (std::uniform_real_distribution<tgir::Real>()(random) < fr) {
      tgir::ComputeReflectedVector(pathVertex.vShadingNormal, nk,
                                   &pathVertex.vOutgoingDirection);  // 反射
      if (hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal) <=
          0) {
        return false;  // 打ち止め
      }
    } else {
      tgir::ComputeRefractedVector(pathVertex.vShadingNormal, gmc,
                                   std::sqrt(ior_sq),
                                   pathVertex.vOutgoingDirection);  // 屈折
      if (hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal) >=
          0) {
        return false;  // 打ち止め
      }
    }
  } else {
    tgir::ComputeReflectedVector(pathVertex.vShadingNormal, nk,
                                 &pathVertex.vOutgoingDirection);  // 完全反射
    if (hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal) <=
        0) {
      return false;  // 打ち止め
    }
  }

  pathVertex.bSpecular = true;
  pathVertex.rSamplingPrev = 1;
  pathVertex.rSamplingNext = 1;

  return true;
}

tgir::Real Glass::GetDensityVariance() const { return 2; }

bool Glass::GetScatterVector(tgir::PathVertex &pathVertex,
                             tgir::PrimarySample const &sample,
                             tgir::Vector3 const &mu, bool const &) const {
  tgir::Real const ior =
      1.400 + 50.0 / ((400 - 230) + (700 - 400) * sample.WavelengthMu());
  tgir::Real const ior_sq =
      hi::square_of(pathVertex.back_side ? hi::rcp(ior) : ior);
  tgir::Real const nk =
      -hi::dot(pathVertex.vOutgoingDirection, pathVertex.vShadingNormal);
  tgir::Real const b = ior_sq + nk * nk - tgir::Real(1);

  if (b > 0) {
    tgir::Real const g = std::sqrt(b);
    tgir::Real const gpc = g + nk;
    tgir::Real const gmc = g - nk;
    tgir::Real const alpha = gmc / gpc;
    tgir::Real const beta = (nk * gpc - 1) / (nk * gmc + 1);
    tgir::Real const fr = tgir::Real(0.5) * alpha * alpha * (beta * beta + 1);

    if (mu[2] < fr) {
      tgir::ComputeReflectedVector(pathVertex.vShadingNormal, nk,
                                   &pathVertex.vOutgoingDirection);  // 反射
      if (hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal) <=
          0) {
        return false;  // 打ち止め
      }

      pathVertex.rSamplingPrev = fr;
      pathVertex.rSamplingNext = fr;
    } else {
      tgir::ComputeRefractedVector(pathVertex.vShadingNormal, gmc,
                                   std::sqrt(ior_sq),
                                   pathVertex.vOutgoingDirection);  // 屈折
      if (hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal) >=
          0) {
        return false;  // 打ち止め
      }

      pathVertex.rSamplingPrev = 1 - fr;
      pathVertex.rSamplingNext = 1 - fr;
    }
  } else {
    tgir::ComputeReflectedVector(pathVertex.vShadingNormal, nk,
                                 &pathVertex.vOutgoingDirection);  // 完全反射
    if (hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal) <=
        0) {
      return false;  // 打ち止め
    }

    pathVertex.rSamplingPrev = 1;
    pathVertex.rSamplingNext = 1;
  }

  pathVertex.bSpecular = true;

  return true;
}

}  // end of namespace tgir
