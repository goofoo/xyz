#include "Mirror.hpp"
#include "core/Scene.hpp"
#include "geom/PathVertex.hpp"
#include "geom/Intersection.hpp"
#include "sampler/PrimarySample.hpp"

namespace tgir {
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

tgir::Spectrum Mirror::CalculateWeightedDirectLighting(
    tgir::PathVertex &, tgir::PathVertex const &, std::size_t const,
    tgir::Scene const &) const {
  return 0;
}

bool Mirror::GetScatterVector(tgir::PathVertex &pathVertex, tgir::Real const,
                              std::mt19937_64 &) const {
  // 入力時のpathVertex.vOutgoingDirectionには，IncomingDirectionが格納されている．
  tgir::Real const cosTheta =
      -hi::dot(pathVertex.vOutgoingDirection,
               pathVertex.vShadingNormal);  // \cos\theta

  // 反射ベクトルを求める
  tgir::ComputeReflectedVector(pathVertex.vShadingNormal, cosTheta,
                               &pathVertex.vOutgoingDirection);
  if (hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal) <=
      0) {
    return false;  // 打ち止め
  }

#if 0
    // 薄膜の厚さ
    static tgir::Real const thinness = 1090;

    // 波長
    tgir::Real const lambda = 400 + (700 - 400) * wavelength;

    // 薄膜の屈折率
    tgir::Real const ior = 1.400 + 50.0 / (lambda - 230);
    tgir::Real const iorSquared = hi::square_of(ior);

    tgir::Real const sinThetaSquared = tgir::Real(1) - hi::square_of(cosTheta); // \sin^2\theta

    tgir::Real const delta = (2*M_PI) * thinness / lambda * std::sqrt(iorSquared - sinThetaSquared); // \mu_{1} = 1 としている(真空なので)

    // 金属反射なので \mu_{1} < \mu_{f} < \mu_{2} よって deletaDash = 0
    // \mu_{1} = 1;
    // \mu_{f} = 1.4 以上
    // \mu_{2} = 金属なので100ぐらい?
    static tgir::Real const deletaDash = 0;

    tgir::Real const cosDeltaSquared = hi::square_of(std::cos(delta+deletaDash));

    pathVertex.sQuantum *= cosDeltaSquared;
#endif

  pathVertex.bSpecular = true;
  pathVertex.rSamplingPrev = 1;
  pathVertex.rSamplingNext = 1;

  return true;
}

tgir::Real Mirror::GetDensityVariance() const { return 1; }

bool Mirror::GetScatterVector(tgir::PathVertex &pathVertex,
                              tgir::PrimarySample const &,
                              tgir::Vector3 const &, bool const &) const {
  // 入力時のpathVertex.vOutgoingDirectionには，IncomingDirectionが格納されている．
  tgir::Real const cosTheta =
      -hi::dot(pathVertex.vOutgoingDirection,
               pathVertex.vShadingNormal);  // \cos\theta

  // 反射ベクトルを求める
  tgir::ComputeReflectedVector(pathVertex.vShadingNormal, cosTheta,
                               &pathVertex.vOutgoingDirection);

  if (hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal) <=
      0) {
    return false;  // 打ち止め
  }

#if 0
    // 薄膜の厚さ
    static tgir::Real const thinness = 1090;

    // 波長
    tgir::Real const lambda = 400 + (700 - 400) * sample.WavelengthMu();

    // 薄膜の屈折率
    tgir::Real const ior = 1.400 + 50.0 / (lambda - 230);
    tgir::Real const iorSquared = hi::square_of(ior);

    tgir::Real const sinThetaSquared = tgir::Real(1) - hi::square_of(cosTheta); // \sin^2\theta

    tgir::Real const delta = (2*M_PI) * thinness / lambda * std::sqrt(iorSquared - sinThetaSquared); // \mu_{1} = 1 としている(真空なので)

    // 金属反射なので \mu_{1} < \mu_{f} < \mu_{2} よって deletaDash = 0
    // \mu_{1} = 1;
    // \mu_{f} = 1.4 以上
    // \mu_{2} = 金属なので100ぐらい?
    static tgir::Real const deletaDash = 0;

    tgir::Real const cosDeltaSquared = hi::square_of(std::cos(delta+deletaDash));

    pathVertex.sQuantum *= cosDeltaSquared;
#endif

  pathVertex.bSpecular = true;
  pathVertex.rSamplingPrev = 1;
  pathVertex.rSamplingNext = 1;

  return true;
}

}  // end of namespace tgir
