#include "Light.hpp"

namespace tgir {
void Light::Render() const {
  hi::basic_vector4<GLfloat> const color(1);
  ::glDisable(GL_LIGHTING);
  ::glDisable(GL_BLEND);
  hi::glMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color.data());
}

tgir::Spectrum Light::CalculateWeightedDirectLighting(
    tgir::PathVertex &, tgir::PathVertex const &, std::size_t const,
    tgir::Scene const &) const {
  return 0;
}

bool Light::GetScatterVector(tgir::PathVertex &, tgir::Real const,
                             std::mt19937_64 &) const {
  return false;
}

tgir::Real Light::GetDensityVariance() const { return 0; }

bool Light::GetScatterVector(tgir::PathVertex &,           ///< [in/out]
                             tgir::PrimarySample const &,  ///< [in]
                             tgir::Vector3 const &,        ///< [in]
                             bool const &                  ///< [in]
                             ) const {
  return false;
}
}  // end of namespace tgir
