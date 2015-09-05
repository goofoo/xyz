#include "light.hpp"

namespace xyz {
void Light::Render() const {
  hi::basic_vector4<GLfloat> const color(1);
  ::glDisable(GL_LIGHTING);
  ::glDisable(GL_BLEND);
  hi::glMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color.data());
}

bool Light::BoundsImportance(PathVertex *const pLightPathVertex,
                             std::mt19937_64 &random) const {
  UNREFERENCED_PARAMETER(pLightPathVertex);
  UNREFERENCED_PARAMETER(random);
  return false;
}

float3_t Light::CalculateWeightedDirectLighting(
    __inout PathVertex *const pEyePathVertex,
    __in PathVertex const &lightPathVertex, __in Scene const &) const {
  UNREFERENCED_PARAMETER(pEyePathVertex);
  UNREFERENCED_PARAMETER(lightPathVertex);
  // UNREFERENCED_PARAMETER(scene);
  return float3_t(0);
}

bool Light::NextScatteringDirection(__inout PathVertex *const,
                                    __inout IPrimarySample &sample) const {
  UNREFERENCED_PARAMETER(sample);
  return false;
}

float_t Light::GetDensityVariance() const { return 0; }

bool Light::NextScatteringDirection(PathVertex &,      ///< [in/out]
                                    float3_t const &,  ///< [in]
                                    bool const &       ///< [in]
                                    ) const {
  return false;
}
}  // end of namespace xyz

namespace xyz {
bool Light::LightsToCameraScatteringDirection(
    __inout PathVertex *const pPathVertex,
    __inout IPrimarySample &sample) const {
  UNREFERENCED_PARAMETER(pPathVertex);
  UNREFERENCED_PARAMETER(sample);

  return false;
}

bool Light::CameraToLightsScatteringDirection(
    __inout PathVertex *const pPathVertex,
    __inout IPrimarySample &sample) const {
  UNREFERENCED_PARAMETER(pPathVertex);
  UNREFERENCED_PARAMETER(sample);

  return false;
}
}  // end of namespace xyz
