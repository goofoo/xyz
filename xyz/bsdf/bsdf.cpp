#include "bsdf.hpp"
#include "../core/scene.hpp"
#include "../geom/pathvertex.hpp"
#include "../geom/intersection.hpp"

namespace xyz {
bool Bsdf::BoundsPhoton(PathVertex &, std::mt19937_64 &) const {
  return false;
};

bool Bsdf::BoundsPhotonWithImportonMap(PathVertex &, ImportonMap const &,
                                       ImportonQuery &,
                                       std::mt19937_64 &) const {
  return false;
}

bool Bsdf::BoundsPhotonWithParticleFilter(PathVertex &, ImportonMap const &,
                                          ParticleQuery &,
                                          std::mt19937_64 &) const {
  return false;
}

float3_t Bsdf::CalculateWeightedDirectLighting(
    __inout PathVertex *const pEyePathVertex,
    __in PathVertex const &lightPathVertex,
    __in ImportanceQuery const &importons, __in Scene const &scene) const {
  UNREFERENCED_PARAMETER(importons);
  return CalculateWeightedDirectLighting(pEyePathVertex, lightPathVertex,
                                         scene);
}

bool Bsdf::NextScatteringDirection(__inout PathVertex *const pEyePathVertex,
                                   __inout ImportanceQuery const &importons,
                                   __inout IPrimarySample &sample) const {
  UNREFERENCED_PARAMETER(importons);
  return NextScatteringDirection(pEyePathVertex, sample);
}

bool Bsdf::CameraToLightsScatteringDirection(
    __inout PathVertex *const pPathVertex,
    __in ImportanceQuery const &importons,
    __inout IPrimarySample &sample) const {
  UNREFERENCED_PARAMETER(importons);
  return CameraToLightsScatteringDirection(pPathVertex, sample);
}

}  // end of namespace xyz
