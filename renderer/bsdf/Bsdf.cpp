#include "Bsdf.hpp"
#include "core/Scene.hpp"
#include "geom/PathVertex.hpp"
#include "geom/Intersection.hpp"
#include "sampler/PrimarySample.hpp"

namespace tgir {
bool Bsdf::BoundsImportance(tgir::PathVertex &, std::mt19937 &) const {
  return false;
};

bool Bsdf::BoundsPhoton(tgir::PathVertex &, std::mt19937_64 &) const {
  return false;
};

bool Bsdf::BoundsPhotonWithImportonMap(tgir::PathVertex &,
                                       tgir::ImportonMap const &,
                                       tgir::ImportonQuery &,
                                       std::mt19937_64 &) const {
  return false;
}

bool Bsdf::BoundsPhotonWithParticleFilter(tgir::PathVertex &,
                                          tgir::ImportonMap const &,
                                          tgir::ParticleQuery &,
                                          std::mt19937_64 &) const {
  return false;
}

bool Bsdf::GetScatterVector(tgir::PathVertex &, tgir::Real const,
                            tgir::ImportonMap const &, tgir::ImportanceQuery &,
                            std::mt19937_64 &) const {
  return false;
}

}  // end of namespace tgir
