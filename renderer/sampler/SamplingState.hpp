#ifndef __TGIR_SAMPLER_SAMPLINGSTATE_HPP__
#define __TGIR_SAMPLER_SAMPLINGSTATE_HPP__

#include "PrimarySample.hpp"
#include "core/Film.hpp"
#include "core/config.hpp"
#include "geom/Path.hpp"

namespace tgir {
template <typename T>
class SamplingState {
  typedef std::vector<tgir::PixelDescriptor> SamplesVector;
  typedef T ReplicaVector;

 public:
  inline SamplingState()
      : sample(TGIR_CONFIG_kMaxPathLength),
        paths(TGIR_CONFIG_kMaxPathLength),
        values(TGIR_CONFIG_kMaxPathLength) {}

  inline SamplingState &operator=(SamplingState const &rhs) {
    if (this != &rhs) {
      sample = rhs.sample;
      paths = rhs.paths;
      values = rhs.values;
      special = rhs.special;
      q = rhs.q;
    }

    return *this;
  }

  inline void swap(SamplingState &rhs) {
    std::swap(sample, rhs.sample);
    std::swap(paths, rhs.paths);
    std::swap(values, rhs.values);
    std::swap(special, rhs.special);
    std::swap(q, rhs.q);
  }

 public:
  tgir::PrimarySample sample;     ///< in primary sample space
  tgir::Path paths;               ///< in path space
  SamplesVector values;           ///< evaluted sample
  tgir::PixelDescriptor special;  ///< evaluted special sample
  ReplicaVector q;                ///< luminance
};

}  // end of namespace tgir

namespace std {
template <typename T>
inline void swap(tgir::SamplingState<T> &lhs, tgir::SamplingState<T> &rhs) {
  lhs.swap(rhs);
}
}

#endif __TGIR_SAMPLER_SAMPLINGSTATE_HPP__
