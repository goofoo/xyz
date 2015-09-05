#ifndef __TGIR_SAMPLER_PRIMARYSAMPLE_HPP__
#define __TGIR_SAMPLER_PRIMARYSAMPLE_HPP__

#include "core/config.hpp"

namespace tgir {
/// <summary>
/// the sample in primary sample space
/// </summary>
class PrimarySample {
 public:
  // n: the number of eye-path vertices
  inline explicit PrimarySample(std::size_t const path_length)
      : muLight_(path_length + 1)
#ifdef CONFIG_PINHOLE_CAMERA_MODEL
        ,
        muEye_(path_length)
#else
        ,
        muEye_(path_length + 1)
#endif CONFIG_PINHOLE_CAMERA_MODEL
  {
  }

  inline explicit PrimarySample(PrimarySample const &rhs)
      : muWavelength_(rhs.muWavelength_),
        muLight_(rhs.muLight_),
        muEye_(rhs.muEye_) {}

  inline PrimarySample &operator=(PrimarySample const &rhs) {
    if (this != &rhs) {
      muWavelength_ = rhs.muWavelength_;
      muLight_ = rhs.muLight_;
      muEye_ = rhs.muEye_;
    }

    return *this;
  }

  inline void swap(PrimarySample &rhs) {
    std::swap(muWavelength_, rhs.muWavelength_);
    std::swap(muLight_, rhs.muLight_);
    std::swap(muEye_, rhs.muEye_);
  }

 public:
  void Init(std::mt19937_64 &);

  // stratification
  inline void Stratify(std::size_t const x, std::size_t const w,
                       std::size_t const y, std::size_t const h) {
    muEye_[0][0] = (x + muEye_[0][0]) / w;
    muEye_[0][1] = (y + muEye_[0][1]) / h;
  }

  void Mutate(tgir::PrimarySample const &, tgir::Real const,
              std::mt19937_64 &random);
  void MutateFocusOnLightSubPath(tgir::PrimarySample const &, tgir::Real const,
                                 std::mt19937_64 &random);
  void MutateFocusOnEyeSubPath(tgir::PrimarySample const &, tgir::Real const,
                               std::mt19937_64 &random);
  void MutateEyeSubPath(tgir::PrimarySample const &, tgir::Real const,
                        std::mt19937_64 &random);
  void MutateLightSubPath(tgir::PrimarySample const &, tgir::Real const,
                          std::mt19937_64 &random);

 public:
  inline std::size_t GetPixelPosition(std::size_t const w,
                                      std::size_t const h) const {
    std::size_t const x = static_cast<std::size_t>(muEye_[0][0] * w);
    std::size_t const y = static_cast<std::size_t>(muEye_[0][1] * h);
    return y * w + x;
  }

  inline tgir::Real WavelengthMu() const { return muWavelength_.first; }
  inline std::size_t const &WavelengthIndex() const {
    return muWavelength_.second;
  }
  inline tgir::Vector3 const &L(std::size_t const &s) const {
    return muLight_[s];
  }
  inline tgir::Vector3 const &E(std::size_t const &t) const {
    return muEye_[t];
  }

 private:
  std::pair<tgir::Real, std::size_t> muWavelength_;
  std::vector<tgir::Vector3> muLight_;
  std::vector<tgir::Vector3> muEye_;
};

}  // end of namespace tgir

namespace std {
template <>
inline void swap(tgir::PrimarySample &lhs, tgir::PrimarySample &rhs) {
  lhs.swap(rhs);
}

}  // end of namespace std

#endif __TGIR_SAMPLER_PRIMARYSAMPLE_HPP__
