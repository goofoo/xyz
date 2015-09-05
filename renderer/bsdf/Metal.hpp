#ifndef __TGIR_BSDF_METAL_HPP__
#define __TGIR_BSDF_METAL_HPP__

#include "Bsdf.hpp"

namespace tgir {
//
// Neumann-Neumann BRDF
//
class Metal : public Bsdf {
  HI_DISALLOW_COPY_AND_ASSIGN(Metal);

 public:
  inline Metal(tgir::Real const shininess, tgir::SpectrumData const &albedo)
      : shininess_(shininess), albedo_(albedo) {}

  virtual void Render() const;

  inline tgir::SpectrumData::value_type Albedo(std::size_t const ior) const {
    return albedo_[ior];
  }
  inline tgir::Real const Sininess() const { return shininess_; }

  //
  // åoòHí«ê’ñ@óp
  //
  virtual Bsdf::type What() const { return Bsdf::METAL; }
  virtual tgir::Spectrum CalculateWeightedDirectLighting(
      tgir::PathVertex &, tgir::PathVertex const &, std::size_t const,
      tgir::Scene const &) const;
  virtual bool GetScatterVector(tgir::PathVertex &, tgir::Real const,
                                std::mt19937_64 &) const;
  virtual tgir::Real GetDensityVariance() const;

  //
  // ëoï˚å¸åoòHí«ê’óp
  //
  virtual bool GetScatterVector(tgir::PathVertex &,           ///< [in/out]
                                tgir::PrimarySample const &,  ///< [in]
                                tgir::Vector3 const &,        ///< [in]
                                bool const &                  ///< [in]
                                ) const;

 public:
  virtual bool HasNoDiracDistribution() const { return true; }
  virtual bool IsLightEmittor() const { return false; }

  virtual tgir::Real BSDF(std::size_t const &is,    ///< [in] index of spectrum
                          tgir::Vector3 const &wi,  ///< [in] incoming direction
                          tgir::Vector3 const &wo,  ///< [in] outgoing direction
                          tgir::Vector3 const &ns   ///< [in] shading normal
                          ) const;

 private:
  tgir::Real const shininess_;
  tgir::SpectrumData const albedo_;
};

}  // end of namespace tgir

#endif __TGIR_BSDF_METAL_HPP__
