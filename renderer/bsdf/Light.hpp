#ifndef __TGIR_BSDF_LIGHT_HPP__
#define __TGIR_BSDF_LIGHT_HPP__

#include "Bsdf.hpp"

namespace tgir {
//
// åıåπ
//
class Light : public Bsdf {
  HI_DISALLOW_COPY_AND_ASSIGN(Light);

 public:
  inline Light() {}

  virtual void Render() const;

  //
  // åoòHí«ê’ñ@óp
  //
  virtual tgir::Spectrum CalculateWeightedDirectLighting(
      tgir::PathVertex &, tgir::PathVertex const &, std::size_t const,
      tgir::Scene const &) const;

  virtual bool GetScatterVector(tgir::PathVertex &, tgir::Real const,
                                std::mt19937_64 &) const;

  virtual type What() const { return Bsdf::LIGHT; }

  virtual tgir::Real GetDensityVariance() const;

  //
  // ëoï˚å¸åoòHí«ê’óp
  //
  bool GetScatterVector(tgir::PathVertex &,           ///< [in/out]
                        tgir::PrimarySample const &,  ///< [in]
                        tgir::Vector3 const &,        ///< [in]
                        bool const &                  ///< [in]
                        ) const;

 public:
  virtual bool HasNoDiracDistribution() const { return false; }
  virtual bool IsLightEmittor() const { return true; }

  virtual tgir::Real BSDF(std::size_t const &is,    ///< [in] index of spectrum
                          tgir::Vector3 const &wi,  ///< [in] incoming direction
                          tgir::Vector3 const &wo,  ///< [in] outgoing direction
                          tgir::Vector3 const &ns   ///< [in] shading normal
                          ) const {
    UNREFERENCED_PARAMETER(is);
    UNREFERENCED_PARAMETER(wi);
    UNREFERENCED_PARAMETER(wo);
    UNREFERENCED_PARAMETER(ns);
    return 0;
  }
};

}  // end of namespace tgir

#endif __TGIR_BSDF_LIGHT_HPP__
