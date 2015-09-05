#ifndef __TGIR_BSDF_MIRROR_HPP__
#define __TGIR_BSDF_MIRROR_HPP__

#include "Bsdf.hpp"

namespace tgir {
//
// 完全鏡面反射
//
class Mirror : public Bsdf {
  HI_DISALLOW_COPY_AND_ASSIGN(Mirror);

 public:
  inline Mirror() {}

  virtual void Render() const;

  //
  // 経路追跡法用
  //
  virtual Bsdf::type What() const { return Bsdf::MIRROR; }
  virtual tgir::Spectrum CalculateWeightedDirectLighting(
      tgir::PathVertex &, tgir::PathVertex const &, std::size_t const,
      tgir::Scene const &) const;
  virtual bool GetScatterVector(tgir::PathVertex &, tgir::Real const,
                                std::mt19937_64 &) const;
  virtual tgir::Real GetDensityVariance() const;

  //
  // 双方向経路追跡用
  //
  virtual bool GetScatterVector(tgir::PathVertex &,           ///< [in/out]
                                tgir::PrimarySample const &,  ///< [in]
                                tgir::Vector3 const &,        ///< [in]
                                bool const &                  ///< [in]
                                ) const;

 public:
  virtual bool HasNoDiracDistribution() const { return false; }
  virtual bool IsLightEmittor() const { return false; }

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

#endif __TGIR_BSDF_MIRROR_HPP__
