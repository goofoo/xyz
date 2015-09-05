#ifndef __TGIR_BSDF_GLASS_HPP__
#define __TGIR_BSDF_GLASS_HPP__

#include "Bsdf.hpp"

namespace tgir {
//
// ���S���˓���
//
class Glass : public Bsdf {
  HI_DISALLOW_COPY_AND_ASSIGN(Glass);

 public:
  inline Glass() {}

  virtual void Render() const;

  //
  // �o�H�ǐՖ@�p
  //
  virtual Bsdf::type What() const { return Bsdf::GLASS; }
  virtual tgir::Spectrum CalculateWeightedDirectLighting(
      tgir::PathVertex &, tgir::PathVertex const &, std::size_t const,
      tgir::Scene const &) const;
  virtual bool GetScatterVector(tgir::PathVertex &, tgir::Real const,
                                std::mt19937_64 &) const;
  virtual tgir::Real GetDensityVariance() const;

  //
  // �o�����o�H�ǐ՗p
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

#endif __TGIR_BSDF_GLASS_HPP__
