#ifndef __TGIR_BSDF_LAMBERT_HPP__
#define __TGIR_BSDF_LAMBERT_HPP__

#include "Bsdf.hpp"

namespace tgir {
//
// 完全拡散反射
//
class Lambert : public Bsdf {
  HI_DISALLOW_COPY_AND_ASSIGN(Lambert);

 public:
  inline Lambert(tgir::SpectrumData const &albedo) : albedo_(albedo) {}

  virtual void Render() const;

  virtual tgir::SpectrumData::value_type Albedo(
      std::size_t const uWavelengthIndex) const {
    return albedo_[uWavelengthIndex];
  }

  //
  // フォトンマップ構築用
  //
  virtual bool BoundsImportance(tgir::PathVertex &, std::mt19937_64 &) const;
  virtual bool BoundsPhoton(tgir::PathVertex &, std::mt19937_64 &) const;
  virtual bool BoundsPhotonWithImportonMap(tgir::PathVertex &,
                                           tgir::ImportonMap const &,
                                           tgir::ImportonQuery &,
                                           std::mt19937_64 &) const;
  virtual bool BoundsPhotonWithParticleFilter(tgir::PathVertex &,
                                              tgir::ImportonMap const &,
                                              tgir::ParticleQuery &,
                                              std::mt19937_64 &) const;

  //
  // 経路追跡法用
  //
  virtual Bsdf::type What() const { return Bsdf::LAMBERT; }
  virtual tgir::Spectrum CalculateWeightedDirectLighting(
      tgir::PathVertex &, tgir::PathVertex const &, std::size_t const,
      tgir::Scene const &) const;
  virtual bool GetScatterVector(tgir::PathVertex &, tgir::Real const,
                                std::mt19937_64 &) const;
  virtual bool GetScatterVector(tgir::PathVertex &, tgir::Real const,
                                tgir::ImportonMap const &,
                                tgir::ImportanceQuery &,
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
  virtual bool HasNoDiracDistribution() const { return true; }

  virtual bool IsLightEmittor() const { return false; }

  virtual tgir::Real BSDF(std::size_t const &is,    ///< [in] index of spectrum
                          tgir::Vector3 const &wi,  ///< [in] incoming direction
                          tgir::Vector3 const &wo,  ///< [in] outgoing direction
                          tgir::Vector3 const &ns   ///< [in] shading normal
                          ) const;

 private:
  tgir::SpectrumData const albedo_;
};

}  // end of namespace tgir

#endif __TGIR_BSDF_LAMBERT_HPP__
