#ifndef XYZ_BSDF_GLASS_HPP_
#define XYZ_BSDF_GLASS_HPP_

#include "bsdf.hpp"

namespace xyz {
//
// 完全反射透過
//
class Glass : public Bsdf {
  HI_DISALLOW_COPY_AND_ASSIGN(Glass);

 public:
  inline Glass() {}

  virtual void Render() const;

  virtual bool BoundsImportance(__inout PathVertex *const pLightPathVertex,
                                __inout std::mt19937_64 &random) const;

  //
  // 経路追跡法用
  //
  virtual Bsdf::type What() const { return Bsdf::GLASS; }
  virtual float3_t CalculateWeightedDirectLighting(__inout PathVertex *const,
                                                   __in PathVertex const &,
                                                   __in Scene const &) const;
  virtual bool NextScatteringDirection(__inout PathVertex *const,
                                       __inout IPrimarySample &sample) const;
  virtual float_t GetDensityVariance() const;

  //
  // 双方向経路追跡法用
  //
  virtual bool LightsToCameraScatteringDirection(
      __inout PathVertex *const pPathVertex,
      __inout IPrimarySample &sample) const;

  virtual bool CameraToLightsScatteringDirection(
      __inout PathVertex *const pPathVertex,
      __inout IPrimarySample &sample) const;

  //
  // 双方向経路追跡用
  //
  virtual bool NextScatteringDirection(PathVertex &,      ///< [in/out]
                                       float3_t const &,  ///< [in]
                                       bool const &       ///< [in]
                                       ) const;

 public:
  virtual bool HasNoDiracDistribution() const { return false; }
  virtual bool IsLightEmittor() const { return false; }

  virtual float_t BSDF(std::size_t const &is,  ///< [in] index of spectrum
                       float3_t const &wi,     ///< [in] incoming direction
                       float3_t const &wo,     ///< [in] outgoing direction
                       float3_t const &ns      ///< [in] shading normal
                       ) const {
    UNREFERENCED_PARAMETER(is);
    UNREFERENCED_PARAMETER(wi);
    UNREFERENCED_PARAMETER(wo);
    UNREFERENCED_PARAMETER(ns);
    return 0;
  }
};

}  // end of namespace xyz

#endif XYZ_BSDF_GLASS_HPP_
