#ifndef XYZ_BSDF_LIGHT_HPP_
#define XYZ_BSDF_LIGHT_HPP_

#include "bsdf.hpp"

namespace xyz {
//
// åıåπ
//
class Light : public Bsdf {
  HI_DISALLOW_COPY_AND_ASSIGN(Light);

 public:
  inline Light() {}

  virtual void Render() const;

  virtual bool BoundsImportance(PathVertex *const pLightPathVertex,
                                std::mt19937_64 &random) const;

  //
  // åoòHí«ê’ñ@óp
  //
  virtual type What() const { return Bsdf::LIGHT; }
  virtual float3_t CalculateWeightedDirectLighting(__inout PathVertex *const,
                                                   __in PathVertex const &,
                                                   __in Scene const &) const;
  virtual bool NextScatteringDirection(__inout PathVertex *const,
                                       __inout IPrimarySample &sample) const;
  virtual float_t GetDensityVariance() const;

  //
  // ëoï˚å¸åoòHí«ê’ñ@óp
  //
  virtual bool LightsToCameraScatteringDirection(
      __inout PathVertex *const pPathVertex,
      __inout IPrimarySample &sample) const;

  virtual bool CameraToLightsScatteringDirection(
      __inout PathVertex *const pPathVertex,
      __inout IPrimarySample &sample) const;

  //
  // ëoï˚å¸åoòHí«ê’óp
  //
  bool NextScatteringDirection(PathVertex &,      ///< [in/out]
                               float3_t const &,  ///< [in]
                               bool const &       ///< [in]
                               ) const;

 public:
  virtual bool HasNoDiracDistribution() const { return false; }
  virtual bool IsLightEmittor() const { return true; }

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

#endif XYZ_BSDF_LIGHT_HPP_
