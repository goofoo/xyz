#ifndef XYZ_BSDF_METAL_HPP_
#define XYZ_BSDF_METAL_HPP_

#include "bsdf.hpp"

namespace xyz {
//
// Neumann-Neumann BRDF
//
class Metal : public Bsdf {
  HI_DISALLOW_COPY_AND_ASSIGN(Metal);

 public:
  inline Metal(float_t const shininess, float3_t const &albedo)
      : shininess_(shininess), albedo_(albedo) {}

  virtual void Render() const;

  inline float_t sininess() const { return shininess_; }
  inline float3_t albedo() const { return albedo_; }

  virtual bool BoundsImportance(PathVertex *const pLightPathVertex,
                                std::mt19937_64 &random) const;

  //
  // åoòHí«ê’ñ@óp
  //
  virtual Bsdf::type What() const { return Bsdf::METAL; }
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
  virtual bool NextScatteringDirection(PathVertex &,      ///< [in/out]
                                       float3_t const &,  ///< [in]
                                       bool const &       ///< [in]
                                       ) const;

 public:
  virtual bool HasNoDiracDistribution() const { return true; }
  virtual bool IsLightEmittor() const { return false; }

  virtual float_t BSDF(std::size_t const &is,  ///< [in] index of spectrum
                       float3_t const &wi,     ///< [in] incoming direction
                       float3_t const &wo,     ///< [in] outgoing direction
                       float3_t const &ns      ///< [in] shading normal
                       ) const;

 private:
  float_t const shininess_;
  float3_t const albedo_;
};

}  // end of namespace xyz

#endif XYZ_BSDF_METAL_HPP_
