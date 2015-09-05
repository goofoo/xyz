#ifndef XYZ_BSDF_LAMBERT_HPP_
#define XYZ_BSDF_LAMBERT_HPP_

#include "bsdf.hpp"

namespace xyz {
//
// 完全拡散反射
//
class Lambert : public Bsdf {
 public:
  inline Lambert(float3_t const &albedo) : albedo_(albedo) {}

  virtual void Render() const;

  virtual float3_t albedo() const { return albedo_; }

  //
  // フォトンマップ構築用
  //
  virtual bool BoundsPhoton(PathVertex &, std::mt19937_64 &) const;
  virtual bool BoundsPhotonWithImportonMap(PathVertex &, ImportonMap const &,
                                           ImportonQuery &,
                                           std::mt19937_64 &) const;
  virtual bool BoundsPhotonWithParticleFilter(PathVertex &, ImportonMap const &,
                                              ParticleQuery &,
                                              std::mt19937_64 &) const;

  //
  // 経路追跡法用
  //
  virtual Bsdf::type What() const { return Bsdf::LAMBERT; }

  virtual float3_t CalculateWeightedDirectLighting(__inout PathVertex *const,
                                                   __in PathVertex const &,
                                                   __in Scene const &) const;

  virtual bool NextScatteringDirection(__inout PathVertex *const,
                                       __inout IPrimarySample &sample) const;

  virtual float3_t CalculateWeightedDirectLighting(
      __inout PathVertex *const pEyePathVertex,
      __in PathVertex const &lightPathVertex,
      __in ImportanceQuery const &importons, __in Scene const &scene) const;

  virtual bool NextScatteringDirection(__inout PathVertex *const pEyePathVertex,
                                       __in ImportanceQuery const &importons,
                                       __inout IPrimarySample &sample) const;

  virtual bool BoundsImportance(__inout PathVertex *const pLightPathVertex,
                                __inout std::mt19937_64 &random) const;

  virtual float_t GetDensityVariance() const;

  //
  // 双方向経路追跡用
  //
  virtual bool LightsToCameraScatteringDirection(
      __inout PathVertex *const pPathVertex,
      __inout IPrimarySample &sample) const;

  virtual bool CameraToLightsScatteringDirection(
      __inout PathVertex *const pPathVertex,
      __inout IPrimarySample &sample) const;

  virtual bool Lambert::CameraToLightsScatteringDirection(
      __inout PathVertex *const pPathVertex,
      __in ImportanceQuery const &importons,
      __inout IPrimarySample &sample) const;

  //
  // 双方向経路追跡用
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
  float3_t const albedo_;

  HI_DISALLOW_COPY_AND_ASSIGN(Lambert);
};

}  // end of namespace xyz

#endif XYZ_BSDF_LAMBERT_HPP_
