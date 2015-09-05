#ifndef XYZ_BSDF_BSDF_HPP_
#define XYZ_BSDF_BSDF_HPP_

#include "../core/PhotonMap.hpp"

namespace xyz {
class PathVertex;
class Scene;

class IPrimarySample {
 public:
  virtual ~IPrimarySample() {}
  virtual float_t next() = 0;
};

//
// base of non-symmetric BSDF
//
class Bsdf {
 public:
  enum type {
    LAMBERT,
    MIRROR,
    GLASS,
    METAL,
    LIGHT,
    COUNT,
  };

  virtual void Render() const = 0;

  //
  // フォトンマップ構築用
  //
  /// 光源または視点からの重要度を追跡する
  /// フォトンを追跡する
  virtual bool BoundsPhoton(PathVertex &lightPathVertex,
                            std::mt19937_64 &random) const;
  /// 視点からの重要度を使ってフォトンを追跡する
  virtual bool BoundsPhotonWithImportonMap(PathVertex &lightPathVertex,
                                           ImportonMap const &importonMap,
                                           ImportonQuery &importons,
                                           std::mt19937_64 &random) const;
  /// 粒子フィルタにより視点からの重要度を使ってフォトンを追跡する
  virtual bool BoundsPhotonWithParticleFilter(PathVertex &lightPathVertex,
                                              ImportonMap const &importonMap,
                                              ParticleQuery &importons,
                                              std::mt19937_64 &random) const;

  //
  // 経路追跡法用
  //
  virtual Bsdf::type What() const = 0;

  virtual float3_t CalculateWeightedDirectLighting(
      __inout PathVertex *const pEyePathVertex,
      __in PathVertex const &lightPathVertex,
      __in Scene const &scene) const = 0;

  virtual bool NextScatteringDirection(
      __inout PathVertex *const pEyePathVertex,
      __inout IPrimarySample &sample) const = 0;

  virtual float3_t CalculateWeightedDirectLighting(
      __inout PathVertex *const pEyePathVertex,
      __in PathVertex const &lightPathVertex,
      __in ImportanceQuery const &importons, __in Scene const &scene) const;

  // 光源からの重要度を使って反射ベクトルを計算する
  virtual bool NextScatteringDirection(__inout PathVertex *const pEyePathVertex,
                                       __in ImportanceQuery const &importons,
                                       __inout IPrimarySample &sample) const;

  // 重要度計算用
  virtual bool BoundsImportance(__inout PathVertex *const pLightPathVertex,
                                __inout std::mt19937_64 &random) const = 0;

  // 反射光線の分割数
  virtual float_t GetDensityVariance() const = 0;

  //
  // 双方向経路追跡法用
  //
  virtual bool LightsToCameraScatteringDirection(
      __inout PathVertex *const pPathVertex,
      __inout IPrimarySample &sample) const = 0;

  virtual bool CameraToLightsScatteringDirection(
      __inout PathVertex *const pPathVertex,
      __inout IPrimarySample &sample) const = 0;

  virtual bool CameraToLightsScatteringDirection(
      __inout PathVertex *const pPathVertex,
      __in ImportanceQuery const &importons,
      __inout IPrimarySample &sample) const;

  //
  // 双方向経路追跡用
  //
  virtual bool NextScatteringDirection(PathVertex &,      ///< [in/out]
                                       float3_t const &,  ///< [in]
                                       bool const &       ///< [in]
                                       ) const = 0;

  //
  // まだ使ってないけど実装しておく
  //
  // Get the density with respect to projected solid angle: p(w)
  // If p'(w) is the density with respect to ordinary solid angle, then p(w) =
  // p'(w) / |cos\theta|,
  //   where \theta is the angle between the outgoing direction and the surface
  //   normal.
  virtual bool HasNoDiracDistribution() const = 0;
  virtual bool IsLightEmittor() const = 0;

  virtual float_t BSDF(std::size_t const &is,  ///< [in] index of spectrum
                       float3_t const &wi,     ///< [in] incoming direction
                       float3_t const &wo,     ///< [in] outgoing direction
                       float3_t const &ns      ///< [in] shading normal
                       ) const = 0;
};

}  // end of namespace xyz

/*
視点から光源方向へのサンプリング(fSamplingPrev)では，
シェーディング法線と入射方向ベクトルの内積(fIncomingCosThetaShading)で，
サンプリング分布を除す

光源から視点方向へのサンプリング(fSamplingNext)では，
幾何学的法線と散乱方向ベクトルの内積(fOutgoingCosThetaGeometric)で，
サンプリング分布を除す
*/

#endif XYZ_BSDF_BSDF_HPP_
