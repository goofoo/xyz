#ifndef __TGIR_BSDF_BSDF_HPP__
#define __TGIR_BSDF_BSDF_HPP__

#include "core/PhotonMap.hpp"
#include <random>

namespace tgir {
class PathVertex;
class PrimarySample;
class Scene;

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
  virtual bool BoundsImportance(tgir::PathVertex &pathVertex,
                                std::mt19937_64 &random) const;
  /// フォトンを追跡する
  virtual bool BoundsPhoton(tgir::PathVertex &lightPathVertex,
                            std::mt19937_64 &random) const;
  /// 視点からの重要度を使ってフォトンを追跡する
  virtual bool BoundsPhotonWithImportonMap(tgir::PathVertex &lightPathVertex,
                                           tgir::ImportonMap const &importonMap,
                                           tgir::ImportonQuery &importons,
                                           std::mt19937_64 &random) const;
  /// 粒子フィルタにより視点からの重要度を使ってフォトンを追跡する
  virtual bool BoundsPhotonWithParticleFilter(
      tgir::PathVertex &lightPathVertex, tgir::ImportonMap const &importonMap,
      tgir::ParticleQuery &importons, std::mt19937_64 &random) const;

  //
  // 経路追跡法用
  //
  virtual Bsdf::type What() const = 0;

  virtual tgir::Spectrum CalculateWeightedDirectLighting(
      tgir::PathVertex &eyePathVertex, tgir::PathVertex const &lightPathVertex,
      std::size_t const uWavelengthIndex, tgir::Scene const &scene) const = 0;

  virtual bool GetScatterVector(tgir::PathVertex &eyePathVertex,
                                tgir::Real const wavelength,
                                std::mt19937_64 &random) const = 0;

  // 光源からの重要度を使って反射ベクトルを計算する
  virtual bool GetScatterVector(tgir::PathVertex &eyePathVertex,
                                tgir::Real const wavelength,
                                tgir::ImportonMap const &importonMap,
                                tgir::ImportanceQuery &importons,
                                std::mt19937_64 &random) const;

  virtual tgir::Real GetDensityVariance() const = 0;

  //
  // 双方向経路追跡用
  //
  virtual bool GetScatterVector(tgir::PathVertex &,           ///< [in/out]
                                tgir::PrimarySample const &,  ///< [in]
                                tgir::Vector3 const &,        ///< [in]
                                bool const &                  ///< [in]
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

  virtual tgir::Real BSDF(std::size_t const &is,    ///< [in] index of spectrum
                          tgir::Vector3 const &wi,  ///< [in] incoming direction
                          tgir::Vector3 const &wo,  ///< [in] outgoing direction
                          tgir::Vector3 const &ns   ///< [in] shading normal
                          ) const = 0;
};

}  // end of namespace tgir

#endif __TGIR_BSDF_BSDF_HPP__
