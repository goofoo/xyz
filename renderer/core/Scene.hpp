#ifndef __TGIR_CORE_SCENE_HPP
#define __TGIR_CORE_SCENE_HPP

#include "core/config.hpp"
#include "core/Camera.hpp"
#include "core/Film.hpp"
#include "core/PhotonMap.hpp"
#include "geom/KDTree.hpp"
#include "geom/Path.hpp"

namespace tgir {
class PathVertex;
class Bsdf;
class Triangle;

struct Particle {
  std::size_t uInfoIndex;
  tgir::PathVertex eyePathVertex;
#ifndef CONFIG_RENDER_CAUSTICS
  bool bDirectVisible;
#endif CONFIG_RENDER_CAUSTICS
};

struct ParticleInfo {
  tgir::Spectrum sRadiance;
  std::pair<tgir::Real, std::size_t> wavelength;
  tgir::PathVertex lightPathVertex;
};

#ifndef HI_REMOVE_PHOTON_MAPS
struct Dvpm : public hi::thread {
  void run();
};

struct Dvim : public hi::thread {
  void run();
};

struct Dvpf : public hi::thread {
  void run();
};
#endif HI_REMOVE_PHOTON_MAPS

enum EMutationType {
  MUTATION_LIGHT = 1,
  MUTATION_EYE = 2,
  MUTATION_FULL = 3,
};

//
// シーン (Singleton)
//
class Scene {
 private:
  Scene();
  Scene(Scene const &rhs);
  Scene operator=(Scene const &rhs);

 public:
  static Scene &GetInstance();

  ~Scene();

  void Clear();
  void Load(std::tistream &in);
  void Build();

 public:
  void Evalute(__in tgir::PrimarySample const &sample,
               __out tgir::Path *const pPaths,
               __in EMutationType type = MUTATION_FULL) const;

  void EvaluteRussianRoulette(std::size_t const x, std::size_t const y,
                              tgir::SpectrumVector &, std::mt19937_64 &) const;
  /*
      void EvaluteGoWithTheWinner(
        std::size_t const x, std::size_t const y,
        tgir::SpectrumVector &,
        std::mt19937_64 &) const;
  */
  void EvaluteParticleFilter(std::size_t const x, std::size_t const y,
                             std::vector<tgir::SpectrumVector> &,
                             std::vector<tgir::ParticleInfo> &,
                             std::vector<tgir::Particle> &,
                             std::vector<tgir::Particle> &,
                             std::vector<tgir::Real> &,
                             std::mt19937_64 &) const;

  void EvaluteParticleFilter(std::size_t const x, std::size_t const y,
                             tgir::SpectrumVector &rad,
                             std::vector<tgir::ParticleInfo> &,
                             std::vector<tgir::Particle> &,
                             std::vector<tgir::Particle> &,
                             std::vector<tgir::Real> &,
                             std::mt19937_64 &) const;
  /*
      void EvaluteImportanceDriven(
        std::size_t const x,
        std::size_t const y,
        tgir::SpectrumVector & sColor,
        tgir::PathVertex & eyePathVertex,
        tgir::ImportanceQuery & importans,
        std::mt19937_64 & random) const;
  */
 public:
  // フォトンマップの構築
  void CreateViewImportanceMap(std::mt19937_64 &random);
  void CreatePhotonMap(std::mt19937_64 &random);
  void CreatePhotonMapWithImportonMap(std::mt19937_64 &random);
  void CreatePhotonMapWithParticleFilter(std::mt19937_64 &random);
  void CreateLightImportanceMap(std::mt19937_64 &random);

 public:
  void SampleLightPosition(
      __in tgir::Real const &s,                       ///<
      __in tgir::Real const &t,                       ///<
      __out tgir::PathVertex *const pLightPathVertex  ///< 光源頂点(s_{1})
      ) const;

  void SampleLensPosition(
      __out tgir::PathVertex *const pEyePathVertex  ///< 視点頂点(t_{1})
      ) const;

  void SampleLensPosition(
      __in tgir::Real const &s,                     ///< [in]
      __in tgir::Real const &t,                     ///< [in]
      __out tgir::Vector2 *const pvLensCoordinate,  ///< [out] レンズ座標
      __out tgir::PathVertex *const pEyePathVertex  ///< [out] 視点頂点(t_{1})
      ) const;

 private:
  void TraceLightPath(__in tgir::PrimarySample const &sample,
                      __out tgir::Path *const pPaths) const;

  void TraceEyePath(__in tgir::PrimarySample const &sample,
                    __out tgir::Path *const pPaths) const;

  void CombinePaths(__in tgir::PrimarySample const &sample,
                    __out tgir::Path *const pPaths) const;
  /*
      tgir::Spectrum TraceGoWithTheWinner(
        std::size_t, bool,
        std::pair<tgir::Real, std::size_t> const &,
        tgir::PathVertex const &,
        tgir::PathVertex,
        std::mt19937_64 &) const;
  */

 public:
  void Render() const;

  inline std::vector<tgir::Bsdf const *> const &BSDFs() const { return bsdfs_; }

  inline tgir::Camera &Camera() { return camera_; }

  inline tgir::Camera const &Camera() const { return camera_; }

#if 0
    // プロパティの定義
    __declspec(property(get=GetWidth)) int Width;
    __declspec(property(get=GetHeight)) int Height;
#endif
  inline int GetWidth() const { return width_; }
  inline int GetHeight() const { return height_; }

  inline tgir::Real GetLightArea() const { return lightsCdf_.back(); }

  inline void SetDimension(int const width, int const height) {
    width_ = width;
    height_ = height;
    camera_.SetAspectRatio(static_cast<tgir::Real>(width) / height);
  }

  // 交差判定
  // public & スレッドセーフじゃないので触る場合は注意
  mutable hi::ulong n_FindIntersection;

  inline tgir::Triangle const *FindIntersection(
      tgir::Vector3 const &vOrigin, tgir::Vector3 const &vDirection,
      tgir::Intersection *const pParam) const;

 private:
  // materials
  std::vector<tgir::Bsdf const *> bsdfs_;

  // geometry
  tgir::KDTree kdtree_;
  std::vector<tgir::Triangle const *> triangles_;

  // preview
  unsigned int displayList_;

  // light system
  std::vector<tgir::Triangle const *> lights_;
  std::vector<tgir::Real> lightsCdf_;
  tgir::Real rSurfaceArea_;

  // camera setting
  tgir::Camera camera_;
  int width_;
  int height_;

#ifdef LIGHT_PROBE
  // HDR
  std::vector<float> light_probe_;
#endif LIGHT_PROBE

#ifndef HI_REMOVE_PHOTON_MAPS
  tgir::ImportonMap importanceMap_;
  tgir::PhotonMap globalPhotonMap_;

 public:
  inline tgir::ImportonMap const &ImportonMap() const { return importanceMap_; }

  inline tgir::PhotonMap const &GlobalPhotonMap() const {
    return globalPhotonMap_;
  }
#endif HI_REMOVE_PHOTON_MAPS
};

inline tgir::Triangle const *Scene::FindIntersection(
    tgir::Vector3 const &vOrigin, tgir::Vector3 const &vDirection,
    tgir::Intersection *const pParam) const {
  ++n_FindIntersection;
#ifdef CONFIG_SPACE_SUBDIVISION
  return kdtree_.FindIntersection(vOrigin, vDirection, pParam);
#else
  param.Init();

  tgir::Triangle *s = hi::null;
  for (std::size_t i = 0, size = triangles_.size(); i < size; ++i) {
    if (triangles_[i]->FindIntersection(vOrigin, vDirection, pParam)) {
      s = triangles_[i];
    }
  }

  return s;
#endif CONFIG_SPACE_SUBDIVISION
}

}  // end of namespace tgir

#endif __TGIR_CORE_SCENE_HPP
