#ifndef XYZ_SCENE_HPP
#define XYZ_SCENE_HPP

#include "../core/config.hpp"
#include "../core/thinlenscamera.hpp"
#include "../core/imagefilm.hpp"
#include "../core/photonmap.hpp"
#include "../geom/kdtree.hpp"
#include "../geom/path.hpp"

#ifndef CONFIG_SPACE_SUBDIVISION
#include "../geom/intersection.hpp"
#include "../geom/triangle.hpp"
#endif CONFIG_SPACE_SUBDIVISION

namespace xyz {
class IPrimarySample;
class PathVertex;
class Bsdf;
class Triangle;

struct Particle {
  PathVertex eyePathVertex;
  std::size_t nPixelIndex;
#ifndef CONFIG_RENDER_CAUSTICS
  bool bDirectVisible;
#endif
};

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
  // (0) Path Tracing; PT
  void EvalutePathTracing(__in float_t const x, __in float_t const y,
                          __out CIE_XYZ_Color *const pColor,
                          __inout IPrimarySample &sample) const;

  // (1) Importance Driven Path Tracing; ID-PT
  void Scene::EvaluteImportanceDrivenPathTracing(
      __in float_t const x, __in float_t const y,
      __out CIE_XYZ_Color *const pColor, __inout IPrimarySample &sample) const;

  // (2) Path Tracing with Go with the Winner Strategy; PT w/ GW strategy
  void EvalutePathTracingWithGoWithTheWinnersStrategy(
      __in float_t const x, __in float_t const y,
      __out CIE_XYZ_Color *const pColor, __inout IPrimarySample &sample) const;

  // (3) Importance Driven Path Tracing with Go with the Winner Strategy; ID-PT
  // w/ GW strategy
  void EvaluteImportanceDrivenPathTracingWithGoWithTheWinnersStrategy(
      __in float_t const x, __in float_t const y,
      __out CIE_XYZ_Color *const pColor, __inout IPrimarySample &sample) const;

  // (4) Bidirectional Path Tracing
  void EvaluateBidirectionalPathTracing(
      __in float_t const x, __in float_t const y,
      __out std::vector<pixel_descriptor_t> *const pColors,
      __inout std::vector<PathVertex> &lightsPathVertices,  // from lights
      __inout std::vector<PathVertex> &theEyePathVertices,  // from the eye
      __inout IPrimarySample &sample) const;

  // (5) Importance Driven Bidirectional Path Tracing
  void EvaluateImportanceDrivenBidirectionalPathTracing(
      __in float_t const x, __in float_t const y,
      __out std::vector<pixel_descriptor_t> *const pColors,
      __inout std::vector<PathVertex> &lightsPathVertices,  // from lights
      __inout std::vector<PathVertex> &theEyePathVertices,  // from the eye
      __inout IPrimarySample &sample) const;

  void Scene::EvaluteResampledPathTracing(
      __in std::size_t const x, __in std::size_t const y,
      __out std::vector<CIE_XYZ_Color> *const pColor,
      __inout std::vector<PathVertex> &lightPathVertices,
      __inout std::vector<Particle> &ray, __inout std::vector<Particle> &tmp,
      __inout std::vector<float_t> &cdf, __inout IPrimarySample &sample) const;

 private:
  // Path Tracing with Go with the Winner Strategy; Recursively
  void EvalutePathTracingWithGoWithTheWinnersStrategyRecursively(
      __in std::size_t const k, __in PathVertex const &lightPathVertex,
      __in PathVertex const &eyePathVertex, __inout CIE_XYZ_Color *const pColor,
      __inout IPrimarySample &sample) const;

  // Importance Driven Path Tracing with Go with the Winner Strategy;
  // Recursively
  void
  EvaluteImportanceDrivenPathTracingWithGoWithTheWinnersStrategyRecursively(
      __in std::size_t const k, __in PathVertex const &lightPathVertex,
      __in PathVertex const &eyePathVertex, __inout CIE_XYZ_Color *const pColor,
      __inout IPrimarySample &sample) const;

  // Bidirectional Path Tracing; Subroutine for building lights path
  void BuildLightsPath(__inout std::vector<PathVertex> &lightsPathVertices,
                       __inout IPrimarySample &sample) const;

  // Bidirectional Path Tracing; Subroutine for building the eye path
  void BuildTheEyePath(__in float_t const x, __in float_t const y,
                       __inout std::vector<PathVertex> &theEyePathVertices,
                       __inout IPrimarySample &sample) const;

  void EvaluatePathS1(__inout std::vector<PathVertex> &lightsPathVertices,
                      __in PathVertex const &theEyePathVertex,
                      __out std::vector<pixel_descriptor_t> *const pColors,
                      __inout IPrimarySample &sample) const;

  void EvaluatePathST(__inout std::vector<PathVertex> &lightsPathVertices,
                      __inout std::vector<PathVertex> &theEyePathVertices,
                      __inout CIE_XYZ_Color *const pColor,
                      __inout IPrimarySample &sample) const;

  // Impotance Driven Bidirectional Path Tracing; Subroutine for building the
  // eye path
  void BuildImportanceDrivenTheEyePath(
      __in float_t const x, __in float_t const y,
      __inout std::vector<PathVertex> &theEyePathVertices,
      __inout IPrimarySample &sample) const;

  void EvaluateImportanceDrivenPathS1(
      __inout std::vector<PathVertex> &lightsPathVertices,
      __in PathVertex const &theEyePathVertex,
      __out std::vector<pixel_descriptor_t> *const pColors,
      __inout IPrimarySample &sample) const;

  void EvaluateImportanceDrivenPathST(
      __inout std::vector<PathVertex> &lightsPathVertices,
      __inout std::vector<PathVertex> &theEyePathVertices,
      __inout CIE_XYZ_Color *const pColor,
      __inout IPrimarySample &sample) const;

  /*
      void EvaluteParticleFilter(
        std::size_t const x, std::size_t const y,
        std::vector<float3_t> &,
        std::vector<ParticleInfo> &,
        std::vector<Particle> &,
        std::vector<Particle> &,
        std::vector<float_t> &,
        std::mt19937_64 &) const;

      void EvaluteParticleFilter(
        std::size_t const x, std::size_t const y,
        float3_t & rad,
        std::vector<ParticleInfo> &,
        std::vector<Particle> &,
        std::vector<Particle> &,
        std::vector<float_t> &,
        std::mt19937_64 &) const;
  */
 public:
  // フォトンマップの構築
  void CreateViewImportanceMap(std::mt19937_64 &random);
  void CreatePhotonMap(std::mt19937_64 &random);
  void CreatePhotonMapWithImportonMap(std::mt19937_64 &random);
  void CreatePhotonMapWithParticleFilter(std::mt19937_64 &random);
  void CreateLightImportanceMap(std::mt19937_64 &random);

 private:
  void SampleLightPosition(
      __in float_t const s,                     ///<
      __in float_t const t,                     ///<
      __out PathVertex *const pLightPathVertex  ///< 光源頂点(s_{1})
      ) const;

 public:
  void SampleLensPosition(
      __out PathVertex *const pEyePathVertex  ///< 視点頂点(t_{1})
      ) const;

  void SampleLensPosition(
      __in float_t const &s,                   ///< [in]
      __in float_t const &t,                   ///< [in]
      __out float2_t *const pvLensCoordinate,  ///< [out] レンズ座標
      __out PathVertex *const pEyePathVertex   ///< [out] 視点頂点(t_{1})
      ) const;

 private:
  void TraceLightPath(__in PrimarySample const &sample,
                      __out Path *const pPaths) const;

  void TraceEyePath(__in PrimarySample const &sample,
                    __out Path *const pPaths) const;

  void CombinePaths(__out Path *const pPaths) const;

 public:
  void Render() const;

  inline std::vector<Bsdf const *> const &BSDFs() const { return bsdfs_; }

  inline ThinLensCamera &ThinLensCamera() { return camera_; }

  inline xyz::ThinLensCamera const &ThinLensCamera() const { return camera_; }

#if 0
    // プロパティの定義
    __declspec(property(get=GetWidth)) int Width;
    __declspec(property(get=GetHeight)) int Height;
#endif
  inline int GetWidth() const { return width_; }
  inline int GetHeight() const { return height_; }

  inline float_t GetLightArea() const { return lightsCdf_.back(); }

  inline void SetDimension(int const width, int const height) {
    width_ = width;
    height_ = height;
    camera_.SetAspectRatio(static_cast<float_t>(width) / height);
  }

  // 交差判定
  // public & スレッドセーフじゃないので触る場合は注意
  mutable hi::ulong n_FindIntersection;

  inline Triangle const *FindIntersection(float3_t const &vOrigin,
                                          float3_t const &vDirection,
                                          Intersection *const pParam) const;

 private:
  // materials
  std::vector<Bsdf const *> bsdfs_;

  // geometry
  KDTree kdtree_;
  std::vector<Triangle const *> triangles_;

  // preview
  unsigned int displayList_;

  // light system
  std::vector<Triangle const *> lights_;
  std::vector<float_t> lightsCdf_;
  float_t rSurfaceArea_;

  // camera setting
  xyz::ThinLensCamera camera_;
  int width_;
  int height_;

#ifdef LIGHT_PROBE
  // HDR
  std::vector<float> light_probe_;
#endif LIGHT_PROBE

#ifndef HI_REMOVE_PHOTON_MAPS
  ImportonMap importanceMap_;
  PhotonMap globalPhotonMap_;

 public:
  inline ImportonMap const &ImportonMap() const { return importanceMap_; }

  inline PhotonMap const &GlobalPhotonMap() const { return globalPhotonMap_; }
#endif HI_REMOVE_PHOTON_MAPS
};

inline Triangle const *Scene::FindIntersection(
    float3_t const &vOrigin, float3_t const &vDirection,
    Intersection *const pParam) const {
  ++n_FindIntersection;
#ifdef CONFIG_SPACE_SUBDIVISION
  return kdtree_.FindIntersection(vOrigin, vDirection, pParam);
#else
  pParam->Init();

  Triangle const *s = nullptr;
  for (std::size_t i = 0, size = triangles_.size(); i < size; ++i) {
    if (triangles_[i]->FindIntersection(vOrigin, vDirection, pParam)) {
      s = triangles_[i];
    }
  }

  return s;
#endif CONFIG_SPACE_SUBDIVISION
}

}  // end of namespace xyz

#endif XYZ_SCENE_HPP
