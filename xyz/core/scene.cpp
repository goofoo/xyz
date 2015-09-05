#include "scene.hpp"
#include "../bsdf/bsdf.hpp"

namespace xyz {
Scene::Scene()
    : displayList_(~0U)
#ifdef LIGHT_PROBE
      ,
      light_probe_(3 * LIGHT_PROBE * LIGHT_PROBE)
#endif
{
  SetDimension(720, 480);
  Clear();
}

Scene::~Scene() { Clear(); }

void Scene::Clear() {
  std::for_each(bsdfs_.begin(), bsdfs_.end(), hi::delete_ptr<Bsdf>());
  bsdfs_.clear();

  kdtree_.Clear();

  std::for_each(triangles_.begin(), triangles_.end(),
                hi::delete_ptr<Triangle>());
  triangles_.clear();

  lights_.clear();
  lightsCdf_.clear();
  lightsCdf_.push_back(0);

  rSurfaceArea_ = 0;
}

void Scene::Build() { kdtree_.Build(triangles_); }

void Scene::Render() const {
  camera_.Render();
  ::glCallList(displayList_);
}

void Scene::SampleLightPosition(
    __in float_t const s,                     ///<
    __in float_t const t,                     ///<
    __out PathVertex *const pLightPathVertex  ///< 光源頂点(s_{1})
    ) const {
  if (lights_.empty()) {
    ::_ftprintf_s(stderr, _TEXT("光源は設定されていません\n"));
    ::_ftprintf_s(stderr, _TEXT(__FUNCTION__) _TEXT("\n"));
    pLightPathVertex->pGeometry = nullptr;
    return;
  }

  float_t const u = s * lightsCdf_.back();

  typedef std::vector<float_t>::const_iterator const_iterator;
  const_iterator const it =
      std::upper_bound(lightsCdf_.begin(), lightsCdf_.end(), u);

  assert(it != lightsCdf_.begin());
  assert(it != lightsCdf_.end());

  std::iterator_traits<const_iterator>::difference_type const n =
      it - lightsCdf_.begin() - 1;
  float_t const v = (u - lightsCdf_[n]) / (lightsCdf_[n + 1] - lightsCdf_[n]);

  pLightPathVertex->pGeometry = lights_[n];
  pLightPathVertex->pGeometry->Sample(
      v, t, &pLightPathVertex->vPosition, &pLightPathVertex->vShadingNormal,
      &pLightPathVertex->vGeometricNormal, &pLightPathVertex->vTangent,
      &pLightPathVertex->vBinormal);

  pLightPathVertex->power = 0;  // 後で設定する
  pLightPathVertex->fSamplingNext = hi::rcp(GetLightArea());
  // pLightPathVertex->fSamplingNext = 1;

  pLightPathVertex->bSpecular = false;
}

void Scene::SampleLensPosition(
    __out PathVertex *const pEyePathVertex  ///< [out]
    ) const {
  pEyePathVertex->pGeometry = nullptr;

  // X = tangent, Y = binormal, Z = normal と読み替える
  // Z 軸の方向を逆にして，binormal = cross(tangent, normal) の関係を保つ
  camera_.GetPrimaryRayOriginAndCameraBasis(
      &pEyePathVertex->vPosition, &pEyePathVertex->vTangent,
      &pEyePathVertex->vBinormal, &pEyePathVertex->vGeometricNormal);

  pEyePathVertex->vShadingNormal = pEyePathVertex->vGeometricNormal =
      -pEyePathVertex->vGeometricNormal;

  pEyePathVertex->fSamplingNext = 0;  // 後で設定する
  pEyePathVertex->power = 1;          // PT用
}

void Scene::SampleLensPosition(
    __in float_t const &s,                   ///< [in]
    __in float_t const &t,                   ///< [in]
    __out float2_t *const pvLensCoordinate,  ///< [out] レンズ座標
    __out PathVertex *const pEyePathVertex   ///< [out] 視点頂点(t_{1})
    ) const {
  pEyePathVertex->pGeometry = nullptr;

  // X = tangent, Y = binormal, Z = normal と読み替える
  // Z 軸の方向を逆にして，binormal = cross(tangent, normal) の関係を保つ
  camera_.GetPrimaryRayOriginAndCameraBasis(
      s, t, pvLensCoordinate, &pEyePathVertex->vPosition,
      &pEyePathVertex->vTangent, &pEyePathVertex->vBinormal,
      &pEyePathVertex->vGeometricNormal);

  pEyePathVertex->vShadingNormal = pEyePathVertex->vGeometricNormal =
      -pEyePathVertex->vGeometricNormal;

  pEyePathVertex->fSamplingNext = 0;  // 後で設定する
  pEyePathVertex->power = 1;          // PT用
}

Scene &Scene::GetInstance() {
  static Scene instance;
  return instance;
}

}  // end of namespace xyz
