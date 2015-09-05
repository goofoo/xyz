#include "Scene.hpp"
#include "bsdf/Bsdf.hpp"

namespace tgir {
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
  std::for_each(bsdfs_.begin(), bsdfs_.end(), hi::delete_ptr<tgir::Bsdf>());
  bsdfs_.clear();

  kdtree_.Clear();

  std::for_each(triangles_.begin(), triangles_.end(),
                hi::delete_ptr<tgir::Triangle>());
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
    __in tgir::Real const &s,                       ///<
    __in tgir::Real const &t,                       ///<
    __out tgir::PathVertex *const pLightPathVertex  ///< 光源頂点(s_{1})
    ) const {
  if (lights_.empty()) {
    // std::tcerr << _TEXT("光源は設定されていません") << std::endl;
    // std::tcerr << _TEXT(__FUNCTION__) << std::endl;
    pLightPathVertex->pGeometry = nullptr;
    return;
  }

  tgir::Real const u = s * lightsCdf_.back();

  typedef std::vector<tgir::Real>::const_iterator const_iterator;
  const_iterator const it =
      std::upper_bound(lightsCdf_.begin(), lightsCdf_.end(), u);

  assert(it != lightsCdf_.begin());
  assert(it != lightsCdf_.end());

  std::iterator_traits<const_iterator>::difference_type const n =
      it - lightsCdf_.begin() - 1;
  tgir::Real const v =
      (u - lightsCdf_[n]) / (lightsCdf_[n + 1] - lightsCdf_[n]);

  pLightPathVertex->pGeometry = lights_[n];
  pLightPathVertex->pGeometry->Sample(
      v, t, &pLightPathVertex->vPosition, &pLightPathVertex->vShadingNormal,
      &pLightPathVertex->vGeometricNormal, &pLightPathVertex->vTangent,
      &pLightPathVertex->vBinormal);

  pLightPathVertex->sQuantum = 0;  // 後で設定する
  // pLightPathVertex->rSamplingNext = hi::rcp(GetLightArea());
  pLightPathVertex->rSamplingNext = 1;
}

void Scene::SampleLensPosition(
    __out tgir::PathVertex *const pEyePathVertex  ///< [out]
    ) const {
  pEyePathVertex->pGeometry = nullptr;

  // X = tangent, Y = binormal, Z = normal と読み替える
  // Z 軸の方向を逆にして，binormal = cross(tangent, normal) の関係を保つ
  camera_.GetPrimaryRayOriginAndCameraBasis(
      &pEyePathVertex->vPosition, &pEyePathVertex->vTangent,
      &pEyePathVertex->vBinormal, &pEyePathVertex->vGeometricNormal);

  pEyePathVertex->vShadingNormal = pEyePathVertex->vGeometricNormal =
      -pEyePathVertex->vGeometricNormal;

  pEyePathVertex->rSamplingNext = 0;  // 後で設定する
  pEyePathVertex->sQuantum = 1;       // PT用
}

void Scene::SampleLensPosition(
    __in tgir::Real const &s,                     ///< [in]
    __in tgir::Real const &t,                     ///< [in]
    __out tgir::Vector2 *const pvLensCoordinate,  ///< [out] レンズ座標
    __out tgir::PathVertex *const pEyePathVertex  ///< [out] 視点頂点(t_{1})
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

  pEyePathVertex->rSamplingNext = 0;  // 後で設定する
  pEyePathVertex->sQuantum = 1;       // PT用
}

Scene &Scene::GetInstance() {
  static Scene instance;
  return instance;
}

}  // end of namespace tgir
