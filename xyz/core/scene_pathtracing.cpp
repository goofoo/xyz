#include "scene.hpp"
#include "../bsdf/bsdf.hpp"
#include "../geom/triangle.hpp"
#include "../geom/intersection.hpp"
#include "../geom/pathvertex.hpp"
using namespace xyz;

/// <summary>
/// path tracing with multiple importance sampling
/// </summary>
void Scene::EvalutePathTracing(__in float_t const x, __in float_t const y,
                               __out CIE_XYZ_Color *const pColor,
                               __inout IPrimarySample &sample) const {
  PathVertex eyePathVertex;
  PathVertex lightPathVertex;

  // 光源位置のサンプリング
  SampleLightPosition(sample.next(), sample.next(), &lightPathVertex);
  lightPathVertex.power = GetLightPowerXYZ();

#ifdef CONFIG_PINHOLE_CAMERA_MODEL
  SampleLensPosition(&eyePathVertex);
#else
  float2_t vLensCoordinate;
  SampleLensPosition(sample.next(), sample.next(), &vLensCoordinate,
                     &eyePathVertex);
#endif CONFIG_PINHOLE_CAMERA_MODEL

  // 初期光線の設定
  {
#ifdef CONFIG_PINHOLE_CAMERA_MODEL
    camera_.GetPrimaryRayDirection(x, y, &eyePathVertex.vOutgoingDirection);
#else
    camera_.GetPrimaryRayDirection(x, y, vLensCoordinate,
                                   &eyePathVertex.vOutgoingDirection);
#endif CONFIG_PINHOLE_CAMERA_MODEL
    eyePathVertex.bSpecular = true;
  }

#ifndef CONFIG_RENDER_CAUSTICS
  // 光源に到達した場合，視点から直接見えているかどうかのフラグ
  bool bDirectVisible = true;
#endif

  *pColor = CIE_XYZ_Color(CIE_XYZ_Color::value_type(0.0));
  for (std::size_t k = 1;;) {  // 再帰的に光線を追跡(末尾再帰)
    // 交差判定
    Intersection param;
    eyePathVertex.pGeometry = FindIntersection(
        eyePathVertex.vPosition, eyePathVertex.vOutgoingDirection, &param);

    if (!eyePathVertex.pGeometry) {
      break;  // through the air
    }

    // 交点とその点の基底ベクトルを計算
    {
      eyePathVertex.SetShadingBasis(param);  // 基底
      eyePathVertex.bBackSide = param.is_back_side();

      // シェーディング法線の裏からあたった場合は，素直に終了
      if (-hi::dot(eyePathVertex.vOutgoingDirection,
                   eyePathVertex.vShadingNormal) <= 0.0) {
        break;
      }

      eyePathVertex.vPosition +=
          eyePathVertex.vOutgoingDirection * param.t_max();  // 交点
    }

#ifndef CONFIG_RENDER_CAUSTICS
    if (!eyePathVertex.bSpecular) {
      bDirectVisible = false;
    }
#endif

    // 陰影処理
    Bsdf const &bsdf = *bsdfs_[eyePathVertex.pGeometry->bsdf()];

// 直接光源が見えている場合
#ifdef CONFIG_DOUBLE_LIGHT
    if (bsdf.What() == Bsdf::LIGHT)
#else
    if ((bsdf.What() == Bsdf::LIGHT) && !eyePathVertex.bBackSide)
#endif
    {
      if (eyePathVertex.bSpecular) {
        // 鏡面反射から光源に到達する経路は，これ一つしか存在しないので重み付けなくてよい
        *pColor += lightPathVertex.power * eyePathVertex.power *
                   (M_1_PI / GetLightArea());
      } else {
        float_t const fCosTheta = std::abs(hi::dot(
            eyePathVertex.vOutgoingDirection, eyePathVertex.vShadingNormal));
        float_t const fDistanceSquared = hi::square_of(param.t_max());

        // 暗黙的サンプリングの経路密度
        float_t const fDensityOfImplicit =
            eyePathVertex.fSamplingNext * fCosTheta / fDistanceSquared;

        // 明示的サンプリングの経路密度
        float_t const fDensityOfExplicit = hi::rcp(GetLightArea());

        // 経路密度の比
        float_t const fDensity = fDensityOfExplicit / fDensityOfImplicit;

        // 重みとして，光源方向のサンプリング確率を除し，光源上の点のサンプリング確率を掛ける
        *pColor += lightPathVertex.power * eyePathVertex.power *
                   (M_1_PI / GetLightArea()) *
                   hi::rcp(1 + hi::square_of(fDensity));  // 重み
      }
      break;
    }

    // NOTE: 反射回数の限界数を超えたら打ち切る
    if (++k >= XYZ_CONFIG_kMaxRandomWalkDepth) {
      break;
    }

    // 明示的サンプリング: Direct Lighting Calculation
    {
      *pColor += bsdf.CalculateWeightedDirectLighting(&eyePathVertex,
                                                      lightPathVertex, *this);
    }

    // 暗黙的サンプリング: Indirect Lighting Calculation
    {
      // 確率的に追跡を打ち切る: Russian-Roulette
      {
        float_t const t = eyePathVertex.power[1];
        if (t < CONFIG_TERMINATE_THRESHOLD) {
          float_t const q = t / CONFIG_TERMINATE_THRESHOLD;
          if (sample.next() >= q) {  // Random Termination
            break;                   // Terminate the Path
          }
          eyePathVertex.power *= hi::rcp(q);
        }
      }

      if (!bsdf.NextScatteringDirection(&eyePathVertex, sample)) {
        break;
      }
    }
  }
}

/// <summary>
/// path tracing with multiple importance sampling & go with the winners
/// strategy
/// </summary>
void Scene::EvalutePathTracingWithGoWithTheWinnersStrategy(
    __in float_t const x, __in float_t const y, __out CIE_XYZ_Color *pColor,
    __inout IPrimarySample &sample) const {
  PathVertex lightPathVertex;
  PathVertex eyePathVertex;

  // 光源位置のサンプリング
  SampleLightPosition(sample.next(), sample.next(), &lightPathVertex);
  lightPathVertex.power = GetLightPowerXYZ();

#ifdef CONFIG_PINHOLE_CAMERA_MODEL
  SampleLensPosition(&eyePathVertex);
#else
  float2_t vLensCoordinate;
  SampleLensPosition(sample.next(), sample.next(), &vLensCoordinate,
                     &eyePathVertex);
#endif CONFIG_PINHOLE_CAMERA_MODEL

  // 初期光線の設定
  {
#ifdef CONFIG_PINHOLE_CAMERA_MODEL
    camera_.GetPrimaryRayDirection(x, y, &eyePathVertex.vOutgoingDirection);
#else
    camera_.GetPrimaryRayDirection(x, y, vLensCoordinate,
                                   &eyePathVertex.vOutgoingDirection);
#endif CONFIG_PINHOLE_CAMERA_MODEL
    eyePathVertex.bSpecular = true;
  }

  *pColor = CIE_XYZ_Color(CIE_XYZ_Color::value_type(0));
  this->EvalutePathTracingWithGoWithTheWinnersStrategyRecursively(
      0, lightPathVertex, eyePathVertex, pColor, sample);
}

void Scene::EvalutePathTracingWithGoWithTheWinnersStrategyRecursively(
    __in std::size_t const k, __in PathVertex const &lightPathVertex,
    __in PathVertex const &eyePathVertex, __inout CIE_XYZ_Color *const pColor,
    __inout IPrimarySample &sample) const {
  PathVertex newEyePathVertex(eyePathVertex);

  // 交差判定
  Intersection param;
  newEyePathVertex.pGeometry = FindIntersection(
      newEyePathVertex.vPosition, newEyePathVertex.vOutgoingDirection, &param);

  if (!newEyePathVertex.pGeometry) {
    return;  // through the air
  }

  // 交点とその点の基底ベクトルを計算
  {
    newEyePathVertex.SetShadingBasis(param);  // 基底
    newEyePathVertex.bBackSide = param.is_back_side();

    // シェーディング法線の裏からあたった場合は，素直に終了
    if (-hi::dot(newEyePathVertex.vOutgoingDirection,
                 newEyePathVertex.vShadingNormal) <= 0) {
      return;
    }

    newEyePathVertex.vPosition +=
        newEyePathVertex.vOutgoingDirection * param.t_max();  // 交点
  }

  // 陰影処理
  Bsdf const &bsdf = *bsdfs_[newEyePathVertex.pGeometry->bsdf()];

// 直接光源が見えている場合
#ifdef CONFIG_DOUBLE_LIGHT
  if (bsdf.What() == Bsdf::LIGHT)
#else
  if ((bsdf.What() == Bsdf::LIGHT) && !newEyePathVertex.bBackSide)
#endif
  {
    if (newEyePathVertex.bSpecular) {
      // 鏡面反射から光源に到達する経路は，これ一つしか存在しないので重み付けなくてよい
      *pColor += lightPathVertex.power * newEyePathVertex.power *
                 (M_1_PI / GetLightArea());
    } else {
      float_t const fCosTheta =
          std::abs(hi::dot(newEyePathVertex.vOutgoingDirection,
                           newEyePathVertex.vShadingNormal));
      float_t const fDistanceSquared = hi::square_of(param.t_max());

      // 暗黙的サンプリングの経路密度
      float_t const fDensityOfImplicit =
          newEyePathVertex.fSamplingNext * fCosTheta / fDistanceSquared;

      // 明示的サンプリングの経路密度
      float_t const fDensityOfExplicit = hi::rcp(GetLightArea());

      // 経路密度の比
      float_t const fDensity = fDensityOfExplicit / fDensityOfImplicit;

      // 重みとして，光源方向のサンプリング確率を除し，光源上の点のサンプリング確率を掛ける
      *pColor += lightPathVertex.power * newEyePathVertex.power *
                 (M_1_PI / GetLightArea()) *
                 hi::rcp(1 + hi::square_of(fDensity));  // 重み
    }
    return;
  }

  // NOTE: 反射回数の限界数を超えたら打ち切る
  if ((k + 1) >= XYZ_CONFIG_kMaxRandomWalkDepth) {
    return;
  }

  // 明示的サンプリング: Direct Lighting Calculation
  {
    *pColor += bsdf.CalculateWeightedDirectLighting(&newEyePathVertex,
                                                    lightPathVertex, *this);
  }

#if 0
  // 確率的に追跡を打ち切る
  // Russian-Roulette
  {
    float_t const t = newEyePathVertex.power[1];
    if (t < CONFIG_TERMINATE_THRESHOLD) {
      float_t const q = t / CONFIG_TERMINATE_THRESHOLD;
      if (sample.next() >= q) {// Random Termination
        return; // Terminate the Path
      }
      newEyePathVertex.power *= hi::rcp(q);
    }
  }
#endif

  // 暗黙的サンプリング: Indirect Lighting Calculation
  {
    // 反射本数を決定: Go with the Winners strategy
    float_t const n = newEyePathVertex.power[1] * bsdf.GetDensityVariance();
    std::size_t m = static_cast<std::size_t>(n);
    if (sample.next() < (n - m)) {
      ++m;
    }

    if (m > 0) {
      float3_t const vOutgoingDirection(newEyePathVertex.vOutgoingDirection);
      float3_t const power(newEyePathVertex.power / n);

      for (std::size_t l = 0; l < m; ++l) {
        // この二つのパラメータは新しく計算するために前の状態から計算する必要がある．
        newEyePathVertex.vOutgoingDirection = vOutgoingDirection;
        newEyePathVertex.power = power;

        if (bsdf.NextScatteringDirection(&newEyePathVertex, sample)) {
          this->EvalutePathTracingWithGoWithTheWinnersStrategyRecursively(
              (k + 1), lightPathVertex, newEyePathVertex, pColor, sample);
        }
      }
    }
  }
}
