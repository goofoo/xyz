#include "Scene.hpp"
#include "bsdf/Bsdf.hpp"
#include "geom/Triangle.hpp"
#include "geom/Intersection.hpp"
#include "sampler/PrimarySample.hpp"
#include "geom/Path.hpp"

#include "bsdf/Lambert.hpp"  // これを消せるようにする

/*
  この経路を生成する戦略数について
  戦略数 = (経路長+2) - 1
                      ~~~
  大きさのある光源に衝突する場合と，
  大きさのあるレンズに衝突する場合の２を足す．
  ピンホールカメラモデルなので，
  後者の場合がなくなり１引く．
*/

namespace tgir {
void Scene::Evalute(__in tgir::PrimarySample const &sample,
                    __out tgir::Path *const pPaths,
                    __in EMutationType type) const {
  pPaths->Begin(type);
  {
    if (type & MUTATION_EYE) {
#ifdef CONFIG_PINHOLE_CAMERA_MODEL
      SampleLensPosition(&pPaths->EyePathVertices[1]);
#else
      tgir::Vector3 const &u = sample.E(1);
      SampleLensPosition(u[0], u[1], &pPaths->vLensCoordinate,
                         &pPaths->EyePathVertices[1]);
#endif CONFIG_PINHOLE_CAMERA_MODEL
    }
    if (type & MUTATION_LIGHT) {
      tgir::Vector3 const &u = sample.L(0);
      SampleLightPosition(u[0], u[1], &pPaths->LightPathVertices[1]);
      pPaths->LightPathVertices[1].sQuantum =
          tgir::GetLightPower()[sample.WavelengthIndex()];
      pPaths->LightPathVertices[2].pGeometry = nullptr;
      TraceLightPath(sample, pPaths);
    }
    if (type & MUTATION_EYE) {
      TraceEyePath(sample, pPaths);
    }
    CombinePaths(sample, pPaths);
  }
  pPaths->End();
}

void Scene::TraceLightPath(__in tgir::PrimarySample const &sample,
                           __out tgir::Path *const pPaths) const {
  if (!pPaths->LightPathVertices[1].pGeometry) {
    return;
  }

  tgir::PathVertex &theEyePathVertex = pPaths->EyePathVertices[1];

  // 放射束から放射輝度への変換係数における定数部分
  tgir::Real const fConstFactor = camera_.GetConstFactor();

  // 射出方向の選択
  {
    tgir::PathVertex &theLightPathVertex = pPaths->LightPathVertices[1];

    tgir::ComputeDiffusedVector(sample.L(1), theLightPathVertex.vTangent,
                                theLightPathVertex.vGeometricNormal,
                                theLightPathVertex.vBinormal,
                                &theLightPathVertex.vOutgoingDirection);
  }

  for (std::size_t s = 1;;) {
    PathVertex const &oldLightPathVertex = pPaths->LightPathVertices[s];
    PathVertex &newLightPathVertex = pPaths->LightPathVertices[s + 1];

    // 交差判定
    tgir::Intersection param;
    newLightPathVertex.pGeometry =
        FindIntersection(oldLightPathVertex.vPosition,
                         oldLightPathVertex.vOutgoingDirection, &param);
    if (!newLightPathVertex.pGeometry) {
      break;  // シーン外へ
    }

    // 交点の座標とその点の基底ベクトルを計算
    {
      newLightPathVertex.SetGeometricBasis(
          param);  // 幾何法線を基準に基底ベクトルを計算する
      newLightPathVertex.back_side = param.is_back_side();

      if (-hi::dot(oldLightPathVertex.vOutgoingDirection,
                   newLightPathVertex.vShadingNormal) <= 0) {
        newLightPathVertex.pGeometry = nullptr;
        return;
      }

      newLightPathVertex.sQuantum = oldLightPathVertex.sQuantum;
      newLightPathVertex.vOutgoingDirection =
          oldLightPathVertex.vOutgoingDirection;
      newLightPathVertex.vPosition = oldLightPathVertex.vOutgoingDirection;
      newLightPathVertex.vPosition *= param.t_max();
      newLightPathVertex.vPosition += oldLightPathVertex.vPosition;
      newLightPathVertex.rGeometricFactor =
          hi::dot(oldLightPathVertex.vOutgoingDirection,
                  oldLightPathVertex.vGeometricNormal) *
          hi::dot(oldLightPathVertex.vOutgoingDirection,
                  newLightPathVertex.vShadingNormal) /
          hi::square_of(param.t_max());
      newLightPathVertex.rGeometricFactor =
          std::abs(newLightPathVertex.rGeometricFactor);  // 非負値化
      newLightPathVertex.bSpecular = false;
    }

    //[TODO]glossyとかshinyなマテリアルを実装したときに、改めて設計を考え直す。
    tgir::Bsdf const &bsdf = *bsdfs_[newLightPathVertex.pGeometry->bsdf()];

    switch (bsdf.What()) {
    case tgir::Bsdf::LAMBERT: {
      newLightPathVertex.sQuantum *=
          hi::dot(oldLightPathVertex.vOutgoingDirection,
                  newLightPathVertex.vShadingNormal) /
          hi::dot(oldLightPathVertex.vOutgoingDirection,
                  newLightPathVertex.vGeometricNormal);  // modify BSDF
      newLightPathVertex.sQuantum *=
          hi::polymorphic_downcast<tgir::Lambert const *>(&bsdf)
              ->Albedo(sample.WavelengthIndex());  //

      tgir::Vector3 const vIncomingDirection = hi::normalize(
          newLightPathVertex.vPosition - theEyePathVertex.vPosition);

      tgir::Real const rCosOutDash =
          -hi::dot(vIncomingDirection, newLightPathVertex.vGeometricNormal);
      if ((rCosOutDash <= 0) ||
          (-hi::dot(vIncomingDirection, newLightPathVertex.vShadingNormal) <=
           0)) {
        break;
      }

      tgir::Real const rCosIn =
          hi::dot(vIncomingDirection, theEyePathVertex.vShadingNormal);
      if (rCosIn <= 0)  // カメラでは vGeometricNormal == vShadingNormal
      {
        break;
      }

      // フィルム上の位置を求める
      tgir::Vector2 vFilmPosition;
#ifdef CONFIG_PINHOLE_CAMERA_MODEL
      if (!camera_.GetFilmPosition(
              newLightPathVertex.vPosition - theEyePathVertex.vPosition,
              vFilmPosition))
#else
      if (!camera_.GetFilmPosition(
              pPaths->vLensCoordinate,
              newLightPathVertex.vPosition - theEyePathVertex.vPosition,
              &vFilmPosition))
#endif CONFIG_PINHOLE_CAMERA_MODEL
      {
        break;
      }

// 可視判定
#ifdef CONFIG_BACKFACE_CULLING
      if (FindIntersectionCW(theEyePathVertex.vPosition, vIncomingDirection,
                             param) != newLightPathVertex.pGeometry)
#else
      if (FindIntersection(theEyePathVertex.vPosition, vIncomingDirection,
                           &param) != newLightPathVertex.pGeometry)
#endif CONFIG_BACKFACE_CULLING
      {
        break;
      }

      // 画素の位置を求める
      std::size_t const x =
          static_cast<std::size_t>(vFilmPosition[0] * GetWidth());
      std::size_t const y =
          static_cast<std::size_t>(vFilmPosition[1] * GetHeight());

      // 放射束から放射輝度への変換係数を求める
      tgir::Real const rFluxToRadianceFactor =
          fConstFactor / hi::fourth_power_of(rCosIn);

      // 幾何学項を求める
      tgir::Real const rGeometricFactor =
          rCosOutDash * rCosIn / hi::square_of(param.t_max());

      // 経路長 s+1 の経路を生成
      tgir::Real wst = 1;

      // 視線部分経路を延長していく場合の経路密度を計算
      {
        // 視点からこの経路頂点の方向をサンプリングする確率を|rCosIn|で割ったもの
        pPaths->LightPathVertices[s + 2].rSamplingPrev = rFluxToRadianceFactor;
        tgir::Real pst = 1;
        for (std::size_t i = s + 1; i > 0; --i) {
          pst *= pPaths->LightPathVertices[i + 1].rSamplingPrev /
                 pPaths->LightPathVertices[i - 1].rSamplingNext;

          if (!pPaths->LightPathVertices[i].bSpecular &&
              !pPaths->LightPathVertices[i - 1].bSpecular) {
            wst += hi::square_of(pst * rGeometricFactor /
                                 pPaths->LightPathVertices[i].rGeometricFactor);
          }
        }
      }
      wst = hi::rcp(wst);

      tgir::Real const cst =
          newLightPathVertex.sQuantum * theEyePathVertex.sQuantum *
          (rFluxToRadianceFactor * rGeometricFactor * M_1_PI);

      pPaths->AccumulateS1(y * GetWidth() + x, sample.WavelengthIndex(), s,
                           wst * cst);
    } break;
    case tgir::Bsdf::LIGHT:
      newLightPathVertex.pGeometry = nullptr;
      return;
    }

    if (++s >= pPaths->GetMaxPathLength()) {
      return;  // 打ち止め(エネルギーは消失する)
    }

    // 反射
    if (!bsdf.GetScatterVector(newLightPathVertex, sample, sample.L(s), true)) {
      pPaths->LightPathVertices[s].pGeometry = nullptr;
      return;
    }
  }
}

void Scene::TraceEyePath(tgir::PrimarySample const &sample,
                         tgir::Path *const pPaths) const {
  // 初期光線の設定
  {
    tgir::PathVertex &theEyePathVertex = pPaths->EyePathVertices[1];

    pPaths->SetFilmPosition(sample.GetPixelPosition(GetWidth(), GetHeight()));
#ifdef CONFIG_PINHOLE_CAMERA_MODEL
    camera_.GetPrimaryRayDirection(sample.E(0)[0], sample.E(0)[1],
                                   &theEyePathVertex.vOutgoingDirection);
#else
    camera_.GetPrimaryRayDirection(sample.E(0)[0], sample.E(0)[1],
                                   pPaths->vLensCoordinate,
                                   &theEyePathVertex.vOutgoingDirection);
#endif CONFIG_PINHOLE_CAMERA_MODEL
    tgir::Real const rCosOutDash = hi::dot(theEyePathVertex.vOutgoingDirection,
                                           theEyePathVertex.vShadingNormal);
    theEyePathVertex.rSamplingNext =
        camera_.GetConstFactor() / hi::fourth_power_of(rCosOutDash);
  }

  enum {
    IGNORE_STATE,
    INITIAL_STATE,
    TRACE1_STATE,
    BOUNDS_STATE,
    TRACE2_STATE,
  } pathState = INITIAL_STATE;

  for (std::size_t t = 1;;)  // 再帰的に光線を追跡(末尾再帰)
  {
    PathVertex const &oldEyePathVertex = pPaths->EyePathVertices[t];
    PathVertex &newEyePathVertex = pPaths->EyePathVertices[t + 1];

    // 交差判定
    tgir::Intersection param;

#ifdef CONFIG_BACKFACE_CULLING
    newEyePathVertex.pGeometry =
        (t == 1)
            ? FindIntersectionCW(oldEyePathVertex.vPosition,
                                 oldEyePathVertex.vOutgoingDirection, param)
            : FindIntersection(oldEyePathVertex.vPosition,
                               oldEyePathVertex.vOutgoingDirection, param);
#else
    newEyePathVertex.pGeometry =
        FindIntersection(oldEyePathVertex.vPosition,
                         oldEyePathVertex.vOutgoingDirection, &param);
#endif CONFIG_BACKFACE_CULLING

    if (!newEyePathVertex.pGeometry) {
      return;
    }

    // 交点の座標とその点の基底ベクトルを計算
    {
      newEyePathVertex.SetShadingBasis(
          param);  // シェーディング法線を基準に基底ベクトルを計算する
      newEyePathVertex.back_side = param.is_back_side();

      if (-hi::dot(oldEyePathVertex.vOutgoingDirection,
                   newEyePathVertex.vShadingNormal) <= 0) {
        newEyePathVertex.pGeometry = nullptr;
        return;
      }

      newEyePathVertex.sQuantum = oldEyePathVertex.sQuantum;
      newEyePathVertex.vOutgoingDirection =
          oldEyePathVertex.vOutgoingDirection;  // set a incoming direction
      newEyePathVertex.vPosition = oldEyePathVertex.vOutgoingDirection;
      newEyePathVertex.vPosition *= param.t_max();
      newEyePathVertex.vPosition += oldEyePathVertex.vPosition;
      newEyePathVertex.rGeometricFactor =
          hi::dot(oldEyePathVertex.vOutgoingDirection,
                  oldEyePathVertex.vShadingNormal) *
          hi::dot(oldEyePathVertex.vOutgoingDirection,
                  newEyePathVertex.vGeometricNormal) /
          hi::square_of(param.t_max());
      newEyePathVertex.rGeometricFactor =
          std::abs(newEyePathVertex.rGeometricFactor);  // 非負値化
      newEyePathVertex.bSpecular = false;
    }

    // 陰影処理
    tgir::Bsdf const &bsdf = *bsdfs_[newEyePathVertex.pGeometry->bsdf()];

    switch (bsdf.What()) {
    case tgir::Bsdf::LAMBERT:
      newEyePathVertex.sQuantum *=
          hi::polymorphic_downcast<tgir::Lambert const *>(&bsdf)
              ->Albedo(sample.WavelengthIndex());
      break;
    case tgir::Bsdf::LIGHT:
#ifndef CONFIG_DOUBLE_LIGHT
      if (!param.is_back_side())
#endif CONFIG_DOUBLE_LIGHT
      {
        // 経路長 t の経路を生成
        tgir::Real wst = 1;

        // 光源部分経路を延長していく場合の経路密度を計算
        if (t > 1)  // 光源のサンプル点は直接可視化しないので
        {
          pPaths->EyePathVertices[t + 2].rSamplingPrev =
              pPaths->LightPathVertices[1].rSamplingNext;  // P_{A}(x_{1})
          {
            tgir::Real pst = 1;
            for (std::size_t i = t + 1; i > 1; --i) {
              pst *= pPaths->EyePathVertices[i + 1].rSamplingPrev /
                     pPaths->EyePathVertices[i - 1].rSamplingNext;

              if (!pPaths->EyePathVertices[i].bSpecular &&
                  !pPaths->EyePathVertices[i - 1].bSpecular) {
                wst += hi::square_of(
                    pst / pPaths->EyePathVertices[i].rGeometricFactor);
              }
            }
          }
        }
        wst = hi::rcp(wst);

        tgir::Real const cst = pPaths->LightPathVertices[1].sQuantum *
                               newEyePathVertex.sQuantum *
                               (M_1_PI / GetLightArea());  // power -> radiance

        pPaths->Accumulate0T(sample.WavelengthIndex(), wst * cst,
                             TRACE2_STATE == pathState);
      }

      newEyePathVertex.pGeometry = nullptr;
      return;
    }

    if (++t >= pPaths->GetMaxPathLength()) {
      pPaths->EyePathVertices[t + 1].pGeometry = nullptr;
      return;  // 経路長VNの鏡面反射からの寄与は計算しない。
    }

// 散乱ベクトルの計算
#ifdef CONFIG_PINHOLE_CAMERA_MODEL
    if (!bsdf.GetScatterVector(newEyePathVertex, sample, sample.E(t - 1),
                               false))
#else
    if (!bsdf.GetScatterVector(newEyePathVertex, sample, sample.E(t), false))
#endif CONFIG_PINHOLE_CAMERA_MODEL
    {
      pPaths->EyePathVertices[t + 1].pGeometry = nullptr;
      return;
    }

    /*
    switch (pathState)
    {
    case INITIAL_STATE: pathState = newEyePathVertex.bSpecular ? TRACE1_STATE :
  IGNORE_STATE; break; // ES
    case TRACE1_STATE : if (!newEyePathVertex.bSpecular) pathState =
  BOUNDS_STATE; break;          // S*D
    case BOUNDS_STATE : pathState = newEyePathVertex.bSpecular ? TRACE2_STATE :
  IGNORE_STATE; break; // S
  //case TRACE2_STATE : if (!newEyePathVertex.bSpecular) pathState =
  IGNORE_STATE; break;          // S*
    }
    */
    if (newEyePathVertex.bSpecular) pathState = TRACE2_STATE;
  }
}

void Scene::CombinePaths(__in tgir::PrimarySample const &sample,
                         __out tgir::Path *const pPaths) const {
  tgir::Intersection param;

  for (std::size_t t = 1; t < pPaths->GetMaxPathLength(); ++t) {
    PathVertex &theEyePathVertex = pPaths->EyePathVertices[t + 1];

    if (!theEyePathVertex.pGeometry) {
      break;
    }

    if (tgir::Bsdf::LAMBERT !=
        bsdfs_[theEyePathVertex.pGeometry->bsdf()]->What()) {
      continue;
    }

    for (std::size_t s = 0; (s + t) < pPaths->GetMaxPathLength(); ++s) {
      PathVertex const &theLightPathVertex = pPaths->LightPathVertices[s + 1];

      if (!theLightPathVertex.pGeometry) {
        break;
      }
      if (theLightPathVertex.bSpecular) {
        continue;
      }

      tgir::Vector3 const vIncomingDirection = hi::normalize(
          theEyePathVertex.vPosition - theLightPathVertex.vPosition);

      tgir::Real const rCosOutDash =
          hi::dot(vIncomingDirection, theLightPathVertex.vGeometricNormal);
      if ((rCosOutDash <= 0) ||
          (hi::dot(vIncomingDirection, theLightPathVertex.vShadingNormal) <=
           0)) {
        continue;
      }

      tgir::Real const rCosIn =
          -hi::dot(vIncomingDirection, theEyePathVertex.vShadingNormal);
      if ((rCosIn <= 0) || (-hi::dot(vIncomingDirection,
                                     theEyePathVertex.vGeometricNormal) <= 0)) {
        continue;
      }

      if (FindIntersection(theLightPathVertex.vPosition, vIncomingDirection,
                           &param) != theEyePathVertex.pGeometry) {
        continue;
      }

      tgir::Real const rGeometricFactor =
          rCosOutDash * rCosIn / hi::square_of(param.t_max());

      // 経路長 t+s+1 の経路を生成
      tgir::Real wst = 1;

      // 視線部分経路を延長していく場合の経路密度を計算
      {
        hi::single_stack<tgir::Real> stack(
            pPaths->LightPathVertices[s + 2].rSamplingPrev,
            theEyePathVertex.rSamplingNext);
        tgir::Real pst = 1;
        for (std::size_t i = s + 1; i > 0; --i) {
          pst *= pPaths->LightPathVertices[i + 1].rSamplingPrev /
                 pPaths->LightPathVertices[i - 1].rSamplingNext;

          if (!pPaths->LightPathVertices[i].bSpecular &&
              !pPaths->LightPathVertices[i - 1].bSpecular) {
            wst += hi::square_of(pst * rGeometricFactor /
                                 pPaths->LightPathVertices[i].rGeometricFactor);
          }
        }
      }

      // 光源部分経路を延長していく場合の経路密度を計算
      {
        hi::single_stack<tgir::Real> stack(
            pPaths->EyePathVertices[t + 2].rSamplingPrev,
            theLightPathVertex.rSamplingNext);
        tgir::Real pst = 1;
        for (std::size_t i = t + 1; i > 1; --i) {
          pst *= pPaths->EyePathVertices[i + 1].rSamplingPrev /
                 pPaths->EyePathVertices[i - 1].rSamplingNext;

          if (!pPaths->EyePathVertices[i].bSpecular &&
              !pPaths->EyePathVertices[i - 1].bSpecular) {
            wst += hi::square_of(pst * rGeometricFactor /
                                 pPaths->EyePathVertices[i].rGeometricFactor);
          }
        }
      }

      wst = hi::rcp(wst);

      tgir::Real const cst = theLightPathVertex.sQuantum *
                             theEyePathVertex.sQuantum *
                             (rGeometricFactor * (M_1_PI * M_1_PI));

      if (s > 0) {
        pPaths->AccumulateST(sample.WavelengthIndex(), wst * cst);
      } else {
        pPaths->Accumulate1T(sample.WavelengthIndex(), wst * cst);
      }
    }
  }
}

}  // end of namespace tgir
