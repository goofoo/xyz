#include "scene.hpp"
#include "../bsdf/bsdf.hpp"
#include "../geom/triangle.hpp"
#include "../geom/intersection.hpp"
#include "../geom/pathvertex.hpp"
using namespace xyz;

void Scene::EvaluateImportanceDrivenBidirectionalPathTracing(
    __in float_t const x, __in float_t const y,
    __out std::vector<pixel_descriptor_t> *const pColors,
    __inout std::vector<PathVertex> &lightsPathVertices,  // from lights
    __inout std::vector<PathVertex> &theEyePathVertices,  // from the eye
    __inout IPrimarySample &sample) const {
  (*pColors)[0].second = CIE_XYZ_Color(0);  // 初期化

  // 光源上の位置をサンプリングする．
  SampleLightPosition(sample.next(), sample.next(), &lightsPathVertices[1]);
  lightsPathVertices[1].power = GetLightPowerXYZ();
  lightsPathVertices[2].pGeometry = nullptr;

  // フィルム位置をサンプリングする．
  SampleLensPosition(&theEyePathVertices[1]);

  // 光源から視点へ向けて経路を生成する．
  BuildLightsPath(lightsPathVertices, sample);
  EvaluateImportanceDrivenPathS1(lightsPathVertices, theEyePathVertices[1],
                                 pColors, sample);

  // 視点から光源へ向けて経路を生成する．
  BuildImportanceDrivenTheEyePath(x, y, theEyePathVertices, sample);

  // 双方の経路を接続して評価する．
  EvaluateImportanceDrivenPathST(lightsPathVertices, theEyePathVertices,
                                 &((*pColors)[0].second), sample);
}

void Scene::EvaluateImportanceDrivenPathS1(
    __inout std::vector<PathVertex> &lightsPathVertices,
    __in PathVertex const &theEyePathVertex,
    __out std::vector<pixel_descriptor_t> *const pColors,
    __inout IPrimarySample &sample) const {
  UNREFERENCED_PARAMETER(sample);

  Intersection param;

  // 放射束から放射輝度への変換するための係数の定数部分
  float_t const fFluxToRadianceCoefficient =
      camera_.FluxToRadianceCoefficient();

  // (t=1)の場合の各頂点について処理する．
  for (std::size_t s = 2; s < XYZ_CONFIG_kMaxRandomWalkDepth; ++s) {
    PathVertex const &theLightPathVertex = lightsPathVertices[s];

    // 交点がなければ終了する．
    if (!theLightPathVertex.pGeometry) {
      break;
    }

    // specular surface の場合は処理をスキップする．
    Bsdf::type const theLightPathVertexType =
        bsdfs_[theLightPathVertex.pGeometry->bsdf()]->What();
    if ((Bsdf::MIRROR == theLightPathVertexType) ||
        (Bsdf::GLASS == theLightPathVertexType)) {
      continue;
    }

    // 視点位置から光源の交点位置へ向かうベクトルを求める．
    float3_t const vTheEyeToLightsDirection = hi::normalize(
        theLightPathVertex.vPosition - theEyePathVertex.vPosition);

    // 幾何学的法線と散乱方向ベクトルの内積を求める．
    float_t const fOutgoingCosThetaGeometric =
        -hi::dot(vTheEyeToLightsDirection, theLightPathVertex.vGeometricNormal);
    if (fOutgoingCosThetaGeometric <= 0) {
      continue;
    }

    // シェーディング法線と散乱方向ベクトルの内積を求める．
    float_t const fOutgoingCosThetaShading =
        -hi::dot(vTheEyeToLightsDirection, theLightPathVertex.vShadingNormal);
    if (fOutgoingCosThetaShading <= 0) {
      continue;
    }

    // 幾何学的法線とシェーディング法線と入射方向ベクトルの内積を求める．
    float_t const fIncomingCosThetaShading =
        hi::dot(vTheEyeToLightsDirection, theEyePathVertex.vShadingNormal);
    // float_t const fIncomingCosThetaGeometric = fIncomingCosThetaShading;
    if (fIncomingCosThetaShading <=
        0)  // カメラでは vGeometricNormal == vShadingNormal
    {
      continue;
    }

    // ピンホールカメラモデルにおけるフィルム上の位置を求める．
    float2_t vFilmPosition;
    if (!camera_.GetFilmPosition(
            theLightPathVertex.vPosition - theEyePathVertex.vPosition,
            &vFilmPosition)) {
      continue;
    }

    // 可視判定を行って遮蔽物がないかチェックする．
    if (FindIntersection(theEyePathVertex.vPosition, vTheEyeToLightsDirection,
                         &param) != theLightPathVertex.pGeometry) {
      continue;
    }

    // 画素の位置を求める．
    {
      std::size_t const jx =
          static_cast<std::size_t>(vFilmPosition[0] * GetWidth());
      std::size_t const jy =
          static_cast<std::size_t>(vFilmPosition[1] * GetHeight());
      std::size_t const j = jy * GetWidth() + jx;
      (*pColors)[s - 1].first = j;  // 正しい画素の位置を設定する．
    }

    // 放射束から放射輝度への変換係数を求める
    float_t const fFluxToRadianceFactor =
        fFluxToRadianceCoefficient /
        hi::fourth_power_of(fIncomingCosThetaShading);

    // 幾何学項を求める
    float_t const fGeometricFactor = fIncomingCosThetaShading *
                                     fOutgoingCosThetaGeometric /
                                     hi::square_of(param.t_max());

    // 経路長がs+t-1(ただしt=1)の経路を生成
    float_t wst = 1;

    // 視線部分経路を延長していく場合の経路密度を計算
    {
      // 視点からこの経路頂点の方向をサンプリングする確率を|rCosIn|で割ったもの
      hi::single_stack<float_t> stack(lightsPathVertices[s + 1].fSamplingPrev,
                                      fFluxToRadianceFactor);

      float_t pst = 1;
      for (std::size_t i = s; i > 0; --i) {
        pst *= lightsPathVertices[i + 1].fSamplingPrev /
               lightsPathVertices[i - 1].fSamplingNext;

        if (!lightsPathVertices[i].bSpecular &&
            !lightsPathVertices[i - 1].bSpecular) {
          wst += hi::square_of(pst * fGeometricFactor /
                               lightsPathVertices[i].fGeometricFactor);
        }
      }
    }
    wst = hi::rcp(wst);

    (*pColors)[s - 1].second =
        theLightPathVertex.power  // NOTE: ランバート面しか扱わないので Lambert
        // クラスでの処理のおかげでつじつまが合っている．
        *
        theEyePathVertex.power *
        (fFluxToRadianceFactor * fGeometricFactor * M_1_PI *
         wst);  // NOTE: M_1_PI はランバート面を考えているため．
  }
}

// 視点から光源に向けて経路を生成する．
void Scene::BuildImportanceDrivenTheEyePath(
    __in float_t const x, __in float_t const y,
    __inout std::vector<PathVertex> &theEyePathVertices,
    __inout IPrimarySample &sample) const {
  // 初期光線の設定
  {
    camera_.GetPrimaryRayDirection(x, y,
                                   &(theEyePathVertices[1].vIncomingDirection));

    theEyePathVertices[1].power = 1;
    theEyePathVertices[1].fIncomingCosThetaShading =
        hi::dot(theEyePathVertices[1].vIncomingDirection,
                theEyePathVertices[1].vShadingNormal);
    theEyePathVertices[1].fOutgoingCosThetaGeometric = 1;  // NOTE: 仮の値
    theEyePathVertices[1].fSamplingPrev =
        camera_.GetConstFactor() /
        hi::fourth_power_of(theEyePathVertices[1].fIncomingCosThetaShading);
    theEyePathVertices[1].fSamplingNext = 1;     // NOTE: 仮の値
    theEyePathVertices[1].fGeometricFactor = 1;  // NOTE: 仮の値
    theEyePathVertices[1].fBSDFxIPDF = 1;
  }

  for (std::size_t t = 1;;)  // 再帰的に光線を追跡(末尾再帰)
  {
    PathVertex const &oldTheEyePathVertex = theEyePathVertices[t];
    PathVertex &newTheEyePathVertex = theEyePathVertices[t + 1];

    // 交差判定
    Intersection param;
    newTheEyePathVertex.pGeometry =
        FindIntersection(oldTheEyePathVertex.vPosition,
                         oldTheEyePathVertex.vIncomingDirection, &param);
    if (!newTheEyePathVertex.pGeometry) {
      return;
    }

    // 交点の座標とその点の基底ベクトルを計算
    {
      // シェーディング法線を基準に基底ベクトルを計算する．
      newTheEyePathVertex.SetShadingBasis(param);
      newTheEyePathVertex.bBackSide = param.is_back_side();

      // 散乱方向ベクトルをコピーする．
      newTheEyePathVertex.vOutgoingDirection =
          -oldTheEyePathVertex.vIncomingDirection;

      // シェーディング法線の裏から当たった場合は終了する．
      float_t const fOutgoingCosThetaShading =
          hi::dot(newTheEyePathVertex.vOutgoingDirection,
                  newTheEyePathVertex.vShadingNormal);
      if (fOutgoingCosThetaShading <= 0) {
        newTheEyePathVertex.pGeometry = nullptr;
        return;
      }

      // 輸送されたエネルギーをコピーする．
      newTheEyePathVertex.power = oldTheEyePathVertex.power;
      newTheEyePathVertex.power *= oldTheEyePathVertex.fBSDFxIPDF;

      // 交点座標を設定する．
      newTheEyePathVertex.vPosition = oldTheEyePathVertex.vIncomingDirection;
      newTheEyePathVertex.vPosition *= param.t_max();
      newTheEyePathVertex.vPosition += oldTheEyePathVertex.vPosition;

      // 幾何学的法線と散乱方向ベクトルの内積を求める．
      newTheEyePathVertex.fOutgoingCosThetaGeometric =
          hi::dot(newTheEyePathVertex.vOutgoingDirection,
                  newTheEyePathVertex.vGeometricNormal);

      // 幾何学項を計算する．
      newTheEyePathVertex.fGeometricFactor =
          newTheEyePathVertex.fOutgoingCosThetaGeometric *
          oldTheEyePathVertex.fIncomingCosThetaShading /
          hi::square_of(param.t_max());
    }

    Bsdf const *const pBSDF = bsdfs_[newTheEyePathVertex.pGeometry->bsdf()];

    // NOTE: 反射回数の限界数を超えたら打ち切る．
    if (++t >= XYZ_CONFIG_kMaxRandomWalkDepth) {
      theEyePathVertices[t + 1].pGeometry = nullptr;
      return;
    }

    // 光源が直接可視の場合の処理．
    Bsdf::type const theEyePathVertexType = pBSDF->What();
    if (Bsdf::LIGHT == pBSDF->What()) {
      theEyePathVertices[t + 1].pGeometry = nullptr;
      return;
    }

    // 散乱ベクトルの計算
    if ((t > 2) && (Bsdf::LAMBERT == theEyePathVertexType) &&
        (oldTheEyePathVertex.bSpecular))  // デルタ反射後の集光模様を描く
    {
      ImportanceQuery importons(newTheEyePathVertex);
      Scene::GetInstance().ImportonMap().kNN_query(
          newTheEyePathVertex.vPosition, importons);
      if (!pBSDF->CameraToLightsScatteringDirection(&newTheEyePathVertex,
                                                    importons, sample)) {
        theEyePathVertices[t + 1].pGeometry = nullptr;
        return;
      }
    } else {
      if (!pBSDF->CameraToLightsScatteringDirection(&newTheEyePathVertex,
                                                    sample)) {
        theEyePathVertices[t + 1].pGeometry = nullptr;
        return;
      }
    }
  }
}

void Scene::EvaluateImportanceDrivenPathST(
    __inout std::vector<PathVertex> &lightsPathVertices,
    __inout std::vector<PathVertex> &theEyePathVertices,
    __inout CIE_XYZ_Color *const pColor, __inout IPrimarySample &sample) const {
  UNREFERENCED_PARAMETER(sample);

  // (t>1)の場合
  Intersection param;
  for (std::size_t t = 2; t < XYZ_CONFIG_kMaxRandomWalkDepth; ++t) {
    PathVertex &theEyePathVertex = theEyePathVertices[t];

    // 交点が何もなければ終了
    if (!theEyePathVertex.pGeometry) {
      break;
    }

    Bsdf::type const theEyePathVertexType =
        bsdfs_[theEyePathVertex.pGeometry->bsdf()]->What();

    if (Bsdf::LIGHT == theEyePathVertexType) {
      // 経路長がs+t-1(ただしs=0)の経路を生成
      float_t wst = 1;

      if (t > 2) {
        // 光源部分経路を延長していく場合の経路密度を計算する．
        // 光源に衝突すると経路が終了するので上書きしてよい．
        theEyePathVertices[t + 1].fSamplingNext =
            lightsPathVertices[0].fSamplingNext;  // P_{A}(x_{1})
        theEyePathVertex.fSamplingNext = M_1_PI;  // 平面光源なので
        {
          float_t pst = 1;
          for (std::size_t i = t; i > 1; --i) {
            pst *= theEyePathVertices[i + 1].fSamplingNext /
                   theEyePathVertices[i - 1].fSamplingPrev;

            if (!theEyePathVertices[i].bSpecular &&
                !theEyePathVertices[i - 1].bSpecular) {
              wst +=
                  hi::square_of(pst / theEyePathVertices[i].fGeometricFactor);
            }
          }
        }
        wst = hi::rcp(wst);
      }

      *pColor += lightsPathVertices[1].power * theEyePathVertex.power *
                 (hi::rcp(GetLightArea()) * M_1_PI *
                  wst);  // NOTE: ランバート反射面のみ扱う．

      break;
    }

    // specular surface の場合は処理をスキップする
    if ((Bsdf::MIRROR == theEyePathVertexType) ||
        (Bsdf::GLASS == theEyePathVertexType)) {
      continue;
    }

    // 各頂点について処理する．
    for (std::size_t s = 1; (s + t - 1) < XYZ_CONFIG_kMaxRandomWalkDepth; ++s) {
      PathVertex const &theLightPathVertex = lightsPathVertices[s];

      // 交点がなければ終了する．
      if (!theLightPathVertex.pGeometry) {
        break;
      }

      // specular surface の場合は処理をスキップする．
      Bsdf::type const theLightPathVertexType =
          bsdfs_[theLightPathVertex.pGeometry->bsdf()]->What();
      if ((Bsdf::MIRROR == theLightPathVertexType) ||
          (Bsdf::GLASS == theLightPathVertexType)) {
        continue;
      }

      // explicit sampling 用の経路を生成してみる．
      float3_t const vLightsToTheEyeDirection = hi::normalize(
          theEyePathVertex.vPosition - theLightPathVertex.vPosition);

      // 幾何学的法線と散乱方向ベクトルの内積を求める．
      float_t const fOutgoingCosThetaGeometric = hi::dot(
          vLightsToTheEyeDirection, theLightPathVertex.vGeometricNormal);
      if (fOutgoingCosThetaGeometric <= 0) {
        continue;
      }

      // シェーディング法線と散乱方向ベクトルの内積を求める．
      float_t const fOutgoingCosThetaShading =
          hi::dot(vLightsToTheEyeDirection, theLightPathVertex.vShadingNormal);
      if (fOutgoingCosThetaShading <= 0) {
        continue;
      }

      // 幾何学的法線と入射方向ベクトルの内積を求める．
      float_t const fIncomingCosThetaGeometric =
          -hi::dot(vLightsToTheEyeDirection, theEyePathVertex.vGeometricNormal);
      if (fIncomingCosThetaGeometric <= 0) {
        continue;
      }

      // シェーディング法線と入射方向ベクトルの内積を求める．
      float_t const fIncomingCosThetaShading =
          -hi::dot(vLightsToTheEyeDirection, theEyePathVertex.vShadingNormal);
      if (fIncomingCosThetaShading <= 0) {
        continue;
      }

      // 交差判定を行い二頂点間に遮蔽物がないか確認する．
      if (FindIntersection(theLightPathVertex.vPosition,
                           vLightsToTheEyeDirection,
                           &param) != theEyePathVertex.pGeometry) {
        continue;
      }

      // 幾何学項を計算する．
      float_t const fGeometricFactor = fIncomingCosThetaShading *
                                       fOutgoingCosThetaGeometric /
                                       hi::square_of(param.t_max());

      // 経路長がs+t-1の経路を生成する．
      float_t wst = 1;  // この経路自身の重み

      // 視点から光源方向へ経路を延長していく場合の経路密度を計算する．
      // つまり，光源から視点方向へ生成した経路を逆方向にたどる．
      {
        // 現在の視点経路頂点から光源経路頂点をサンプリングする密度を設定する．
        hi::single_stack<float_t> stack(
            lightsPathVertices[s + 1].fSamplingPrev,
            M_1_PI);  // TODO: 要修正(これはランバート面の場合だけ)

        // 相対的な重みを蓄えていく．
        float_t pst = 1;
        for (std::size_t i = s; i > 0; --i) {
          pst *= lightsPathVertices[i + 1].fSamplingPrev /
                 lightsPathVertices[i - 1]
                     .fSamplingNext;  // 境界条件を後でチェックする

          if (!lightsPathVertices[i].bSpecular &&
              !lightsPathVertices[i - 1].bSpecular) {
            wst += hi::square_of(pst * fGeometricFactor /
                                 lightsPathVertices[i].fGeometricFactor);
          }
        }
      }

      // 光源部分経路を延長していく場合の経路密度を計算
      {
        hi::single_stack<float_t> stack(theEyePathVertices[t + 1].fSamplingNext,
                                        M_1_PI);  // TODO: ランバート面なので

        // 相対的な重みを蓄えていく．
        float_t pst = 1;
        for (std::size_t i = t; i > 1;
             --i)  // NOTE: ピンホールカメラモデルなのでiは2まで
        {
          pst *= theEyePathVertices[i + 1].fSamplingNext /
                 theEyePathVertices[i - 1].fSamplingPrev;

          if (!theEyePathVertices[i].bSpecular &&
              !theEyePathVertices[i - 1].bSpecular) {
            wst += hi::square_of(pst * fGeometricFactor /
                                 theEyePathVertices[i].fGeometricFactor);
          }
        }
      }
      wst = hi::rcp(wst);

      *pColor += theLightPathVertex.power  // NOTE:
                 // ランバート面しか扱わないのでつじつまが合っている．
                 *
                 theEyePathVertex.power *
                 (fGeometricFactor * M_1_PI * M_1_PI *
                  wst);  // TODO: M_1_PI はランバート面だけしか扱わないため．
    }
  }
}
