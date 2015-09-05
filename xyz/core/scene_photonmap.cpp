#ifndef HI_REMOVE_PHOTON_MAPS

#include "scene.hpp"
#include "../bsdf/bsdf.hpp"
#include "../geom/triangle.hpp"
#include "../geom/intersection.hpp"
#include "../geom/pathvertex.hpp"

namespace {
std::size_t const EMITTED_BUNDLES = CONFIG_EMITED_SIZE / CONFIG_BUNDLE_SIZE;
xyz::float_t const PHOTON_POWER = CONFIG_LIGHT_POWER / CONFIG_EMITED_SIZE;

}  // end of unnamed namespace

namespace xyz {
/// <summary>
/// 視点からの重要度のマップを構築
/// </summary>
void Scene::CreateViewImportanceMap(std::mt19937_64 &random) {
  ::_ftprintf_s(stderr, _TEXT("> CreateViewImportanceMap\n"));
  importanceMap_.clear(CONFIG_IMPORTON_SIZE * 10);

  PathVertex eyePathVertex;

  for (std::size_t y = 0, h = CONFIG_IMPORTON_STRATIFY; y < h; ++y) {
    for (std::size_t x = 0, w = CONFIG_IMPORTON_STRATIFY; x < w; ++x) {
      SampleLensPosition(&eyePathVertex);
      camera_.GetPrimaryRayDirection((x + random.next<float_t>()) / w,
                                     (y + random.next<float_t>()) / h,
                                     &eyePathVertex.vOutgoingDirection);

      for (std::size_t depth = 0; depth < XYZ_CONFIG_kMaxRandomWalkDepth;
           ++depth)  // 再帰的に光線を追跡(末尾再帰)
      {
        // 交差判定
        Intersection param;
        eyePathVertex.pGeometry = FindIntersection(
            eyePathVertex.vPosition, eyePathVertex.vOutgoingDirection, &param);

        if (!eyePathVertex.pGeometry) {
          break;  // through the air
        }

        // 交点とその点の基底ベクトルを計算
        {
          eyePathVertex.SetGeometricBasis(param);  // 基底
          eyePathVertex.bBackSide = param.is_back_side();
          eyePathVertex.vPosition +=
              eyePathVertex.vOutgoingDirection * param.t_max();  // 交点
        }

        Bsdf const *const bsdf = bsdfs_[eyePathVertex.pGeometry->bsdf()];

        // インポータンスマップに追加
        if (bsdf->What() == Bsdf::LAMBERT) {
          importanceMap_.push_back(Importance(
              eyePathVertex.vPosition, eyePathVertex.vOutgoingDirection));
        }

        // 反射処理
        if (!bsdf->BoundsImportance(&eyePathVertex, random)) {
          break;
        }
      }
    }
  }

  ::_ftprintf_s(stderr, _TEXT("> BuildImportanceMap\n"));
  importanceMap_.build();
}

/// <summary>
/// フォトンマップを構築
/// </summary>
void Scene::CreatePhotonMap(std::mt19937_64 &random) {
  ::_ftprintf_s(stderr, _TEXT("> CreatePhotonMap\n"));
  globalPhotonMap_.clear(CONFIG_EMITED_SIZE * 10);

  PathVertex lightPathVertex;

  for (std::size_t i = 0; i < CONFIG_EMITED_SIZE; ++i) {
    SampleLightPosition(random.next<float_t>(), random.next<float_t>(),
                        &lightPathVertex);
    lightPathVertex.power = PHOTON_POWER;

    ComputeDiffusedVector(
        random.next<float_t>(), random.next<float_t>(),
        lightPathVertex.vTangent, lightPathVertex.vGeometricNormal,
        lightPathVertex.vBinormal, &lightPathVertex.vOutgoingDirection);

    for (std::size_t depth = 0; depth < XYZ_CONFIG_kMaxRandomWalkDepth;
         ++depth)  // 再帰的に光線を追跡(末尾再帰)
    {
      // 交差判定
      Intersection param;
      lightPathVertex.pGeometry =
          FindIntersection(lightPathVertex.vPosition,
                           lightPathVertex.vOutgoingDirection, &param);

      if (!lightPathVertex.pGeometry) {
        break;  // through the air
      }

      // 交点とその点の基底ベクトルを計算
      {
        lightPathVertex.SetGeometricBasis(
            param);  // 基底(幾何学法線しか使用しない)
        lightPathVertex.bBackSide = param.is_back_side();
        lightPathVertex.vPosition +=
            lightPathVertex.vOutgoingDirection * param.t_max();  // 交点
      }

      // フォトンマップに追加
      globalPhotonMap_.push_back(Photon(lightPathVertex.vPosition,
                                        lightPathVertex.vOutgoingDirection,
                                        lightPathVertex.power));

      // 反射処理
      if (!bsdfs_[lightPathVertex.pGeometry->bsdf()]->BoundsPhoton(
              lightPathVertex, random)) {
        break;
      }
    }
  }

  ::_ftprintf_s(stderr, _TEXT("> BuildPhotonMap\n"));
  globalPhotonMap_.build();
}

/// <summary>
/// インポータンスマップを使用してフォトンマップを構築
/// </summary>
void Scene::CreatePhotonMapWithImportonMap(std::mt19937_64 &random) {
  ::_ftprintf_s(stderr, _TEXT("> CreatePhotonMap\n"));
  globalPhotonMap_.clear(CONFIG_EMITED_SIZE * 10);

  PathVertex lightPathVertex;
  ImportonQuery importons(lightPathVertex);

  for (std::size_t i = 0; i < CONFIG_EMITED_SIZE; ++i) {
    SampleLightPosition(random.next<float_t>(), random.next<float_t>(),
                        &lightPathVertex);
    lightPathVertex.power = PHOTON_POWER;

    importanceMap_.kNN_query(lightPathVertex.vPosition, importons);
    importons.SetDiffusedVector(random);

    for (std::size_t depth = 0; depth < XYZ_CONFIG_kMaxRandomWalkDepth;
         ++depth)  // 再帰的に光線を追跡(末尾再帰)
    {
      // 交差判定
      Intersection param;
      lightPathVertex.pGeometry =
          FindIntersection(lightPathVertex.vPosition,
                           lightPathVertex.vOutgoingDirection, &param);

      if (!lightPathVertex.pGeometry) {
        break;  // through the air
      }

      // 交点とその点の基底ベクトルを計算
      {
        lightPathVertex.SetGeometricBasis(
            param);  // 基底(幾何学法線しか使用しない)
        lightPathVertex.bBackSide = param.is_back_side();
        lightPathVertex.vPosition +=
            lightPathVertex.vOutgoingDirection * param.t_max();  // 交点
      }

      // フォトンマップに追加
      globalPhotonMap_.push_back(Photon(lightPathVertex.vPosition,
                                        lightPathVertex.vOutgoingDirection,
                                        lightPathVertex.power));

      // 反射処理
      if (!bsdfs_[lightPathVertex.pGeometry->bsdf()]
               ->BoundsPhotonWithImportonMap(lightPathVertex, importanceMap_,
                                             importons, random)) {
        break;
      }
    }
  }

  ::_ftprintf_s(stderr, _TEXT("> BuildPhotonMap\n"));
  globalPhotonMap_.build();
}

/// <summary>
/// インポータンスマップと粒子フィルタを使用してフォトンマップを構築
/// </summary>
void Scene::CreatePhotonMapWithParticleFilter(std::mt19937_64 &random) {
  ::_ftprintf_s(stderr, _TEXT("> CreatePhotonMap\n"));
  globalPhotonMap_.clear(CONFIG_EMITED_SIZE * 10);

  std::vector<PathVertex> lightPathVertices(CONFIG_BUNDLE_SIZE);
  std::vector<PathVertex> tempPathVertices(CONFIG_BUNDLE_SIZE);
  ParticleQuery importons(lightPathVertices, tempPathVertices);

  for (std::size_t i = 0; i < EMITTED_BUNDLES; ++i) {
    for (std::size_t j = 0; j < CONFIG_BUNDLE_SIZE; ++j) {
      SampleLightPosition(random.next<float_t>(), random.next<float_t>(),
                          &lightPathVertices[j]);
      lightPathVertices[j].power = PHOTON_POWER;

      importons.SetVertexIndex(j, j);
      importanceMap_.kNN_query(lightPathVertices[j].vPosition, importons);

      tempPathVertices[j] = lightPathVertices[j];
    }

    importons.SetDiffusedVector(CONFIG_BUNDLE_SIZE, random);

    std::size_t nRayCount = CONFIG_BUNDLE_SIZE;

    for (std::size_t depth = 0; depth < XYZ_CONFIG_kMaxRandomWalkDepth;
         ++depth)  // 再帰的に光線を追跡(末尾再帰)
    {
      std::size_t nDiffuseCount = 0;

      for (std::size_t j = 0; j < nRayCount; ++j) {
        // 交差判定
        Intersection param;
        lightPathVertices[j].pGeometry =
            FindIntersection(lightPathVertices[j].vPosition,
                             lightPathVertices[j].vOutgoingDirection, &param);

        if (!lightPathVertices[j].pGeometry) {
          continue;  // through the air
        }

        // 交点とその点の基底ベクトルを計算
        {
          lightPathVertices[j].SetGeometricBasis(
              param);  // 基底(幾何学法線しか使用しない)
          lightPathVertices[j].bBackSide = param.is_back_side();
          lightPathVertices[j].vPosition +=
              lightPathVertices[j].vOutgoingDirection * param.t_max();  // 交点
        }

        // フォトンマップに追加
        globalPhotonMap_.push_back(
            Photon(lightPathVertices[j].vPosition,
                   lightPathVertices[j].vOutgoingDirection,
                   lightPathVertices[j].power));

        // 反射処理
        importons.SetVertexIndex(j, nDiffuseCount);
        if (!bsdfs_[lightPathVertices[j].pGeometry->bsdf()]
                 ->BoundsPhotonWithParticleFilter(
                     lightPathVertices[j], importanceMap_, importons, random)) {
          continue;
        }

        tempPathVertices[nDiffuseCount++] = lightPathVertices[j];
      }

      if (nDiffuseCount <= 0) {
        break;  // 次の束へ
      }

      importons.SetDiffusedVector(nRayCount = nDiffuseCount, random);
    }
  }

  ::_ftprintf_s(stderr, _TEXT("> BuildPhotonMap\n"));
  globalPhotonMap_.build();
}

/// <summary>
/// 光源からの重要度のマップを構築
/// </summary>
void Scene::CreateLightImportanceMap(std::mt19937_64 &random) {
  ::_ftprintf_s(stderr, _TEXT("> CreateLightImportanceMap\n"));
  importanceMap_.clear(CONFIG_IMPORTON_SIZE * XYZ_CONFIG_kMaxRandomWalkDepth);

  PathVertex lightPathVertex;

  for (std::size_t i = 0; i < CONFIG_EMITED_SIZE; ++i) {
    SampleLightPosition(random.next<float_t>(), random.next<float_t>(),
                        &lightPathVertex);

    ComputeDiffusedVector(
        random.next<float_t>(), random.next<float_t>(),
        lightPathVertex.vTangent, lightPathVertex.vGeometricNormal,
        lightPathVertex.vBinormal, &lightPathVertex.vOutgoingDirection);

    for (std::size_t depth = 0; depth < XYZ_CONFIG_kMaxRandomWalkDepth;
         ++depth)  // 再帰的に光線を追跡(末尾再帰)
    {
      // 交差判定
      Intersection param;
      lightPathVertex.pGeometry =
          this->FindIntersection(lightPathVertex.vPosition,
                                 lightPathVertex.vOutgoingDirection, &param);

      if (!lightPathVertex.pGeometry) {
        break;  // through the air
      }

      Bsdf const *const bsdf = bsdfs_[lightPathVertex.pGeometry->bsdf()];

      // インポータンスマップに格納
      if (bsdf->What() == Bsdf::LAMBERT) {
        // 交点とその点の基底ベクトルを計算
        lightPathVertex.SetGeometricBasis(param);  // 基底
        lightPathVertex.bBackSide = param.is_back_side();
        lightPathVertex.vPosition +=
            lightPathVertex.vOutgoingDirection * param.t_max();  // 交点

        /// TEST: 一度以上の反射を経た鏡面反射後のフォトンしか格納しない
        if ((depth > 0) && (lightPathVertex.bSpecular))
        // if (depth > 0)
        {
          importanceMap_.push_back(Importance(
              lightPathVertex.vPosition, lightPathVertex.vOutgoingDirection));
        }
      } else {
        // 交点とその点の基底ベクトルを計算
        lightPathVertex.SetShadingBasis(param);  // 基底
        lightPathVertex.bBackSide = param.is_back_side();
        lightPathVertex.vPosition +=
            lightPathVertex.vOutgoingDirection * param.t_max();  // 交点
      }

      // 反射処理
      if (!bsdf->BoundsImportance(&lightPathVertex, random)) {
        break;
      }
    }
  }

  ::_ftprintf_s(stderr, _TEXT("> BuildLightImportanceMap\n"));
  importanceMap_.build();
}
}  // end of namespace xyz

#endif  // HI_REMOVE_PHOTON_MAPS
