#ifndef HI_REMOVE_PHOTON_MAPS

#include "Scene.hpp"
#include "bsdf/Bsdf.hpp"
#include "geom/Triangle.hpp"
#include "geom/Intersection.hpp"
#include "geom/PathVertex.hpp"

namespace {
std::size_t const EMITTED_BUNDLES = CONFIG_EMITED_SIZE / CONFIG_BUNDLE_SIZE;
tgir::Real const PHOTON_POWER = CONFIG_LIGHT_POWER / CONFIG_EMITED_SIZE;

}  // end of unnamed namespace

namespace tgir {
/// <summary>
/// 視点からの重要度のマップを構築
/// </summary>
void Scene::CreateViewImportanceMap(std::mt19937_64 &random) {
  std::tcerr << _TEXT("> CreateLightImportanceMap") << std::endl;
  importanceMap_.clear(CONFIG_IMPORTON_SIZE * 10);

  tgir::PathVertex eyePathVertex;

  for (std::size_t y = 0, h = CONFIG_IMPORTON_STRATIFY; y < h; ++y) {
    for (std::size_t x = 0, w = CONFIG_IMPORTON_STRATIFY; x < w; ++x) {
      SampleLensPosition(&eyePathVertex);
      camera_.GetPrimaryRayDirection(
          (x + std::uniform_real_distribution<tgir::Real>()(random)) / w,
          (y + std::uniform_real_distribution<tgir::Real>()(random)) / h,
          &eyePathVertex.vOutgoingDirection);

      for (std::size_t depth = 0; depth < TGIR_CONFIG_kMaxRandomWalkDepth;
           ++depth)  // 再帰的に光線を追跡(末尾再帰)
      {
        // 交差判定
        tgir::Intersection param;
        eyePathVertex.pGeometry = FindIntersection(
            eyePathVertex.vPosition, eyePathVertex.vOutgoingDirection, &param);

        if (!eyePathVertex.pGeometry) {
          break;  // through the air
        }

        // 交点とその点の基底ベクトルを計算
        {
          eyePathVertex.SetGeometricBasis(param);  // 基底
          eyePathVertex.back_side = param.is_back_side();
          eyePathVertex.vPosition +=
              eyePathVertex.vOutgoingDirection * param.t_max();  // 交点
        }

        // インポータンスマップに追加
        if (depth > 0) {
          importanceMap_.push_back(Importance(
              eyePathVertex.vPosition, eyePathVertex.vOutgoingDirection));
        }

        // 反射処理
        if (!bsdfs_[eyePathVertex.pGeometry->bsdf()]->BoundsImportance(
                eyePathVertex, random)) {
          break;
        }
      }
    }
  }

  std::tcerr << _TEXT("> BuildImportanceMap") << std::endl;
  importanceMap_.build();
}

/// <summary>
/// フォトンマップを構築
/// </summary>
void Scene::CreatePhotonMap(std::mt19937_64 &random) {
  std::tcerr << _TEXT("> CreatePhotonMap") << std::endl;
  globalPhotonMap_.clear(CONFIG_EMITED_SIZE * 10);

  tgir::PathVertex lightPathVertex;

  for (std::size_t i = 0; i < CONFIG_EMITED_SIZE; ++i) {
    SampleLightPosition(std::uniform_real_distribution<tgir::Real>()(random),
                        std::uniform_real_distribution<tgir::Real>()(random),
                        &lightPathVertex);
    lightPathVertex.sQuantum = PHOTON_POWER;

    tgir::ComputeDiffusedVector(
        std::uniform_real_distribution<tgir::Real>()(random),
        std::uniform_real_distribution<tgir::Real>()(random),
        lightPathVertex.vTangent, lightPathVertex.vGeometricNormal,
        lightPathVertex.vBinormal, &lightPathVertex.vOutgoingDirection);

    for (std::size_t depth = 0; depth < TGIR_CONFIG_kMaxRandomWalkDepth;
         ++depth)  // 再帰的に光線を追跡(末尾再帰)
    {
      // 交差判定
      tgir::Intersection param;
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
        lightPathVertex.back_side = param.is_back_side();
        lightPathVertex.vPosition +=
            lightPathVertex.vOutgoingDirection * param.t_max();  // 交点
      }

      // フォトンマップに追加
      globalPhotonMap_.push_back(Photon(lightPathVertex.vPosition,
                                        lightPathVertex.vOutgoingDirection,
                                        lightPathVertex.sQuantum));

      // 反射処理
      if (!bsdfs_[lightPathVertex.pGeometry->bsdf()]->BoundsPhoton(
              lightPathVertex, random)) {
        break;
      }
    }
  }

  std::tcerr << _TEXT("> BuildPhotonMap") << std::endl;
  globalPhotonMap_.build();
}

/// <summary>
/// インポータンスマップを使用してフォトンマップを構築
/// </summary>
void Scene::CreatePhotonMapWithImportonMap(std::mt19937_64 &random) {
  std::tcerr << _TEXT("> CreatePhotonMap") << std::endl;
  globalPhotonMap_.clear(CONFIG_EMITED_SIZE * 10);

  tgir::PathVertex lightPathVertex;
  tgir::ImportonQuery importons(lightPathVertex);

  for (std::size_t i = 0; i < CONFIG_EMITED_SIZE; ++i) {
    SampleLightPosition(std::uniform_real_distribution<tgir::Real>()(random),
                        std::uniform_real_distribution<tgir::Real>()(random),
                        &lightPathVertex);
    lightPathVertex.sQuantum = PHOTON_POWER;

    importanceMap_.kNN_query(lightPathVertex.vPosition, importons);
    importons.SetDiffusedVector(random);

    for (std::size_t depth = 0; depth < TGIR_CONFIG_kMaxRandomWalkDepth;
         ++depth)  // 再帰的に光線を追跡(末尾再帰)
    {
      // 交差判定
      tgir::Intersection param;
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
        lightPathVertex.back_side = param.is_back_side();
        lightPathVertex.vPosition +=
            lightPathVertex.vOutgoingDirection * param.t_max();  // 交点
      }

      // フォトンマップに追加
      globalPhotonMap_.push_back(Photon(lightPathVertex.vPosition,
                                        lightPathVertex.vOutgoingDirection,
                                        lightPathVertex.sQuantum));

      // 反射処理
      if (!bsdfs_[lightPathVertex.pGeometry->bsdf()]
               ->BoundsPhotonWithImportonMap(lightPathVertex, importanceMap_,
                                             importons, random)) {
        break;
      }
    }
  }

  std::tcerr << _TEXT("> BuildPhotonMap") << std::endl;
  globalPhotonMap_.build();
}

/// <summary>
/// インポータンスマップと粒子フィルタを使用してフォトンマップを構築
/// </summary>
void Scene::CreatePhotonMapWithParticleFilter(std::mt19937_64 &random) {
  std::tcerr << _TEXT("> CreatePhotonMap") << std::endl;
  globalPhotonMap_.clear(CONFIG_EMITED_SIZE * 10);

  std::vector<tgir::PathVertex> lightPathVertices(CONFIG_BUNDLE_SIZE);
  std::vector<tgir::PathVertex> tempPathVertices(CONFIG_BUNDLE_SIZE);
  tgir::ParticleQuery importons(lightPathVertices, tempPathVertices);

  for (std::size_t i = 0; i < EMITTED_BUNDLES; ++i) {
    for (std::size_t j = 0; j < CONFIG_BUNDLE_SIZE; ++j) {
      SampleLightPosition(std::uniform_real_distribution<tgir::Real>()(random),
                          std::uniform_real_distribution<tgir::Real>()(random),
                          &lightPathVertices[j]);
      lightPathVertices[j].sQuantum = PHOTON_POWER;

      importons.SetVertexIndex(j, j);
      importanceMap_.kNN_query(lightPathVertices[j].vPosition, importons);

      tempPathVertices[j] = lightPathVertices[j];
    }

    importons.SetDiffusedVector(CONFIG_BUNDLE_SIZE, random);

    std::size_t nRayCount = CONFIG_BUNDLE_SIZE;

    for (std::size_t depth = 0; depth < TGIR_CONFIG_kMaxRandomWalkDepth;
         ++depth)  // 再帰的に光線を追跡(末尾再帰)
    {
      std::size_t nDiffuseCount = 0;

      for (std::size_t j = 0; j < nRayCount; ++j) {
        // 交差判定
        tgir::Intersection param;
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
          lightPathVertices[j].back_side = param.is_back_side();
          lightPathVertices[j].vPosition +=
              lightPathVertices[j].vOutgoingDirection * param.t_max();  // 交点
        }

        // フォトンマップに追加
        globalPhotonMap_.push_back(
            Photon(lightPathVertices[j].vPosition,
                   lightPathVertices[j].vOutgoingDirection,
                   lightPathVertices[j].sQuantum));

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

  std::tcerr << _TEXT("> BuildPhotonMap") << std::endl;
  globalPhotonMap_.build();
}

/// <summary>
/// 光源からの重要度のマップを構築
/// </summary>
void Scene::CreateLightImportanceMap(std::mt19937_64 &random) {
  std::tcerr << _TEXT("> CreateLightImportanceMap") << std::endl;
  importanceMap_.clear(CONFIG_IMPORTON_SIZE * 10);

  tgir::PathVertex lightPathVertex;

  for (std::size_t i = 0; i < CONFIG_EMITED_SIZE; ++i) {
    SampleLightPosition(std::uniform_real_distribution<tgir::Real>()(random),
                        std::uniform_real_distribution<tgir::Real>()(random),
                        &lightPathVertex);

    tgir::ComputeDiffusedVector(
        std::uniform_real_distribution<tgir::Real>()(random),
        std::uniform_real_distribution<tgir::Real>()(random),
        lightPathVertex.vTangent, lightPathVertex.vGeometricNormal,
        lightPathVertex.vBinormal, &lightPathVertex.vOutgoingDirection);

    for (std::size_t depth = 0; depth < TGIR_CONFIG_kMaxRandomWalkDepth;
         ++depth)  // 再帰的に光線を追跡(末尾再帰)
    {
      // 交差判定
      tgir::Intersection param;
      lightPathVertex.pGeometry =
          FindIntersection(lightPathVertex.vPosition,
                           lightPathVertex.vOutgoingDirection, &param);

      if (!lightPathVertex.pGeometry) {
        break;  // through the air
      }

      // 交点とその点の基底ベクトルを計算
      {
        lightPathVertex.SetGeometricBasis(param);  // 基底
        lightPathVertex.back_side = param.is_back_side();
        lightPathVertex.vPosition +=
            lightPathVertex.vOutgoingDirection * param.t_max();  // 交点
      }

      // インポータンスマップに追加
      if (depth > 0) {
        importanceMap_.push_back(Importance(
            lightPathVertex.vPosition, lightPathVertex.vOutgoingDirection));
      }

      // 反射処理
      if (!bsdfs_[lightPathVertex.pGeometry->bsdf()]->BoundsImportance(
              lightPathVertex, random)) {
        break;
      }
    }
  }

  std::tcerr << _TEXT("> BuildLightImportanceMap") << std::endl;
  importanceMap_.build();
}
}  // end of namespace tgir

namespace tgir {
void Dvpm::run() {
  std::tcerr << _TEXT("** Photon Mapping (0) **") << std::endl;

  tgir::Scene const &scene = tgir::Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  tgir::Film film(W, H);
  std::mt19937_64 random;

  std::tcerr << _TEXT("> Create Global Photon Map") << std::endl;
  tgir::Scene::GetInstance().CreatePhotonMap(random);
  {
    tgir::PhotonMapProperty photonMapProperty =
        std::for_each(scene.GlobalPhotonMap().begin(),
                      scene.GlobalPhotonMap().end(), tgir::PhotonMapProperty());
    photonMapProperty.print();
  }

  tgir::PathVertex eyePathVertex;
  tgir::PhotonQuery photons(eyePathVertex);

  std::tcerr << _TEXT("> Do Direct Visualization") << std::endl;
  for (int y = 0; y < H; ++y) {
    std::tcerr << _TEXT("\r> ") << y;

    for (int x = 0; x < W; ++x) {
      scene.SampleLensPosition(&eyePathVertex);
      scene.Camera().GetPrimaryRayDirection((x + 0.5) / W, (y + 0.5) / H,
                                            &eyePathVertex.vOutgoingDirection);

      // 交差判定
      tgir::Intersection param;
      eyePathVertex.pGeometry = scene.FindIntersection(
          eyePathVertex.vPosition, eyePathVertex.vOutgoingDirection, &param);

      if (!eyePathVertex.pGeometry) {
        continue;  // through the air
      }

      // 交点とその点の基底ベクトルを計算
      {
        eyePathVertex.SetGeometricBasis(param);  // 基底
        eyePathVertex.back_side = param.is_back_side();
        eyePathVertex.vPosition +=
            eyePathVertex.vOutgoingDirection * param.t_max();  // 交点
      }

      scene.GlobalPhotonMap().kNN_query(eyePathVertex.vPosition, photons);

      // deposite
      film[y * W + x] =
          photons.GetEstimatedRadiance() * (CONFIG_ALBEDO * M_1_PI);
    }
  }

  // 定期的な保存
  film.SaveAsPfmRgb(_TEXT("0pm.pfm"), 1);
}
}  // end of namespace tgir

namespace tgir {
void Dvim::run() {
  std::tcerr << _TEXT("** Photon Mapping (1) **") << std::endl;

  tgir::Scene const &scene = tgir::Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  tgir::Film film(W, H);
  std::mt19937_64 random;

  std::tcerr << _TEXT("> Create Inportance Map") << std::endl;
  tgir::Scene::GetInstance().CreateViewImportanceMap(random);

  std::tcerr << _TEXT("> Create Global Photon Map") << std::endl;
  tgir::Scene::GetInstance().CreatePhotonMapWithImportonMap(random);
  {
    tgir::PhotonMapProperty photonMapProperty =
        std::for_each(scene.GlobalPhotonMap().begin(),
                      scene.GlobalPhotonMap().end(), tgir::PhotonMapProperty());
    photonMapProperty.print();
  }

  tgir::PathVertex eyePathVertex;
  tgir::PhotonQuery photons(eyePathVertex);

  std::tcerr << _TEXT("> Do Direct Visualization") << std::endl;
  for (int y = 0; y < H; ++y) {
    std::tcerr << _TEXT("\r> ") << y;

    for (int x = 0; x < W; ++x) {
      scene.SampleLensPosition(&eyePathVertex);
      scene.Camera().GetPrimaryRayDirection((x + 0.5) / W, (y + 0.5) / H,
                                            &eyePathVertex.vOutgoingDirection);

      // 交差判定
      tgir::Intersection param;
      eyePathVertex.pGeometry = scene.FindIntersection(
          eyePathVertex.vPosition, eyePathVertex.vOutgoingDirection, &param);

      if (!eyePathVertex.pGeometry) {
        continue;  // through the air
      }

      // 交点とその点の基底ベクトルを計算
      {
        eyePathVertex.SetGeometricBasis(param);  // 基底
        eyePathVertex.back_side = param.is_back_side();
        eyePathVertex.vPosition +=
            eyePathVertex.vOutgoingDirection * param.t_max();  // 交点
      }

      scene.GlobalPhotonMap().kNN_query(eyePathVertex.vPosition, photons);

      // deposite
      film[y * W + x] =
          photons.GetEstimatedRadiance() * (CONFIG_ALBEDO * M_1_PI);
    }
  }

  // 定期的な保存
  film.SaveAsPfmRgb(_TEXT("1pm.pfm"), 1);
}
}  // end of namespace tgir

namespace tgir {
void Dvpf::run() {
  std::tcerr << _TEXT("** Photon Mapping (2) **") << std::endl;

  tgir::Scene const &scene = tgir::Scene::GetInstance();

  int const W = scene.GetWidth();
  int const H = scene.GetHeight();

  tgir::Film film(W, H);
  std::mt19937_64 random;

  std::tcerr << _TEXT("> Create Inportance Map") << std::endl;
  tgir::Scene::GetInstance().CreateViewImportanceMap(random);

  std::tcerr << _TEXT("> Create Global Photon Map") << std::endl;
  tgir::Scene::GetInstance().CreatePhotonMapWithParticleFilter(random);
  {
    tgir::PhotonMapProperty photonMapProperty =
        std::for_each(scene.GlobalPhotonMap().begin(),
                      scene.GlobalPhotonMap().end(), tgir::PhotonMapProperty());
    photonMapProperty.print();
  }

  tgir::PathVertex eyePathVertex;
  tgir::PhotonQuery photons(eyePathVertex);

  std::tcerr << _TEXT("> Do Direct Visualization") << std::endl;
  for (int y = 0; y < H; ++y) {
    std::tcerr << _TEXT("\r> ") << y;

    for (int x = 0; x < W; ++x) {
      scene.SampleLensPosition(&eyePathVertex);
      scene.Camera().GetPrimaryRayDirection((x + 0.5) / W, (y + 0.5) / H,
                                            &eyePathVertex.vOutgoingDirection);

      // 交差判定
      tgir::Intersection param;
      eyePathVertex.pGeometry = scene.FindIntersection(
          eyePathVertex.vPosition, eyePathVertex.vOutgoingDirection, &param);

      if (!eyePathVertex.pGeometry) {
        continue;  // through the air
      }

      // 交点とその点の基底ベクトルを計算
      {
        eyePathVertex.SetGeometricBasis(param);  // 基底
        eyePathVertex.back_side = param.is_back_side();
        eyePathVertex.vPosition +=
            eyePathVertex.vOutgoingDirection * param.t_max();  // 交点
      }

      scene.GlobalPhotonMap().kNN_query(eyePathVertex.vPosition, photons);

      // deposite
      film[y * W + x] =
          photons.GetEstimatedRadiance() * (CONFIG_ALBEDO * M_1_PI);
    }
  }

  // 定期的な保存
  film.SaveAsPfmRgb(_TEXT("2pm.pfm"), 1);
}
}  // end of namespace tgir

#endif
