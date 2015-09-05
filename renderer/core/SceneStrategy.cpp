#include "Scene.hpp"
#include "bsdf/Bsdf.hpp"
#include "geom/Triangle.hpp"
#include "geom/Intersection.hpp"
#include "geom/PathVertex.hpp"

namespace tgir {
/// <summary>
/// path tracing with multiple importance sampling
/// </summary>
void Scene::EvaluteRussianRoulette(std::size_t const x, std::size_t const y,
                                   tgir::SpectrumVector &sColor,
                                   std::mt19937_64 &random) const {
  tgir::PathVertex eyePathVertex;
  tgir::PathVertex lightPathVertex;

  std::pair<tgir::Real, std::size_t> wavelength;
  wavelength.first = std::uniform_real_distribution<tgir::Real>()(random);
  wavelength.second =
      static_cast<std::size_t>(wavelength.first * tgir::SpectrumData::Size);

  SampleLightPosition(std::uniform_real_distribution<tgir::Real>()(random),
                      std::uniform_real_distribution<tgir::Real>()(random),
                      &lightPathVertex);
  lightPathVertex.sQuantum = tgir::GetLightPower()[wavelength.second];

#ifdef CONFIG_PINHOLE_CAMERA_MODEL
  SampleLensPosition(eyePathVertex);
#else
  tgir::Vector2 vLensCoordinate;
  SampleLensPosition(std::uniform_real_distribution<tgir::Real>()(random),
                     std::uniform_real_distribution<tgir::Real>()(random),
                     &vLensCoordinate, &eyePathVertex);
#endif CONFIG_PINHOLE_CAMERA_MODEL

  // 初期光線の設定
  {
#ifdef CONFIG_PINHOLE_CAMERA_MODEL
    camera_.GetPrimaryRayDirection(
        (x + std::uniform_real_distribution<tgir::Real>()(random)) / GetWidth(),
        (y + std::uniform_real_distribution<tgir::Real>()(random)) /
            GetHeight(),
        eyePathVertex.vOutgoingDirection);
#else
    camera_.GetPrimaryRayDirection(
        (x + std::uniform_real_distribution<tgir::Real>()(random)) / GetWidth(),
        (y + std::uniform_real_distribution<tgir::Real>()(random)) /
            GetHeight(),
        vLensCoordinate, &eyePathVertex.vOutgoingDirection);
#endif CONFIG_PINHOLE_CAMERA_MODEL
    eyePathVertex.bSpecular = true;
  }

  tgir::Spectrum sRadiance = 0;

#ifndef CONFIG_RENDER_CAUSTICS
  bool bDirectVisible = true;
#endif

  for (std::size_t k = 1;;)  // 再帰的に光線を追跡(末尾再帰)
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
      eyePathVertex.SetShadingBasis(param);  // 基底
      eyePathVertex.back_side = param.is_back_side();

      if (-hi::dot(eyePathVertex.vOutgoingDirection,
                   eyePathVertex.vShadingNormal) <= 0) {
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
    tgir::Bsdf const &bsdf = *bsdfs_[eyePathVertex.pGeometry->bsdf()];

// 直接光源が見えている場合
#ifdef CONFIG_DOUBLE_LIGHT
    if (bsdf.What() == tgir::Bsdf::LIGHT)
#else
    if ((bsdf.What() == tgir::Bsdf::LIGHT) && !eyePathVertex.back_side)
#endif
    {
      if (eyePathVertex.bSpecular) {
        // 鏡面反射から光源に到達する経路は，これ一つしか存在しないので重み付けなくてよい
        sRadiance += lightPathVertex.sQuantum * eyePathVertex.sQuantum *
                     (M_1_PI / GetLightArea());
      } else {
        tgir::Real const rCosTheta = std::abs(hi::dot(
            eyePathVertex.vOutgoingDirection, eyePathVertex.vShadingNormal));
        tgir::Real const rDistanceSquared = hi::square_of(param.t_max());
        tgir::Real const rDensity =
            rDistanceSquared *
            hi::rcp(eyePathVertex.rSamplingNext * rCosTheta * GetLightArea());

        // 重みとして，光源方向のサンプリング確率を除し，光源上の点のサンプリング確率を掛ける
        sRadiance += lightPathVertex.sQuantum * eyePathVertex.sQuantum *
                     (M_1_PI / GetLightArea()) *
                     hi::rcp(1 + hi::square_of(rDensity));  // 重み
      }
      break;
    }

    // Direct Lighting Calculation
    sRadiance += bsdf.CalculateWeightedDirectLighting(
        eyePathVertex, lightPathVertex, wavelength.second, *this);

    // NOTE: 反射回数の限界数を超えたら打ち切る
    if (++k >= TGIR_CONFIG_kMaxRandomWalkDepth) {
      break;
    }

    // Russian-Roulette
    tgir::Real const t = eyePathVertex.sQuantum;
    if (t < CONFIG_TERMINATE_THRESHOLD) {
      tgir::Real const q = t / CONFIG_TERMINATE_THRESHOLD;
      if (std::uniform_real_distribution<tgir::Real>()(random) >=
          q)  // Random Termination
      {
        break;  // Terminate the Path
      }
      eyePathVertex.sQuantum *= hi::rcp(q);
    }

    // Indirect Lighting Calculation
    if (!bsdf.GetScatterVector(eyePathVertex, wavelength.first, random)) {
      break;
    }
  }

  hi::spd2xyz(wavelength.second, sRadiance, sColor);
}
}  // end of namespace tgir

/*
namespace tgir
{
  void Scene::EvaluteGoWithTheWinner(
    std::size_t const x, std::size_t const y,
    tgir::SpectrumVector & sColor,
    std::mt19937_64 & random) const
  {
    tgir::PathVertex eyePathVertex;
    tgir::PathVertex lightPathVertex;

    std::pair<tgir::Real, std::size_t> wavelength;
    wavelength.first = std::uniform_real_distribution<tgir::Real>()(random);
    wavelength.second = static_cast<std::size_t>(wavelength.first *
tgir::SpectrumData::Size);

    SampleLightPosition(tgir::Vector2(
      std::uniform_real_distribution<tgir::Real>()(random),
      std::uniform_real_distribution<tgir::Real>()(random)), lightPathVertex);
    lightPathVertex.sQuantum = tgir::GetLightPower()[wavelength.second];
    SampleLensPosition(eyePathVertex);

    // 初期光線の設定
    {
      camera_.GetPrimaryRayDirection(
        (x + std::uniform_real_distribution<tgir::Real>()(random)) / GetWidth(),
        (y + std::uniform_real_distribution<tgir::Real>()(random)) /
GetHeight(),
        eyePathVertex.vOutgoingDirection);
      eyePathVertex.bSpecular = true;
    }

    hi::spd2xyz(wavelength.second,
      TraceGoWithTheWinner(1, true, wavelength, lightPathVertex, eyePathVertex,
random), sColor);
  }

  tgir::Spectrum Scene::TraceGoWithTheWinner(
    std::size_t depth, bool bDirectVisible,
    std::pair<tgir::Real, std::size_t> const & wavelength,
    tgir::PathVertex const & lightPathVertex,
    tgir::PathVertex eyePathVertex,
    std::mt19937_64 & random) const
  {
    // 交差判定
    tgir::Intersection param;
    eyePathVertex.pGeometry = FindIntersection(
      eyePathVertex.vPosition,
      eyePathVertex.vOutgoingDirection, param);
    if (!eyePathVertex.pGeometry)
    {
      return 0;
    }

    // 交点とその点の基底ベクトルを計算
    {
      eyePathVertex.SetShadingBasis(param); // 基底
      eyePathVertex.back_side = param.is_back_side();

      if (-hi::dot(eyePathVertex.vOutgoingDirection,
eyePathVertex.vShadingNormal) <= 0)
      {
        return 0;
      }

      eyePathVertex.vPosition += eyePathVertex.vOutgoingDirection *
param.t_max(); // 交点
    }

#ifndef CONFIG_RENDER_CAUSTICS
    if (!eyePathVertex.bSpecular)
    {
      bDirectVisible = false;
    }
#endif

    ++depth;

    // 陰影処理
    tgir::Bsdf const & bsdf = *bsdfs_[eyePathVertex.pGeometry->bsdf()];

    // Direct Lighting Calculation
    tgir::Real sRadiance = bsdf.CalculateDirectLighting(eyePathVertex,
lightPathVertex, wavelength.second, *this);

    if (bsdf.What() == tgir::Bsdf::LIGHT)
    {
#ifndef CONFIG_RENDER_CAUSTICS
      return ((eyePathVertex.bSpecular && bDirectVisible) ?
#else
      return ((eyePathVertex.bSpecular) ?
#endif
        lightPathVertex.sQuantum * eyePathVertex.sQuantum *
(M_1_PI/lightsCdf_.back()) : 0) + sRadiance;
    }

    // NOTE: 反射回数の限界数を超えたら打ち切る
    if (depth >= TGIR_CONFIG_kMaxRandomWalkDepth)
    {
      return sRadiance;
    }

    tgir::Real const n = eyePathVertex.sQuantum * bsdf.GetDensityVariance();
    if (n <= 0)
    {
      return sRadiance;
    }
    eyePathVertex.sQuantum /= n;

    std::size_t const m = static_cast<std::size_t>(n);
    for (std::size_t i = 0; i < m; ++i)
    {
      if (bsdf.GetScatterVector(eyePathVertex, wavelength.first, random))
      {
        sRadiance += TraceGoWithTheWinner(depth, bDirectVisible, wavelength,
lightPathVertex, eyePathVertex, random);
      }
    }
    if (std::uniform_real_distribution<tgir::Real>()(random) < (n-m))
    {
      if (bsdf.GetScatterVector(eyePathVertex, wavelength.first, random))
      {
        sRadiance += TraceGoWithTheWinner(depth, bDirectVisible, wavelength,
lightPathVertex, eyePathVertex, random);
      }
    }
    return sRadiance;
  }
} // end of namespace tgir
*/

namespace tgir {
void Scene::EvaluteParticleFilter(std::size_t const x, std::size_t const y,
                                  std::vector<tgir::SpectrumVector> &rad,
                                  std::vector<tgir::ParticleInfo> &info,
                                  std::vector<tgir::Particle> &ray,
                                  std::vector<tgir::Particle> &tmp,
                                  std::vector<tgir::Real> &cdf,
                                  std::mt19937_64 &random) const {
  {
    std::size_t i = 0;
    for (std::size_t dy = 0; dy < CONFIG_BLOCK_SIZE; ++dy)
      for (std::size_t dx = 0; dx < CONFIG_BLOCK_SIZE; ++dx)
        for (std::size_t da = 0; da < CONFIG_LIGHT_SIZE; ++da)
          for (std::size_t db = 0; db < CONFIG_LIGHT_SIZE; ++db)
            for (std::size_t ds = 0; ds < CONFIG_COLOR_SIZE; ++ds, ++i) {
              info[i].sRadiance = 0;

              // 波長の選択
              info[i].wavelength.first =
                  (ds + std::uniform_real_distribution<tgir::Real>()(random)) /
                  CONFIG_COLOR_SIZE;
              info[i].wavelength.second = static_cast<std::size_t>(
                  info[i].wavelength.first * tgir::SpectrumData::Size);

              // 光源の選択
              SampleLightPosition(
                  (da + std::uniform_real_distribution<tgir::Real>()(random)) /
                      CONFIG_LIGHT_SIZE,
                  (db + std::uniform_real_distribution<tgir::Real>()(random)) /
                      CONFIG_LIGHT_SIZE,
                  &info[i].lightPathVertex);
              info[i].lightPathVertex.sQuantum =
                  tgir::GetLightPower()[info[i].wavelength.second];

              // 初期光線の設定
              ray[i].uInfoIndex = i;
              SampleLensPosition(&ray[i].eyePathVertex);
              camera_.GetPrimaryRayDirection(
                  ((x + dx) +
                   std::uniform_real_distribution<tgir::Real>()(random)) /
                      GetWidth(),
                  ((y + dy) +
                   std::uniform_real_distribution<tgir::Real>()(random)) /
                      GetHeight(),
                  &ray[i].eyePathVertex.vOutgoingDirection);
              ray[i].eyePathVertex.bSpecular = true;
#ifndef CONFIG_RENDER_CAUSTICS
              ray[i].bDirectVisible = true;
#endif
            }
  }

  std::size_t n = CONFIG_FILTER_SIZE;
  for (std::size_t k = 1; n > 0;) {
    // 交点を計算して，反射ポテンシャルを求める．
    std::size_t nPathVertex = 0;  // 反射する頂点の個数

    for (std::size_t i = 0; i < n; ++i) {
      tgir::Particle &p = ray[i];
      tgir::Intersection param;
      p.eyePathVertex.pGeometry =
          FindIntersection(p.eyePathVertex.vPosition,
                           p.eyePathVertex.vOutgoingDirection, &param);
      if (!p.eyePathVertex.pGeometry) {
        continue;
      }

      // 交点の座標とその点の基底ベクトルを計算
      p.eyePathVertex.SetShadingBasis(param);  // 基底
      p.eyePathVertex.back_side = param.is_back_side();
      if (-hi::dot(p.eyePathVertex.vOutgoingDirection,
                   p.eyePathVertex.vShadingNormal) <= 0) {
        continue;
      }
      p.eyePathVertex.vPosition +=
          p.eyePathVertex.vOutgoingDirection * param.t_max();

#ifndef CONFIG_RENDER_CAUSTICS
      if (!p.eyePathVertex.bSpecular) {
        p.bDirectVisible = false;
      }
#endif

      tgir::Bsdf const &bsdf = *bsdfs_[p.eyePathVertex.pGeometry->bsdf()];
      tgir::PathVertex const &lightPathVertex =
          info[p.uInfoIndex].lightPathVertex;

// 直接光源が見えている場合
#ifdef CONFIG_DOUBLE_LIGHT
      if (bsdf.What() == tgir::Bsdf::LIGHT)
#else
      if ((bsdf.What() == tgir::Bsdf::LIGHT) && !p.eyePathVertex.back_side)
#endif
      {
        if (p.eyePathVertex.bSpecular) {
          // 鏡面反射から光源に到達する経路は，これ一つしか存在しないので重み付けなくてよい
          info[p.uInfoIndex].sRadiance += lightPathVertex.sQuantum *
                                          p.eyePathVertex.sQuantum *
                                          (M_1_PI / GetLightArea());
        } else {
          /*
            光源上の点を選択した経路を生成する確率: p(x) = 1/A
            偶然に光源に到達する経路を生成する確率: p(y) = p(w) cos\theta / d^2
            p^2(x) / (p^2(x) + p^2(y)) = 1 / (1 + (p(y)/p(x))^2)
          */
          tgir::Real const rCosTheta =
              std::abs(hi::dot(p.eyePathVertex.vOutgoingDirection,
                               p.eyePathVertex.vShadingNormal));
          tgir::Real const rDistanceSquared = hi::square_of(param.t_max());
          tgir::Real const rDensity =
              rDistanceSquared * hi::rcp(p.eyePathVertex.rSamplingNext *
                                         rCosTheta * GetLightArea());

          // 重みとして，光源方向のサンプリング確率を除し，光源上の点のサンプリング確率を掛ける
          info[p.uInfoIndex].sRadiance +=
              lightPathVertex.sQuantum * p.eyePathVertex.sQuantum *
              (M_1_PI / GetLightArea()) *
              hi::rcp(1 + hi::square_of(rDensity));  // 重み
        }
        continue;
      }

      // 直接照明計算
      info[p.uInfoIndex].sRadiance += bsdf.CalculateWeightedDirectLighting(
          p.eyePathVertex, lightPathVertex,
          info[p.uInfoIndex].wavelength.second, *this);
      tmp[nPathVertex++] = p;
    }

    if ((nPathVertex <= 0) || (++k >= TGIR_CONFIG_kMaxRandomWalkDepth)) {
      break;
    }

    // 累積質量分布を計算
    cdf[0] = 0;
    for (std::size_t i = 0; i < nPathVertex; ++i) {
#if ((CONFIG_SS_TYPE == 1) || (CONFIG_SS_TYPE == 2))
      tgir::Real const t =
          std::min(tmp[i].eyePathVertex.sQuantum / CONFIG_TERMINATE_THRESHOLD,
                   tgir::Real(1));
#else
      tgir::Real const t = tmp[i].eyePathVertex.sQuantum;
#endif
      cdf[i + 1] = cdf[i] + t;
    }
    if (cdf[nPathVertex] <= 0) {
      break;
    }
#if (CONFIG_SS_TYPE == 1)
    tgir::Real const rRayCount = cdf[nPathVertex];
#else
    tgir::Real ess = 0;  // effective sample size
    {
      tgir::Real const nf = hi::rcp(cdf[nPathVertex]);
      for (std::size_t i = 0; i < nPathVertex; ++i) {
        ess += hi::square_of((cdf[i + 1] - cdf[i]) * nf);
      }
    }
    tgir::Real const rRayCount = hi::rcp(ess);
#endif

    // 粒子フィルタリング
    n = 0;
    {
      std::size_t nCount = static_cast<std::size_t>(rRayCount);
      if (std::uniform_real_distribution<tgir::Real>()(random) <
          (rRayCount - nCount)) {
        ++nCount;
      }

      // 再サンプリング
      std::vector<tgir::Real>::const_iterator const begin = cdf.begin();
      std::vector<tgir::Real>::const_iterator const end =
          cdf.begin() + nPathVertex + 1;
      for (std::size_t i = 0; i < nCount; ++i) {
        tgir::Real const value =
            cdf[nPathVertex] *
            (i + std::uniform_real_distribution<tgir::Real>()(random)) / nCount;
        std::vector<tgir::Real>::const_iterator const it =
            std::upper_bound(begin, end, value);

        std::size_t const uIndex = it - begin - 1;
        tgir::Particle &particle = ray[n];
        particle = tmp[uIndex];
        particle.eyePathVertex.sQuantum *=
            cdf[nPathVertex] / (rRayCount * (cdf[uIndex + 1] - cdf[uIndex]));

        if (bsdfs_[particle.eyePathVertex.pGeometry->bsdf()]->GetScatterVector(
                particle.eyePathVertex,
                info[particle.uInfoIndex].wavelength.first, random)) {
          ++n;
        }
      }
    }
  }

  for (std::size_t i = 0; i < CONFIG_FILTER_SIZE; ++i) {
    hi::spd2xyz(info[i].wavelength.second, info[i].sRadiance, rad[i]);
  }
}
}  // end of namespace tgir

namespace tgir {
void Scene::EvaluteParticleFilter(std::size_t const x, std::size_t const y,
                                  tgir::SpectrumVector &sColor,
                                  std::vector<tgir::ParticleInfo> &info,
                                  std::vector<tgir::Particle> &ray,
                                  std::vector<tgir::Particle> &tmp,
                                  std::vector<tgir::Real> &cdf,
                                  std::mt19937_64 &random) const {
  {
    std::size_t i = 0;
    for (std::size_t da = 0; da < CONFIG_LIGHT_SIZE; ++da)
      for (std::size_t db = 0; db < CONFIG_LIGHT_SIZE; ++db)
        for (std::size_t dy = 0; dy < CONFIG_BLOCK_SIZE; ++dy)
          for (std::size_t dx = 0; dx < CONFIG_BLOCK_SIZE; ++dx, ++i) {
            info[i].sRadiance = 0;

            // 光源の選択
            SampleLightPosition(
                (da + std::uniform_real_distribution<tgir::Real>()(random)) /
                    CONFIG_LIGHT_SIZE,
                (db + std::uniform_real_distribution<tgir::Real>()(random)) /
                    CONFIG_LIGHT_SIZE,
                &info[i].lightPathVertex);
            info[i].lightPathVertex.sQuantum =
                tgir::GetLightPower()[info[i].wavelength.second];

            // 初期光線の設定
            ray[i].uInfoIndex = i;
            SampleLensPosition(&ray[i].eyePathVertex);
            camera_.GetPrimaryRayDirection(
                (x + std::uniform_real_distribution<tgir::Real>()(random)) /
                    GetWidth(),
                (y + std::uniform_real_distribution<tgir::Real>()(random)) /
                    GetHeight(),
                &ray[i].eyePathVertex.vOutgoingDirection);
            ray[i].eyePathVertex.bSpecular = true;
#ifndef CONFIG_RENDER_CAUSTICS
            ray[i].bDirectVisible = true;
#endif
          }
  }

  std::size_t n = CONFIG_FILTER_SIZE;
  for (std::size_t k = 1; n > 0;) {
    // 交点を計算して，反射ポテンシャルを求める．
    std::size_t nPathVertex = 0;  // 反射する頂点の個数

    for (std::size_t i = 0; i < n; ++i) {
      tgir::Particle &p = ray[i];
      tgir::Intersection param;
      p.eyePathVertex.pGeometry =
          FindIntersection(p.eyePathVertex.vPosition,
                           p.eyePathVertex.vOutgoingDirection, &param);
      if (!p.eyePathVertex.pGeometry) {
        continue;
      }

      // 交点の座標とその点の基底ベクトルを計算
      p.eyePathVertex.SetShadingBasis(param);  // 基底
      p.eyePathVertex.back_side = param.is_back_side();
      if (-hi::dot(p.eyePathVertex.vOutgoingDirection,
                   p.eyePathVertex.vShadingNormal) <= 0) {
        continue;
      }
      p.eyePathVertex.vPosition +=
          p.eyePathVertex.vOutgoingDirection * param.t_max();

#ifndef CONFIG_RENDER_CAUSTICS
      if (!p.eyePathVertex.bSpecular) {
        p.bDirectVisible = false;
      }
#endif

      tgir::Bsdf const &bsdf = *bsdfs_[p.eyePathVertex.pGeometry->bsdf()];
      tgir::PathVertex const &lightPathVertex =
          info[p.uInfoIndex].lightPathVertex;

// 直接光源が見えている場合
#ifdef CONFIG_DOUBLE_LIGHT
      if (bsdf.What() == tgir::Bsdf::LIGHT)
#else
      if ((bsdf.What() == tgir::Bsdf::LIGHT) && !p.eyePathVertex.back_side)
#endif
      {
        if (p.eyePathVertex.bSpecular) {
          // 鏡面反射から光源に到達する経路は，これ一つしか存在しないので重み付けなくてよい
          info[p.uInfoIndex].sRadiance += lightPathVertex.sQuantum *
                                          p.eyePathVertex.sQuantum *
                                          (M_1_PI / GetLightArea());
        } else {
          tgir::Real const rCosTheta =
              std::abs(hi::dot(p.eyePathVertex.vOutgoingDirection,
                               p.eyePathVertex.vShadingNormal));
          tgir::Real const rDistanceSquared = hi::square_of(param.t_max());
          tgir::Real const rDensity =
              rDistanceSquared * hi::rcp(p.eyePathVertex.rSamplingNext *
                                         rCosTheta * GetLightArea());

          // 重みとして，光源方向のサンプリング確率を除し，光源上の点のサンプリング確率を掛ける
          info[p.uInfoIndex].sRadiance +=
              lightPathVertex.sQuantum * p.eyePathVertex.sQuantum *
              (M_1_PI / GetLightArea()) *
              hi::rcp(1 + hi::square_of(rDensity));  // 重み
        }
        continue;
      }

      // 直接照明計算
      info[p.uInfoIndex].sRadiance += bsdf.CalculateWeightedDirectLighting(
          p.eyePathVertex, lightPathVertex,
          info[p.uInfoIndex].wavelength.second, *this);
      tmp[nPathVertex++] = p;
    }

    if ((nPathVertex <= 0) || (++k >= TGIR_CONFIG_kMaxRandomWalkDepth)) {
      break;
    }

    // 累積質量分布を計算
    cdf[0] = 0;
    for (std::size_t i = 0; i < nPathVertex; ++i) {
#if ((CONFIG_SS_TYPE == 1) || (CONFIG_SS_TYPE == 2))
      tgir::Real const t =
          std::min(tmp[i].eyePathVertex.sQuantum / CONFIG_TERMINATE_THRESHOLD,
                   tgir::Real(1));
#else
      tgir::Real const t = tmp[i].eyePathVertex.sQuantum;
#endif
      cdf[i + 1] = cdf[i] + t;
    }
    if (cdf[nPathVertex] <= 0) {
      break;
    }
#if (CONFIG_SS_TYPE == 1)
    tgir::Real const rRayCount = cdf[nPathVertex];
#else
    tgir::Real ess = 0;  // effective sample size
    {
      tgir::Real const nf = hi::rcp(cdf[nPathVertex]);
      for (std::size_t i = 0; i < nPathVertex; ++i) {
        ess += hi::square_of((cdf[i + 1] - cdf[i]) * nf);
      }
    }
    tgir::Real const rRayCount = hi::rcp(ess);
#endif

    // 粒子フィルタリング
    n = 0;
    {
      std::size_t nCount = static_cast<std::size_t>(rRayCount);
      if (std::uniform_real_distribution<tgir::Real>()(random) <
          (rRayCount - nCount)) {
        ++nCount;
      }

      // 再サンプリング
      std::vector<tgir::Real>::const_iterator const begin = cdf.begin();
      std::vector<tgir::Real>::const_iterator const end =
          begin + nPathVertex + 1;
      for (std::size_t i = 0; i < nCount; ++i) {
        tgir::Real const value =
            cdf[nPathVertex] *
            (i + std::uniform_real_distribution<tgir::Real>()(random)) / nCount;
        std::vector<tgir::Real>::const_iterator const it =
            std::upper_bound(begin, end, value);

        std::size_t const uIndex = it - begin - 1;
        tgir::Particle &particle = ray[n];
        particle = tmp[uIndex];
        particle.eyePathVertex.sQuantum *=
            cdf[nPathVertex] / (rRayCount * (cdf[uIndex + 1] - cdf[uIndex]));

        if (bsdfs_[particle.eyePathVertex.pGeometry->bsdf()]->GetScatterVector(
                particle.eyePathVertex,
                info[particle.uInfoIndex].wavelength.first, random)) {
          ++n;
        }
      }
    }
  }

  sColor = 0;
  for (std::size_t i = 0; i < CONFIG_FILTER_SIZE; ++i) {
    tgir::SpectrumVector sRadiance;
    hi::spd2xyz(info[i].wavelength.second, info[i].sRadiance, sRadiance);
    sColor += sRadiance;
  }
}
}  // end of namespace tgir

/*
namespace tgir
{
  /// <summary>
  /// Importance Driven Path Tracing
  /// </summary>
  void Scene::EvaluteImportanceDriven(
    std::size_t const x,
    std::size_t const y,
    tgir::SpectrumVector & sColor,
    tgir::PathVertex & eyePathVertex,
    tgir::ImportanceQuery & query,
    std::mt19937_64 & random) const
  {
    tgir::PathVertex lightPathVertex;

    std::pair<tgir::Real, std::size_t> wavelength;
    wavelength.first = std::uniform_real_distribution<tgir::Real>()(random);
    wavelength.second = static_cast<std::size_t>(wavelength.first *
tgir::SpectrumData::Size);

    SampleLightPosition(tgir::Vector2(
      std::uniform_real_distribution<tgir::Real>()(random),
      std::uniform_real_distribution<tgir::Real>()(random)), lightPathVertex);
    lightPathVertex.sQuantum = tgir::GetLightPower()[wavelength.second];

    // 初期光線の設定
    {
      SampleLensPosition(eyePathVertex);
      camera_.GetPrimaryRayDirection(
        (x + std::uniform_real_distribution<tgir::Real>()(random)) / GetWidth(),
        (y + std::uniform_real_distribution<tgir::Real>()(random)) /
GetHeight(),
        eyePathVertex.vOutgoingDirection);
      eyePathVertex.bSpecular = true;
    }

    tgir::Spectrum sRadiance = 0;

#ifndef CONFIG_RENDER_CAUSTICS
    bool bDirectVisible = true;
#endif

    for (std::size_t k = 1; ;) // 再帰的に光線を追跡(末尾再帰)
    {
      // 交差判定
      tgir::Intersection param;
      eyePathVertex.pGeometry = FindIntersection(
        eyePathVertex.vPosition,
        eyePathVertex.vOutgoingDirection, param);

      if (!eyePathVertex.pGeometry)
      {
        break; // through the air
      }

      // 交点とその点の基底ベクトルを計算
      {
        eyePathVertex.SetShadingBasis(param); // 基底
        eyePathVertex.back_side = param.is_back_side();

        if (-hi::dot(eyePathVertex.vOutgoingDirection,
eyePathVertex.vShadingNormal) <= 0)
        {
          break;
        }

        eyePathVertex.vPosition += eyePathVertex.vOutgoingDirection *
param.t_max(); // 交点
      }

#ifndef CONFIG_RENDER_CAUSTICS
      if (!eyePathVertex.bSpecular)
      {
        bDirectVisible = false;
      }
#endif

      // 陰影処理
      tgir::Bsdf const & bsdf = *bsdfs_[eyePathVertex.pGeometry->bsdf()];

      if (bsdf.What() == tgir::Bsdf::LIGHT)
      {
#ifdef CONFIG_RENDER_CAUSTICS
        if (eyePathVertex.bSpecular)
#else
        if (eyePathVertex.bSpecular && bDirectVisible)
#endif
        {
          sRadiance +=
            lightPathVertex.sQuantum
            * eyePathVertex.sQuantum
            * (M_1_PI/lightsCdf_.back());
        }
        break;
      }

      // Direct Lighting Calculation
      sRadiance += bsdf.CalculateDirectLighting(eyePathVertex, lightPathVertex,
wavelength.second, *this);

      // NOTE: 反射回数の限界数を超えたら打ち切る
      if (++k >= TGIR_CONFIG_kMaxRandomWalkDepth)
      {
        break;
      }

      // Russian-Roulette
      tgir::Real const t = eyePathVertex.sQuantum;
      if (t < CONFIG_TERMINATE_THRESHOLD)
      {
        tgir::Real const q = t / CONFIG_TERMINATE_THRESHOLD;
        if (std::uniform_real_distribution<tgir::Real>()(random) >= q) // Random
Termination
        {
          break; // Terminate the Path
        }
        eyePathVertex.sQuantum *= hi::rcp(q);
      }

      // Indirect Lighting Calculation
      if (!bsdf.GetScatterVector(eyePathVertex, wavelength.first,
importanceMap_, query, random))
      {
        break;
      }
    }

    hi::spd2xyz(wavelength.second, sRadiance, sColor);
  }

} // end of namespace tgir
*/