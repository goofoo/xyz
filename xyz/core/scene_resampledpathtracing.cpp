#include "scene.hpp"
#include "../bsdf/bsdf.hpp"
#include "../geom/triangle.hpp"
#include "../geom/intersection.hpp"
#include "../geom/pathvertex.hpp"
using namespace xyz;

void Scene::EvaluteResampledPathTracing(
    __in std::size_t const x, __in std::size_t const y,
    __out std::vector<CIE_XYZ_Color> *const pColor,
    __inout std::vector<PathVertex> &lightPathVertices,
    __inout std::vector<Particle> &ray, __inout std::vector<Particle> &tmp,
    __inout std::vector<float_t> &cdf, __inout IPrimarySample &sample) const {
  // 初期光線の設定
  {
    int const W = GetWidth();
    int const H = GetHeight();

    for (std::size_t i = 0, dy = 0; dy < CONFIG_BLOCK_SIZE; ++dy) {
      for (std::size_t dx = 0; dx < CONFIG_BLOCK_SIZE; ++dx, ++i) {
        // 光源位置のサンプリング
        PathVertex &theLightPathVertex = lightPathVertices[i];
        SampleLightPosition(sample.next(), sample.next(), &theLightPathVertex);
        theLightPathVertex.power = GetLightPowerXYZ();

        // 初期光線の設定
        PathVertex &theEyePathVertex = ray[i].eyePathVertex;
        SampleLensPosition(&theEyePathVertex);
        camera_.GetPrimaryRayDirection((x + dx + sample.next()) / W,
                                       (y + dy + sample.next()) / H,
                                       &theEyePathVertex.vOutgoingDirection);
        theEyePathVertex.bSpecular = true;
#ifndef CONFIG_RENDER_CAUSTICS
        theEyePathVertex.bDirectVisible = true;
#endif
        ray[i].nPixelIndex = i;

        (*pColor)[i] = CIE_XYZ_Color(CIE_XYZ_Color::value_type(0));
      }
    }
  }

  std::size_t n = CONFIG_BLOCK_SIZE * CONFIG_BLOCK_SIZE;  // 初期粒子数
  for (std::size_t k = 1; n > 0;)  // 再帰的に光線を追跡(末尾再帰)
  {
    std::size_t m = 0;  // 反射する光線数
    for (std::size_t i = 0; i < n; ++i) {
      PathVertex const &theLightPathVertex =
          lightPathVertices[ray[i].nPixelIndex];
      PathVertex &theEyePathVertex = ray[i].eyePathVertex;

      // 交差判定
      Intersection param;
      theEyePathVertex.pGeometry =
          FindIntersection(theEyePathVertex.vPosition,
                           theEyePathVertex.vOutgoingDirection, &param);

      if (!theEyePathVertex.pGeometry) {
        continue;
      }

      // 交点の座標とその点の基底ベクトルを計算
      theEyePathVertex.SetShadingBasis(param);  // 基底
      theEyePathVertex.bBackSide = param.is_back_side();

      // シェーディング法線の裏からあたった場合は，素直に終了
      if (-hi::dot(theEyePathVertex.vOutgoingDirection,
                   theEyePathVertex.vShadingNormal) <= 0) {
        continue;
      }

      theEyePathVertex.vPosition +=
          theEyePathVertex.vOutgoingDirection * param.t_max();  // 交点

      // 陰影処理
      Bsdf const &bsdf = *bsdfs_[theEyePathVertex.pGeometry->bsdf()];

// 直接光源が見えている場合
#ifdef CONFIG_DOUBLE_LIGHT
      if (bsdf.What() == Bsdf::LIGHT)
#else
      if ((bsdf.What() == Bsdf::LIGHT) && !theEyePathVertex.bBackSide)
#endif
      {
        if (theEyePathVertex.bSpecular) {
          // 鏡面反射から光源に到達する経路は，これ一つしか存在しないので重み付けなくてよい
          (*pColor)[ray[i].nPixelIndex] += theLightPathVertex.power *
                                           theEyePathVertex.power *
                                           (M_1_PI / GetLightArea());
        } else {
          float_t const fCosTheta =
              std::abs(hi::dot(theEyePathVertex.vOutgoingDirection,
                               theEyePathVertex.vShadingNormal));
          float_t const fDistanceSquared = hi::square_of(param.t_max());

          // 暗黙的サンプリングの経路密度
          float_t const fDensityOfImplicit =
              theEyePathVertex.fSamplingNext * fCosTheta / fDistanceSquared;

          // 明示的サンプリングの経路密度
          float_t const fDensityOfExplicit = hi::rcp(GetLightArea());

          // 経路密度の比
          float_t const fDensity = fDensityOfExplicit / fDensityOfImplicit;

          // 重みとして，光源方向のサンプリング確率を除し，光源上の点のサンプリング確率を掛ける
          (*pColor)[ray[i].nPixelIndex] +=
              theLightPathVertex.power * theEyePathVertex.power *
              (M_1_PI / GetLightArea()) *
              hi::rcp(1 + hi::square_of(fDensity));  // 重み
        }
        continue;
      }

      tmp[m++] = ray[i];
    }

    // NOTE: 反射回数の限界数を超えたら打ち切る
    if ((m <= 0) || (++k >= XYZ_CONFIG_kMaxRandomWalkDepth)) {
      break;
    }

    // 明示的サンプリング: Direct Lighting Calculation
    // 累積質量分布を計算
    cdf[0] = 0;
    for (std::size_t i = 0; i < m; ++i) {
      PathVertex const &theLightPathVertex =
          lightPathVertices[ray[i].nPixelIndex];
      PathVertex &theEyePathVertex = tmp[i].eyePathVertex;
      Bsdf const &bsdf = *bsdfs_[theEyePathVertex.pGeometry->bsdf()];
      (*pColor)[tmp[i].nPixelIndex] += bsdf.CalculateWeightedDirectLighting(
          &theEyePathVertex, theLightPathVertex, *this);

      // 累積質量分布
      float_t const t = std::min(
          float_t(1), theEyePathVertex.power[1] / CONFIG_TERMINATE_THRESHOLD);
      cdf[i + 1] = cdf[i] + t;
    }
    if (cdf[m] <= 0) {
      break;
    }

    // 粒子フィルタリング
    n = 0;  // 再帰的に追跡する光線数
    {
      std::size_t c = static_cast<std::size_t>(cdf[m]);
      if (sample.next() < (cdf[m] - c)) {
        ++c;
      }

      // 再サンプリング
      typedef std::vector<float_t>::const_iterator const_iterator_t;
      const_iterator_t const begin = cdf.begin();
      const_iterator_t const end = begin + m + 1;
      for (std::size_t i = 0; i < c; ++i) {  // c個のサンプルを生成する
        float_t const value = cdf[m] * (i + sample.next()) / c;
        const_iterator_t const it = std::upper_bound(begin, end, value);

        std::size_t const j = it - begin - 1;
        ray[n] = tmp[j];

        PathVertex &theEyePathVertex = ray[n].eyePathVertex;
        theEyePathVertex.power *=
            cdf[m] / (c * (cdf[j + 1] - cdf[j]));  // 重みを調整
        Bsdf const &bsdf = *bsdfs_[theEyePathVertex.pGeometry->bsdf()];
        if (bsdf.NextScatteringDirection(&theEyePathVertex, sample)) {
          ++n;
        }
      }
    }
  }
}
