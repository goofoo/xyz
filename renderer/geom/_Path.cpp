#include "Path.hpp"
#include "../core/Scene.hpp"

namespace tgir {
Path::Path(Path const &rhs)
    : maxPathLength_(rhs.maxPathLength_)
#if !defined(NDEBUG)
      ,
      lightPathVertices_(rhs.lightPathVertices_),
      eyePathVertices_(rhs.eyePathVertices_),
      pdf_(rhs.pdf_)
#else
      ,
      LightPathVertices(rhs.LightPathVertices),
      EyePathVertices(rhs.EyePathVertices)
#endif
      ,
      contributions_(rhs.contributions_) {
}

Path::Path(std::size_t const maxPathLength)
    : maxPathLength_(maxPathLength)
#if !defined(NDEBUG)
      ,
      lightPathVertices_(maxPathLength + 1),
      eyePathVertices_(maxPathLength + 1)
#else
      ,
      LightPathVertices(maxPathLength + 1 + 2),
      EyePathVertices(maxPathLength + 1 + 2)
#endif
      ,
      contributions_(maxPathLength) {
  /// サンプリング中で変化しない値を設定

  // 光源経路の仮想頂点
  LightPathVertices[0].back_side = false;
  LightPathVertices[0].bSpecular = false;
  LightPathVertices[0].rSamplingPrev = 0;
  LightPathVertices[0].rSamplingNext = 1;
  LightPathVertices[0].rGeometricFactor = 1;

  // 光源頂点
  LightPathVertices[1].back_side = false;
  LightPathVertices[1].bSpecular = false;
  LightPathVertices[1].rSamplingPrev = 0;
  LightPathVertices[1].rGeometricFactor = 1;

  // 視点経路の仮想頂点
  EyePathVertices[0].back_side = false;
  EyePathVertices[0].bSpecular = false;
  EyePathVertices[0].rSamplingPrev = 0;
  EyePathVertices[0].rSamplingNext = 1;
  EyePathVertices[0].rGeometricFactor = 1;

  // 視点頂点
  EyePathVertices[1].sQuantum = 1;
  EyePathVertices[1].back_side = false;
  EyePathVertices[1].bSpecular = false;
  EyePathVertices[1].rSamplingPrev = 0;
  EyePathVertices[1].rGeometricFactor = 0;
}

void Path::Begin(EMutationType type) {
  contribution_from_indirect_path_ = 0;

  switch (type) {
  case MUTATION_FULL:
    for (std::size_t k = 1, N = maxPathLength_; k < N; ++k) {
      contributions_[k].first = INVALID_CONSTRIBUTION;
    }
    contribution_from_particle_path_ = 0;
    contribution_from_explicit_path_ = 0;
    contribution_from_implicit_path_ = 0;
    contribution_from_specific_path_ = 0;
    break;
  case MUTATION_LIGHT:
    contribution_from_particle_path_ = 0;
    for (std::size_t k = 1, N = maxPathLength_; k < N; ++k) {
      contributions_[k].first = INVALID_CONSTRIBUTION;
    }
    break;
  case MUTATION_EYE:
    contribution_from_explicit_path_ = 0;
    contribution_from_implicit_path_ = 0;
    contribution_from_specific_path_ = 0;
    break;
  }
}

void Path::SetFilmPosition(std::size_t const p) { contributions_[0].first = p; }

void Path::End(std::vector<tgir::PixelDescriptor> &f) {
  contributions_[0].second =
      contribution_from_indirect_path_ + contribution_from_explicit_path_ +
      contribution_from_implicit_path_ + contribution_from_specific_path_;

  for (std::size_t k = 0, N = maxPathLength_; k < N; ++k) {
    f[k] = contributions_[k];
  }
}

void Path::AccumulateST(std::size_t const uWavelengthIndex,
                        tgir::Spectrum const cst) {
  tgir::SpectrumVector value;
  hi::spd2xyz(uWavelengthIndex, cst, value);
  contribution_from_indirect_path_ += value;
}

void Path::AccumulateS1(
    std::size_t const p,                 ///< [in] 画素の位置
    std::size_t const uWavelengthIndex,  ///< [in] 波長の位置
    std::size_t const k,                 ///< [in] 経路長-1
    tgir::Spectrum const cs1)            ///< [in]
{
  tgir::SpectrumVector value;
  hi::spd2xyz(uWavelengthIndex, cs1, value);
  contribution_from_particle_path_ += value;
  contributions_[k] = std::make_pair(p, value);
}

void Path::Accumulate1T(std::size_t const uWavelengthIndex,
                        tgir::Spectrum const c1t) {
  tgir::SpectrumVector value;
  hi::spd2xyz(uWavelengthIndex, c1t, value);
  contribution_from_explicit_path_ += value;
}

void Path::Accumulate0T(std::size_t const uWavelengthIndex,
                        tgir::Spectrum const c0t, bool const bSpecial) {
  tgir::SpectrumVector value;
  hi::spd2xyz(uWavelengthIndex, c0t, value);

  if (bSpecial) {
    contribution_from_specific_path_ += value;
  } else {
    contribution_from_implicit_path_ += value;
  }
}

}  // end of namespace tgir
