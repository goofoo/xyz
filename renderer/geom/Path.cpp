#include "Path.hpp"
#include "../core/Scene.hpp"

namespace tgir {
Path::Path(__in std::size_t const maxPathLength)
    : LightPathVertices(maxPathLength + 1 + 2),
      EyePathVertices(maxPathLength + 1 + 2),
      maxPathLength_(maxPathLength),
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

Path::Path(__in Path const &rhs)
    : vLensCoordinate(rhs.vLensCoordinate),
      LightPathVertices(rhs.LightPathVertices),
      EyePathVertices(rhs.EyePathVertices),
      maxPathLength_(rhs.maxPathLength_),
      contributions_(rhs.contributions_),
      contribution_from_particle_path_(rhs.contribution_from_particle_path_),
      contribution_from_indirect_path_(rhs.contribution_from_indirect_path_),
      contribution_from_explicit_path_(rhs.contribution_from_explicit_path_),
      contribution_from_implicit_path_(rhs.contribution_from_implicit_path_),
      contribution_from_specific_path_(rhs.contribution_from_specific_path_) {}

Path &Path::operator=(__in Path const &rhs) {
  if (this != &rhs) {
    vLensCoordinate = rhs.vLensCoordinate;
    LightPathVertices = rhs.LightPathVertices;
    EyePathVertices = rhs.EyePathVertices;

    maxPathLength_ = rhs.maxPathLength_;
    contributions_ = rhs.contributions_;
    contribution_from_particle_path_ = rhs.contribution_from_particle_path_;
    contribution_from_indirect_path_ = rhs.contribution_from_indirect_path_;
    contribution_from_explicit_path_ = rhs.contribution_from_explicit_path_;
    contribution_from_implicit_path_ = rhs.contribution_from_implicit_path_;
    contribution_from_specific_path_ = rhs.contribution_from_specific_path_;
  }

  return *this;
}

void Path::swap(__inout Path &rhs) {
  std::swap(vLensCoordinate, rhs.vLensCoordinate);
  std::swap(LightPathVertices, rhs.LightPathVertices);
  std::swap(EyePathVertices, rhs.EyePathVertices);

  std::swap(maxPathLength_, rhs.maxPathLength_);
  std::swap(contributions_, rhs.contributions_);
  std::swap(contribution_from_particle_path_,
            rhs.contribution_from_particle_path_);
  std::swap(contribution_from_indirect_path_,
            rhs.contribution_from_indirect_path_);
  std::swap(contribution_from_explicit_path_,
            rhs.contribution_from_explicit_path_);
  std::swap(contribution_from_implicit_path_,
            rhs.contribution_from_implicit_path_);
  std::swap(contribution_from_specific_path_,
            rhs.contribution_from_specific_path_);
}

void Path::Begin(__in EMutationType type) {
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
    for (std::size_t k = 1, N = maxPathLength_; k < N; ++k) {
      contributions_[k].first = INVALID_CONSTRIBUTION;
    }
    contribution_from_particle_path_ = 0;
    break;
  case MUTATION_EYE:
    contribution_from_explicit_path_ = 0;
    contribution_from_implicit_path_ = 0;
    contribution_from_specific_path_ = 0;
    break;
  }
}

void Path::End() {}

void Path::SetFilmPosition(__in std::size_t const &p) {
  contributions_[0].first = p;
}

void Path::AccumulateST(
    __in std::size_t const &uWavelengthIndex,  ///< [in] 波長の位置
    __in tgir::Spectrum const &cst)            ///< [in]
{
  tgir::SpectrumVector value;
  hi::spd2xyz(uWavelengthIndex, cst, value);
  contribution_from_indirect_path_ += value;
}

void Path::AccumulateS1(
    __in std::size_t const &p,                 ///< [in] 画素の位置
    __in std::size_t const &uWavelengthIndex,  ///< [in] 波長の位置
    __in std::size_t const &k,                 ///< [in] 経路長-1
    __in tgir::Spectrum const &cs1)            ///< [in]
{
  tgir::SpectrumVector value;
  hi::spd2xyz(uWavelengthIndex, cs1, value);
  contribution_from_particle_path_ += value;
  contributions_[k] = std::make_pair(p, value);  // NOTE: 経路長をいじってる
}

void Path::Accumulate1T(
    __in std::size_t const &uWavelengthIndex,  ///< [in] 波長の位置
    __in tgir::Spectrum const &c1t)            ///< [in]
{
  tgir::SpectrumVector value;
  hi::spd2xyz(uWavelengthIndex, c1t, value);
  contribution_from_explicit_path_ += value;
}

void Path::Accumulate0T(
    __in std::size_t const &uWavelengthIndex,  ///< [in] 波長の位置
    __in tgir::Spectrum const &c0t,            ///< [in]
    __in bool const &bSpecial)                 ///< [in]
{
  tgir::SpectrumVector value;
  hi::spd2xyz(uWavelengthIndex, c0t, value);

  if (bSpecial) {
    contribution_from_specific_path_ += value;
  } else {
    contribution_from_implicit_path_ += value;
  }
}

void Path::Values(
    __out std::vector<tgir::PixelDescriptor> *const pValues) const {
  for (std::size_t k = 0, size = maxPathLength_; k < size; ++k) {
    (*pValues)[k] = contributions_[k];
  }
  (*pValues)[0].second =
      contribution_from_indirect_path_ + contribution_from_explicit_path_ +
      contribution_from_implicit_path_ + contribution_from_specific_path_;
}

void Path::ValuesAndSpecial(
    __out std::vector<tgir::PixelDescriptor> *const pValues,
    __out tgir::PixelDescriptor *const pSpecial) const {
  for (std::size_t k = 0, size = maxPathLength_; k < size; ++k) {
    (*pValues)[k] = contributions_[k];
  }
  (*pValues)[0].second = contribution_from_indirect_path_ +
                         contribution_from_explicit_path_ +
                         contribution_from_implicit_path_;
  pSpecial->first = (*pValues)[0].first;
  pSpecial->second = contribution_from_specific_path_;
}

tgir::Real Path::PDFSimplifiedMLT() const {
#ifdef TGIR_CONFIG_FLG_USE_LUMINANCE_PDFS
  return contribution_from_particle_path_[1] +
         contribution_from_indirect_path_[1] +
         contribution_from_explicit_path_[1] +
         contribution_from_implicit_path_[1] +
         contribution_from_specific_path_[1];
#else TGIR_CONFIG_FLG_USE_LUMINANCE_PDFS
  return hi::sum(contribution_from_particle_path_) +
         hi::sum(contribution_from_indirect_path_) +
         hi::sum(contribution_from_explicit_path_) +
         hi::sum(contribution_from_implicit_path_) +
         hi::sum(contribution_from_specific_path_);
#endif TGIR_CONFIG_FLG_USE_LUMINANCE_PDFS
}

tgir::Real Path::PDFFullMLT() const {
  /*
  #ifdef TGIR_CONFIG_FLG_USE_LUMINANCE_PDFS
      return contribution_from_particle_path_[1]
           + contribution_from_indirect_path_[1]
           + contribution_from_implicit_path_[1]
           + contribution_from_specific_path_[1]
           ;
  #else  TGIR_CONFIG_FLG_USE_LUMINANCE_PDFS
      return hi::sum(contribution_from_particle_path_)
           + hi::sum(contribution_from_indirect_path_)
           + hi::sum(contribution_from_implicit_path_)
           + hi::sum(contribution_from_specific_path_)
           ;
  #endif TGIR_CONFIG_FLG_USE_LUMINANCE_PDFS
  */
  return PDFSimplifiedMLT();
}

void Path::PDFsFullRELT(__out tgir::Vector4 *const pPDFs) const {
#ifdef TGIR_CONFIG_FLG_USE_LUMINANCE_PDFS
  tgir::Real const a = contribution_from_particle_path_[1];
  tgir::Real const b = contribution_from_indirect_path_[1];
  tgir::Real const c = contribution_from_explicit_path_[1];
  tgir::Real const d = contribution_from_implicit_path_[1];
  tgir::Real const e = contribution_from_specific_path_[1];
#else TGIR_CONFIG_FLG_USE_LUMINANCE_PDFS
  tgir::Real const a = hi::sum(contribution_from_particle_path_);
  tgir::Real const b = hi::sum(contribution_from_indirect_path_);
  tgir::Real const c = hi::sum(contribution_from_explicit_path_);
  tgir::Real const d = hi::sum(contribution_from_implicit_path_);
  tgir::Real const e = hi::sum(contribution_from_specific_path_);
#endif TGIR_CONFIG_FLG_USE_LUMINANCE_PDFS
  /*
      static tgir::Real const m0a = 0.1, m0b = 0.3, m0c = 0.1, m0d = 0.4, m0e =
     0.8;
      static tgir::Real const m1a = 0.2, m1b = 0.3, m1c = 0.2, m1d = 0.3, m1e =
     0.1;
      static tgir::Real const m2a = 0.7, m2b = 0.4, m2c = 0.7, m2d = 0.3, m2e =
     0.1;
      (*pPDFs)[0] = m0a * a + m0b * b + m0c * c + m0d * d + m0e * e;
      (*pPDFs)[1] = m1a * a + m1b * b + m1c * c + m1d * d + m1e * e;
      (*pPDFs)[2] = m2a * a + m2b * b + m2c * c + m2d * d + m2e * e;
      (*pPDFs)[3] = 1;
  */
  (*pPDFs)[0] = +e + b * 1e-3;
  (*pPDFs)[1] = a + d + e + b;
  (*pPDFs)[2] = a + c + d + e + b;
  (*pPDFs)[3] = 1;
}
}  // end of namespace tgir
