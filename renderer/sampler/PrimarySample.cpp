#include "PrimarySample.hpp"

namespace {
inline tgir::Real Mutate(tgir::Real const fPrevMu, tgir::Real const r2,
                         tgir::Real const r3, std::mt19937_64 &random) {
  tgir::Real const fDv =
      r2 * std::exp(r3 * std::uniform_real_distribution<tgir::Real>()(random));
  tgir::Real const fMu =
      fPrevMu + ((std::uniform_real_distribution<tgir::Real>()(random) < 0.5)
                     ? -fDv
                     : fDv);
  return fMu - std::floor(fMu);
}
}

namespace tgir {
void PrimarySample::Init(std::mt19937_64 &random) {
  {
    muWavelength_.first = std::uniform_real_distribution<tgir::Real>()(random);
    muWavelength_.second = static_cast<std::size_t>(muWavelength_.first *
                                                    tgir::SpectrumData::Size);
  }
  for (std::size_t i = 0, size = muLight_.size(); i < size; ++i) {
    muLight_[i][0] = std::uniform_real_distribution<tgir::Real>()(random);
    muLight_[i][1] = std::uniform_real_distribution<tgir::Real>()(random);
    muLight_[i][2] = std::uniform_real_distribution<tgir::Real>()(random);
  }
  for (std::size_t i = 0, size = muEye_.size(); i < size; ++i) {
    muEye_[i][0] = std::uniform_real_distribution<tgir::Real>()(random);
    muEye_[i][1] = std::uniform_real_distribution<tgir::Real>()(random);
    muEye_[i][2] = std::uniform_real_distribution<tgir::Real>()(random);
  }
}

void PrimarySample::Mutate(tgir::PrimarySample const &prev, tgir::Real const r2,
                           std::mt19937_64 &random) {
  static tgir::Real const r3 = -std::log(tgir::Real(16));  // r1 = r2 / 16
  static tgir::Real const r4 = hi::rcp<tgir::Real>(32);  // インデックス選択用
  {
    muWavelength_.first =
        ::Mutate(prev.muWavelength_.first, hi::rcp<tgir::Real>(16), r3, random);
    muWavelength_.second = static_cast<std::size_t>(muWavelength_.first *
                                                    tgir::SpectrumData::Size);
  }
  for (std::size_t i = 0, size = muLight_.size(); i < size; ++i) {
    muLight_[i][0] = ::Mutate(prev.muLight_[i][0], r2, r3, random);
    muLight_[i][1] = ::Mutate(prev.muLight_[i][1], r2, r3, random);
    muLight_[i][2] = ::Mutate(prev.muLight_[i][2], r4, r3, random);
  }
  for (std::size_t i = 0, size = muEye_.size(); i < size; ++i) {
    muEye_[i][0] = ::Mutate(prev.muEye_[i][0], r2, r3, random);
    muEye_[i][1] = ::Mutate(prev.muEye_[i][1], r2, r3, random);
    muEye_[i][2] = ::Mutate(prev.muEye_[i][2], r4, r3, random);
  }
}

void PrimarySample::MutateFocusOnLightSubPath(tgir::PrimarySample const &prev,
                                              tgir::Real const r2,
                                              std::mt19937_64 &random) {
  static tgir::Real const r3 = -std::log(tgir::Real(16));  // r1 = r2 / 16
  static tgir::Real const r4 = hi::rcp<tgir::Real>(32);  // インデックス選択用
  {
    muWavelength_.first =
        ::Mutate(prev.muWavelength_.first, hi::rcp<tgir::Real>(16), r3, random);
    muWavelength_.second = static_cast<std::size_t>(muWavelength_.first *
                                                    tgir::SpectrumData::Size);
  }
  for (std::size_t i = 0, size = muLight_.size(); i < size; ++i) {
    muLight_[i][0] = ::Mutate(prev.muLight_[i][0], r2, r3, random);
    muLight_[i][1] = ::Mutate(prev.muLight_[i][1], r2, r3, random);
    muLight_[i][2] = ::Mutate(prev.muLight_[i][2], r4, r3, random);
  }
  for (std::size_t i = 0, size = muEye_.size(); i < size; ++i) {
    muEye_[i][0] = std::uniform_real_distribution<tgir::Real>()(random);
    muEye_[i][1] = std::uniform_real_distribution<tgir::Real>()(random);
    muEye_[i][2] = std::uniform_real_distribution<tgir::Real>()(random);
  }
}

void PrimarySample::MutateFocusOnEyeSubPath(tgir::PrimarySample const &prev,
                                            tgir::Real const r2,
                                            std::mt19937_64 &random) {
  static tgir::Real const r3 = -std::log(tgir::Real(16));  // r1 = r2 / 16
  static tgir::Real const r4 = hi::rcp<tgir::Real>(32);  // インデックス選択用
  {
    muWavelength_.first =
        ::Mutate(prev.muWavelength_.first, hi::rcp<tgir::Real>(16), r3, random);
    muWavelength_.second = static_cast<std::size_t>(muWavelength_.first *
                                                    tgir::SpectrumData::Size);
  }
  for (std::size_t i = 0, size = muLight_.size(); i < size; ++i) {
    muLight_[i][0] = std::uniform_real_distribution<tgir::Real>()(random);
    muLight_[i][1] = std::uniform_real_distribution<tgir::Real>()(random);
    muLight_[i][2] = std::uniform_real_distribution<tgir::Real>()(random);
  }
  for (std::size_t i = 0, size = muEye_.size(); i < size; ++i) {
    muEye_[i][0] = ::Mutate(prev.muEye_[i][0], r2, r3, random);
    muEye_[i][1] = ::Mutate(prev.muEye_[i][1], r2, r3, random);
    muEye_[i][2] = ::Mutate(prev.muEye_[i][2], r4, r3, random);
  }
}

#if 1
void PrimarySample::MutateLightSubPath(tgir::PrimarySample const &prev,
                                       tgir::Real const r2,
                                       std::mt19937_64 &random) {
  static tgir::Real const r3 = -std::log(tgir::Real(16));  // r1 = r2 / 16
  static tgir::Real const r4 = hi::rcp<tgir::Real>(32);  // インデックス選択用
  { muWavelength_ = prev.muWavelength_; }
  for (std::size_t i = 0, size = muLight_.size(); i < size; ++i) {
    muLight_[i][0] = ::Mutate(prev.muLight_[i][0], r2, r3, random);
    muLight_[i][1] = ::Mutate(prev.muLight_[i][1], r2, r3, random);
    muLight_[i][2] = ::Mutate(prev.muLight_[i][2], r4, r3, random);
  }
  for (std::size_t i = 0, size = muEye_.size(); i < size; ++i) {
    muEye_[i] = prev.muEye_[i];
  }
}

void PrimarySample::MutateEyeSubPath(tgir::PrimarySample const &prev,
                                     tgir::Real const r2,
                                     std::mt19937_64 &random) {
  static tgir::Real const r3 = -std::log(tgir::Real(16));  // r1 = r2 / 16
  static tgir::Real const r4 = hi::rcp<tgir::Real>(32);  // インデックス選択用
  { muWavelength_ = prev.muWavelength_; }
  for (std::size_t i = 0, size = muLight_.size(); i < size; ++i) {
    muLight_[i] = prev.muLight_[i];
  }
  for (std::size_t i = 0, size = muEye_.size(); i < size; ++i) {
    muEye_[i][0] = ::Mutate(prev.muEye_[i][0], r2, r3, random);
    muEye_[i][1] = ::Mutate(prev.muEye_[i][1], r2, r3, random);
    muEye_[i][2] = ::Mutate(prev.muEye_[i][2], r4, r3, random);
  }
}
#endif
}  // end of namespace tgir
