#pragma once

#include "core/config.hpp"
#include "geom/Triangle.hpp"
#include "core/Film.hpp"
#include "geom/PathVertex.hpp"

namespace tgir {
class Path {
 public:
  static std::size_t const INVALID_CONSTRIBUTION = ~0U;

  explicit Path(Path const &);
  explicit Path(std::size_t const maxPathLength);

  inline std::size_t GetMaxPathLength() const { return maxPathLength_; }

  inline tgir::Real PDFSimplifiedMLT() const {
    return hi::sum(contribution_from_particle_path_) +
           hi::sum(contribution_from_indirect_path_) +
           hi::sum(contribution_from_explicit_path_) +
           hi::sum(contribution_from_implicit_path_) +
           hi::sum(contribution_from_specific_path_);
  }

  inline tgir::Real PDFFullMLT() const {
    return hi::sum(contribution_from_indirect_path_) +
           hi::sum(contribution_from_implicit_path_) +
           hi::sum(contribution_from_specific_path_);
  }

#ifdef TGIR_CONFIG_TRIPLE_REPLICS
  inline void PDFsFullRELT(__out tgir::Vector3 *const pPDFs) const
#else
  inline void PDFsFullRELT(__out tgir::Vector4 *const pPDFs) const
#endif
  {
    tgir::Real const a = hi::sum(contribution_from_particle_path_);
    tgir::Real const b = hi::sum(contribution_from_indirect_path_);
    tgir::Real const c = hi::sum(contribution_from_explicit_path_);
    tgir::Real const d = hi::sum(contribution_from_implicit_path_);
    tgir::Real const e = hi::sum(contribution_from_specific_path_);

#ifdef TGIR_CONFIG_TRIPLE_REPLICS
    (*pPDFs)[0] = b + d + e;
    (*pPDFs)[1] = a + b + c + d + e;
    (*pPDFs)[2] = 1;
#else
    (*pPDFs)[0] = b + e;
    (*pPDFs)[1] = b + d + e;
    (*pPDFs)[2] = a + b + c + d + e;
    (*pPDFs)[3] = 1;
#endif
  }

  void Begin(enum EMutationType);
  void End(std::vector<tgir::PixelDescriptor> &f);
  void SetFilmPosition(std::size_t const p);

  // ���q�ǐՖ@�ɂ���^(C_{s,1})��~�ς���
  void AccumulateS1(std::size_t const, std::size_t const, std::size_t const,
                    tgir::Spectrum const);

  // ���ڏƖ��v�Z���s�����ꍇ�̊�^(C_{1,t})��~�ς���
  void Accumulate1T(std::size_t const, tgir::Spectrum const);

  // �o�����o�H�ǐՂ��s�����ꍇ�̊�^(C_{s,t})��~�ς���
  void AccumulateST(std::size_t const, tgir::Spectrum const);

  // ���ڌ����������Ă���ꍇ�̊�^(C_{0,t})��~�ς���
  void Accumulate0T(std::size_t const, tgir::Spectrum const, bool const);

#if !defined(CONFIG_PINHOLE_CAMERA_MODEL)
  tgir::Vector2 vLensCoordinate;
#endif  // !defined(CONFIG_PINHOLE_CAMERA_MODEL)

 public:
  std::vector<tgir::PathVertex> LightPathVertices;  // �����o�H���_
  std::vector<tgir::PathVertex> EyePathVertices;    // ���_�o�H���_

 private:
  std::size_t const maxPathLength_;

  std::vector<tgir::PixelDescriptor>
      contributions_;  // �e�o�H���̌o�H����̊�^

  tgir::SpectrumVector contribution_from_particle_path_;  // S1
  tgir::SpectrumVector contribution_from_indirect_path_;  // ST
  tgir::SpectrumVector contribution_from_explicit_path_;  // 1T
  tgir::SpectrumVector contribution_from_implicit_path_;  // 0T
  tgir::SpectrumVector contribution_from_specific_path_;  // 0T (special)

  // ����`�̊֐�
  Path &operator=(Path const &);
};

}  // end of namespace tgir
