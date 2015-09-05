#ifndef __TGIR_GEOM_PATH_HPP__
#define __TGIR_GEOM_PATH_HPP__

#include "core/config.hpp"
#include "geom/Triangle.hpp"
#include "core/Film.hpp"
#include "geom/PathVertex.hpp"

namespace tgir {
class Path {
 public:
  static std::size_t const INVALID_CONSTRIBUTION = ~0U;

  explicit Path(__in std::size_t const maxPathLength);
  explicit Path(__in Path const &);

  Path &operator=(__in Path const &rhs);

  void swap(__inout Path &rhs);

  inline std::size_t GetMaxPathLength() const { return maxPathLength_; }

  void Begin(__in enum EMutationType);
  void End();
  void SetFilmPosition(__in std::size_t const &p);

  void Values(__out std::vector<tgir::PixelDescriptor> *const pValues) const;
  void ValuesAndSpecial(__out std::vector<tgir::PixelDescriptor> *const pValues,
                        __out tgir::PixelDescriptor *const pSpecial) const;

  tgir::Real PDFSimplifiedMLT() const;
  tgir::Real PDFFullMLT() const;
  void PDFsFullRELT(__out tgir::Vector4 *const pPDFs) const;

  /// ���q�ǐՖ@�ɂ���^(C_{s,1})��~�ς���.
  void AccumulateS1(
      __in std::size_t const &p,                 ///< [in] ��f�̈ʒu
      __in std::size_t const &uWavelengthIndex,  ///< [in] �g���̈ʒu
      __in std::size_t const &k,                 ///< [in] �o�H��-1
      __in tgir::Spectrum const &cs1);           ///< [in]

  /// ���ڏƖ��v�Z���s�����ꍇ�̊�^(C_{1,t})��~�ς���.
  void Accumulate1T(__in std::size_t const &uWavelengthIndex,
                    __in tgir::Spectrum const &c1t);

  /// �o�����o�H�ǐՂ��s�����ꍇ�̊�^(C_{s,t})��~�ς���.
  void AccumulateST(__in std::size_t const &uWavelengthIndex,
                    __in tgir::Spectrum const &cst);

  /// ���ڌ����������Ă���ꍇ�̊�^(C_{0,t})��~�ς���.
  void Accumulate0T(__in std::size_t const &uWavelengthIndex,
                    __in tgir::Spectrum const &c0t, __in bool const &bSpecial);

 public:
#ifndef CONFIG_PINHOLE_CAMERA_MODEL
  tgir::Vector2 vLensCoordinate;
#endif CONFIG_PINHOLE_CAMERA_MODEL

  std::vector<tgir::PathVertex> LightPathVertices;  ///< �����o�H���_
  std::vector<tgir::PathVertex> EyePathVertices;    ///< ���_�o�H���_

 private:
  std::size_t maxPathLength_;

  std::vector<tgir::PixelDescriptor>
      contributions_;  ///< �e�o�H���̌o�H����̊�^

  tgir::SpectrumVector contribution_from_particle_path_;  ///< S1
  tgir::SpectrumVector contribution_from_indirect_path_;  ///< ST
  tgir::SpectrumVector contribution_from_explicit_path_;  ///< 1T
  tgir::SpectrumVector contribution_from_implicit_path_;  ///< 0T
  tgir::SpectrumVector contribution_from_specific_path_;  ///< 0T (special)
};

}  // end of namespace tgir

namespace std {
template <>
inline void swap(tgir::Path &a, tgir::Path &b) {
  a.swap(b);
}
}

#endif __TGIR_GEOM_PATH_HPP__
