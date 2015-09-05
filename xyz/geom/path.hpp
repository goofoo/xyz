#ifndef XYZ_GEOM_PATH_HPP_
#define XYZ_GEOM_PATH_HPP_

#include "../core/config.hpp"
#include "../core/imagefilm.hpp"
#include "../geom/triangle.hpp"
#include "../geom/pathvertex.hpp"

namespace xyz {
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

  void Values(__out std::vector<pixel_descriptor_t> *const pValues) const;
  void ValuesAndSpecial(__out std::vector<pixel_descriptor_t> *const pValues,
                        __out pixel_descriptor_t *const pSpecial) const;

  float_t PDFSimplifiedMLT() const;
  float_t PDFFullMLT() const;
  void PDFsFullRELT(__out float4_t *const pPDFs) const;

  /// ���q�ǐՖ@�ɂ���^(C_{s,1})��~�ς���.
  void AccumulateS1(__in std::size_t const &p,  ///< [in] ��f�̈ʒu
                    __in std::size_t const &k,  ///< [in] �o�H��-1
                    __in float3_t const &cs1);  ///< [in]

  /// ���ڏƖ��v�Z���s�����ꍇ�̊�^(C_{1,t})��~�ς���.
  void Accumulate1T(__in float3_t const &c1t);

  /// �o�����o�H�ǐՂ��s�����ꍇ�̊�^(C_{s,t})��~�ς���.
  void AccumulateST(__in float3_t const &cst);

  /// ���ڌ����������Ă���ꍇ�̊�^(C_{0,t})��~�ς���.
  void Accumulate0T(__in float3_t const &c0t, __in bool const &bSpecial);

 public:
#ifndef CONFIG_PINHOLE_CAMERA_MODEL
  float2_t vLensCoordinate;
#endif CONFIG_PINHOLE_CAMERA_MODEL

  std::vector<PathVertex> LightPathVertices;  ///< �����o�H���_
  std::vector<PathVertex> EyePathVertices;    ///< ���_�o�H���_

 private:
  std::size_t maxPathLength_;

  std::vector<pixel_descriptor_t> contributions_;  ///< �e�o�H���̌o�H����̊�^

  float3_t contribution_from_particle_path_;  ///< S1
  float3_t contribution_from_indirect_path_;  ///< ST
  float3_t contribution_from_explicit_path_;  ///< 1T
  float3_t contribution_from_implicit_path_;  ///< 0T
  float3_t contribution_from_specific_path_;  ///< 0T (special)
};

}  // end of namespace xyz

namespace std {
template <>
inline void swap(xyz::Path &a, xyz::Path &b) {
  a.swap(b);
}
}

#endif XYZ_GEOM_PATH_HPP_
