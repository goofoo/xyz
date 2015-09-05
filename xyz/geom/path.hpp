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

  /// ó±éqí«ê’ñ@Ç…ÇÊÇÈäÒó^(C_{s,1})Çí~êœÇ∑ÇÈ.
  void AccumulateS1(__in std::size_t const &p,  ///< [in] âÊëfÇÃà íu
                    __in std::size_t const &k,  ///< [in] åoòHí∑-1
                    __in float3_t const &cs1);  ///< [in]

  /// íºê⁄è∆ñæåvéZÇçsÇ¡ÇΩèÍçáÇÃäÒó^(C_{1,t})Çí~êœÇ∑ÇÈ.
  void Accumulate1T(__in float3_t const &c1t);

  /// ëoï˚å¸åoòHí«ê’ÇçsÇ¡ÇΩèÍçáÇÃäÒó^(C_{s,t})Çí~êœÇ∑ÇÈ.
  void AccumulateST(__in float3_t const &cst);

  /// íºê⁄åıåπÇ™å©Ç¶ÇƒÇ¢ÇÈèÍçáÇÃäÒó^(C_{0,t})Çí~êœÇ∑ÇÈ.
  void Accumulate0T(__in float3_t const &c0t, __in bool const &bSpecial);

 public:
#ifndef CONFIG_PINHOLE_CAMERA_MODEL
  float2_t vLensCoordinate;
#endif CONFIG_PINHOLE_CAMERA_MODEL

  std::vector<PathVertex> LightPathVertices;  ///< åıåπåoòHí∏ì_
  std::vector<PathVertex> EyePathVertices;    ///< éãì_åoòHí∏ì_

 private:
  std::size_t maxPathLength_;

  std::vector<pixel_descriptor_t> contributions_;  ///< äeåoòHí∑ÇÃåoòHÇ©ÇÁÇÃäÒó^

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
