#ifndef XYZ_BSDF_BSDF_HPP_
#define XYZ_BSDF_BSDF_HPP_

#include "../core/PhotonMap.hpp"

namespace xyz {
class PathVertex;
class Scene;

class IPrimarySample {
 public:
  virtual ~IPrimarySample() {}
  virtual float_t next() = 0;
};

//
// base of non-symmetric BSDF
//
class Bsdf {
 public:
  enum type {
    LAMBERT,
    MIRROR,
    GLASS,
    METAL,
    LIGHT,
    COUNT,
  };

  virtual void Render() const = 0;

  //
  // �t�H�g���}�b�v�\�z�p
  //
  /// �����܂��͎��_����̏d�v�x��ǐՂ���
  /// �t�H�g����ǐՂ���
  virtual bool BoundsPhoton(PathVertex &lightPathVertex,
                            std::mt19937_64 &random) const;
  /// ���_����̏d�v�x���g���ăt�H�g����ǐՂ���
  virtual bool BoundsPhotonWithImportonMap(PathVertex &lightPathVertex,
                                           ImportonMap const &importonMap,
                                           ImportonQuery &importons,
                                           std::mt19937_64 &random) const;
  /// ���q�t�B���^�ɂ�莋�_����̏d�v�x���g���ăt�H�g����ǐՂ���
  virtual bool BoundsPhotonWithParticleFilter(PathVertex &lightPathVertex,
                                              ImportonMap const &importonMap,
                                              ParticleQuery &importons,
                                              std::mt19937_64 &random) const;

  //
  // �o�H�ǐՖ@�p
  //
  virtual Bsdf::type What() const = 0;

  virtual float3_t CalculateWeightedDirectLighting(
      __inout PathVertex *const pEyePathVertex,
      __in PathVertex const &lightPathVertex,
      __in Scene const &scene) const = 0;

  virtual bool NextScatteringDirection(
      __inout PathVertex *const pEyePathVertex,
      __inout IPrimarySample &sample) const = 0;

  virtual float3_t CalculateWeightedDirectLighting(
      __inout PathVertex *const pEyePathVertex,
      __in PathVertex const &lightPathVertex,
      __in ImportanceQuery const &importons, __in Scene const &scene) const;

  // ��������̏d�v�x���g���Ĕ��˃x�N�g�����v�Z����
  virtual bool NextScatteringDirection(__inout PathVertex *const pEyePathVertex,
                                       __in ImportanceQuery const &importons,
                                       __inout IPrimarySample &sample) const;

  // �d�v�x�v�Z�p
  virtual bool BoundsImportance(__inout PathVertex *const pLightPathVertex,
                                __inout std::mt19937_64 &random) const = 0;

  // ���ˌ����̕�����
  virtual float_t GetDensityVariance() const = 0;

  //
  // �o�����o�H�ǐՖ@�p
  //
  virtual bool LightsToCameraScatteringDirection(
      __inout PathVertex *const pPathVertex,
      __inout IPrimarySample &sample) const = 0;

  virtual bool CameraToLightsScatteringDirection(
      __inout PathVertex *const pPathVertex,
      __inout IPrimarySample &sample) const = 0;

  virtual bool CameraToLightsScatteringDirection(
      __inout PathVertex *const pPathVertex,
      __in ImportanceQuery const &importons,
      __inout IPrimarySample &sample) const;

  //
  // �o�����o�H�ǐ՗p
  //
  virtual bool NextScatteringDirection(PathVertex &,      ///< [in/out]
                                       float3_t const &,  ///< [in]
                                       bool const &       ///< [in]
                                       ) const = 0;

  //
  // �܂��g���ĂȂ����ǎ������Ă���
  //
  // Get the density with respect to projected solid angle: p(w)
  // If p'(w) is the density with respect to ordinary solid angle, then p(w) =
  // p'(w) / |cos\theta|,
  //   where \theta is the angle between the outgoing direction and the surface
  //   normal.
  virtual bool HasNoDiracDistribution() const = 0;
  virtual bool IsLightEmittor() const = 0;

  virtual float_t BSDF(std::size_t const &is,  ///< [in] index of spectrum
                       float3_t const &wi,     ///< [in] incoming direction
                       float3_t const &wo,     ///< [in] outgoing direction
                       float3_t const &ns      ///< [in] shading normal
                       ) const = 0;
};

}  // end of namespace xyz

/*
���_������������ւ̃T���v�����O(fSamplingPrev)�ł́C
�V�F�[�f�B���O�@���Ɠ��˕����x�N�g���̓���(fIncomingCosThetaShading)�ŁC
�T���v�����O���z������

�������王�_�����ւ̃T���v�����O(fSamplingNext)�ł́C
�􉽊w�I�@���ƎU�������x�N�g���̓���(fOutgoingCosThetaGeometric)�ŁC
�T���v�����O���z������
*/

#endif XYZ_BSDF_BSDF_HPP_
