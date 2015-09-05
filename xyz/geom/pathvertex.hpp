#ifndef XYZ_GEOM_PATHVERTEX_HPP_
#define XYZ_GEOM_PATHVERTEX_HPP_

#include "../core/config.hpp"
#include "../geom/triangle.hpp"

namespace xyz {
class Intersection;

class PathVertex {
 public:
  PathVertex() : pGeometry(reinterpret_cast<Triangle const *>(0xDEADBEAF)) {}

  Triangle const *pGeometry;  ///< �O�p�`

  bool bBackSide;  ///< ��_�͖ʂ̗���
  bool bSpecular;  ///< �f���^�֐��ɂ�锽�˂�
  // bool bModified;

  float3_t vPosition;           ///< ��_�̍��W
  float3_t vIncomingDirection;  ///< ���˕���
  float3_t vOutgoingDirection;  ///< �ˏo����
  float3_t vGeometricNormal;    ///< �􉽊w�I�@��(geometric normal)
  float3_t vShadingNormal;      ///< �V�F�[�f�B���O�@��(shading noraml)
  float3_t vTangent;            ///< �ڃx�N�g��
  float3_t vBinormal;           ///< �]�@���x�N�g��

  float3_t power;  ///< importance (potential) or power

  float_t fBSDFxIPDF;  ///< BSDF/PDF �̒l

  float_t fSamplingPrev;  ///<
                          ///���̌o�H���_�����O�̌o�H���_���T���v�����O����ˉe���̊p�𑪓x�Ƃ����m��(���_��������̕���)
  float_t fSamplingNext;  ///<
                          ///���̌o�H���_�������̌o�H���_���T���v�����O����ˉe���̊p�𑪓x�Ƃ����m��(�������王�_�̕���)

  float_t fIncomingCosThetaShading;
  float_t fOutgoingCosThetaGeometric;

  float_t fGeometricFactor;  ///< ���̌o�H���_�ƈ�O�̌o�H���_�Ԃ̊􉽊w��
                             ///\frac{\cos\theta \cdot \cos\theta}{r^{2}}

  inline void SetGeometricBasis(Intersection const &param) {
    pGeometry->GeometricBasis(param, &vShadingNormal, &vGeometricNormal,
                              &vTangent, &vBinormal);
  }

  inline void SetShadingBasis(Intersection const &param) {
    pGeometry->ShadingBasis(param, &vShadingNormal, &vGeometricNormal,
                            &vTangent, &vBinormal);
  }

  bool SetGeometricBasis(Intersection const &param,
                         PathVertex const &oldLightPathVertex);
  bool SetShadingBasis(Intersection const &param,
                       PathVertex const &oldEyePathVertex);
};
}
// end of namespace xyz

#endif XYZ_GEOM_PATHVERTEX_HPP_
