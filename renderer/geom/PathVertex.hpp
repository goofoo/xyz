#ifndef __TGIR_GEOM_PATHVERTEX_HPP__
#define __TGIR_GEOM_PATHVERTEX_HPP__

#include "core/config.hpp"
#include "geom/Triangle.hpp"

namespace tgir {
class Intersection;

class PathVertex {
 public:
  tgir::Triangle const *pGeometry;  ///< �O�p�`

  tgir::Vector3 vPosition;           ///< ��_�̍��W
  tgir::Vector3 vIncomingDirection;  ///< ���˕���
  tgir::Vector3 vOutgoingDirection;  ///< �ˏo����
  tgir::Vector3 vGeometricNormal;    ///< �􉽊w�I�@��(geometric normal)
  tgir::Vector3 vShadingNormal;      ///< �V�F�[�f�B���O�@��(shading noraml)
  tgir::Vector3 vTangent;            ///< �ڃx�N�g��
  tgir::Vector3 vBinormal;           ///< �]�@���x�N�g��

  tgir::Spectrum sQuantum;  ///< importance (potential) or power

  bool back_side;  ///< ��_�͖ʂ̗���
  bool bSpecular;  ///< �f���^�֐��ɂ�锽�˂�
  // bool bModified;

  tgir::Real rSamplingPrev;  ///<
  ///���̌o�H���_�����O�̌o�H���_���T���v�����O����ˉe���̊p�𑪓x�Ƃ����m��
  tgir::Real rSamplingNext;  ///<
  ///���̌o�H���_�������̌o�H���_���T���v�����O����ˉe���̊p�𑪓x�Ƃ����m��
  tgir::Real rGeometricFactor;  ///< ���̌o�H���_�ƈ�O�̌o�H���_�Ԃ̊􉽊w��
                                ///\frac{\cos\theta \cdot \cos\theta}{r^{2}}

  inline void SetGeometricBasis(tgir::Intersection const &param) {
    pGeometry->GeometricBasis(param, &vShadingNormal, &vGeometricNormal,
                              &vTangent, &vBinormal);
  }

  inline void SetShadingBasis(tgir::Intersection const &param) {
    pGeometry->ShadingBasis(param, &vShadingNormal, &vGeometricNormal,
                            &vTangent, &vBinormal);
  }

  bool SetGeometricBasis(tgir::Intersection const &param,
                         tgir::PathVertex const &oldLightPathVertex);
  bool SetShadingBasis(tgir::Intersection const &param,
                       tgir::PathVertex const &oldEyePathVertex);
};
}
// end of namespace tgir

#endif __TGIR_GEOM_PATHVERTEX_HPP__
