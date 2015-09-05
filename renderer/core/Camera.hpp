#ifndef __TGIR_CORE_CAMERA_HPP__
#define __TGIR_CORE_CAMERA_HPP__

namespace tgir {
class PrimarySample;

/**
 * The thin-lens Camera
 * ���w�藝: 1/L + 1/V = 1/f
 *   L: �����Y����s���g�������ʂ܂ł̋���
 *   V: �t�B�����ʂ��烌���Y�܂ł̋���
 *   f: �œ_����(�����Y�̒��S�ʒu���畽�s�������œ_�����Ԉʒu�܂ł̋���)
 * F�l: F = f/2R
 *   R: �����Y�̔��a(=f/2F)
 */
class Camera {
 public:
  Camera();

  //
  // GUI
  //
 public:
  void Render() const;

  inline tgir::Vector3 Eye() const { return vPrincipalPoint_; }
  inline tgir::Vector3 At() const { return vGazingPoint_; }

  void MoveTo(tgir::Vector3 const &v);
  void MoveBy(tgir::Vector3 const &v);
  void LookAt(tgir::Vector3 const &v);
  void LookBy(tgir::Vector3 const &v);

  inline void SetAspectRatio(tgir::Real const &aspect_ratio) {
    fAspectRatio_ = aspect_ratio;
  }

  inline tgir::Real GetLensArea() const {
    return hi::square_of(fLensRadius_) * M_PI;
  }

  /**
   * �s���z�[���J�������f���ɂ����鏉�������x�N�g����
   * �T���v�����O�p�̊m�����x���z(���xds(w)) p(w) �ɂ��āD
   *
   * �t�B�������ʏ�̃T���v�����O�p�̊m�����x���z�֐�(���xdA(x,y))��
   *   p(x,y) = 1/(W*H) ... (�t�B�������ʏ�̈�l���z�̏ꍇ)
   * �ł���D�����ŁCw �͐��E���W�n�ɂ�����t�B�����̕�(m�P��)�ŁC
   * ���l�� h �̓t�B�����̍���(������m�P��)�ł���D
   * �܂�C�t�B�����ʏ�̂�������̈� dA(x,y) ���T���v�����O�����m����
   *   p(x,y) dA(x,y) = 1/(W*H) dA(x,y)
   * �ƂȂ�D
   *
   * ���ʂɂ���������̃T���v�����O�̑��x(ds(w))����
   * ���ʂɂ�����_�̃T���v�����O�̑��x(dA(x))�ւ̕ϊ����́C
   * ���̊p�̒�`���� ds(w) = cos(theta)/r^2 dA(x) �ƂȂ�D
   * �����ŁCtheta �́C���� w �ƕ��ʂ̖@���Ƃ̊p�x�ł���D
   * �܂��CD �́C�����̃T���v�����O�̎n�_����
   * ���ʏ�̃T���v�����O�_�܂ł̋����ł���D
   *
   * �܂�C�t�B�������ʏ�ň�l���z�́C�����Ɋւ���m�����x�֐��ɕϊ������
   *   p(x,y) dA(x,y)
   *     = 1/(W*H) dA(x,y)                 ...
   *�T���v�����O�p�̊m�����x�֐��͕��ʏ�̈�l���z
   *     = 1/(W*H) D^2/cos(theta) ds(w)    ... ���x��ʐϑ��x����������x�ɕϊ�
   *     = V^2/(W*H) 1/cos^3(theta) ds(w)  ... V = D * cos(theta)
   *�������āC�ϐ����ƒ萔�����܂Ƃ߂�
   * �����ŁCV �̓����Y�E�t�B�����Ԃ̋����ł���D
   *
   * p(x,y) dA(x,y) = p(w) ds(w) ���C
   * p(w) = V^2/(W*H) 1/cos^3(theta) �ł���C
   * �ˉe�m�����x�֐��� p��(w) = d^2/(W*H) 1/cos^4(theta)
   * �ƂȂ�D
   *
   * ���̊֐��́C���̒萔�� V^2/(W*H) ��Ԃ��D
   *
   * ���R�C�t�H�[�J�X�������ʂŃT���v�����O�����ꍇ�������l�ɂȂ�D
   * W' = (L/V)*W  ... �t�H�[�J�X�������ʂ̕�
   * H' = (L/V)*H  ... �t�H�[�J�X�������ʂ̍���
   *
   * p(w) = L^2/(W'*H') 1/cos^3(theta)
   *      = L^2/((L/V)^2*W*H) 1/cos^3(theta)
   *      = V^2/(W*H) 1/cos^3(theta)
   */
  inline tgir::Real GetConstFactor() const {
    return hi::square_of(fFocalLength_ / fFilmSize_) / fAspectRatio_;
  }

 private:
  void Reset();

  //
  // renderer
  //
 public:
  //==========================================================================
  //
  // Pinhole camera model.
  //
  //==========================================================================

  void GetPrimaryRayOriginAndCameraBasis(
      __out tgir::Vector3
          *const pvRayOrigin,  ///< �����̎n�_(�O������ԓ��̃����Y��̈ʒu)
      __out tgir::Vector3 *const pvRightwardVector,  ///< �J�����̉E�����x�N�g��
      __out tgir::Vector3 *const pvUpwardVector,     ///< �J�����̏�����x�N�g��
      __out tgir::Vector3 *const pvBackwardVector    ///< �J�����̌�����x�N�g��
      ) const;

  void GetPrimaryRayDirection(
      __in tgir::Real const &s,  ///< �t�B�����ʏ�̈ʒu(X)�D�l���[0,1]
      __in tgir::Real const &t,  ///< �t�B�����ʏ�̈ʒu(Y)�D�l���[0,1]
      __out tgir::Vector3 *const pvRayDirection  ///< �ꎟ�����̕���
      ) const;

  bool GetFilmPosition(
      __in tgir::Vector3 const &vDirection,  ///< [in]
      ///�����Y��̍��W����O������ԓ��̕��̕\�ʏ�̓_�ւ̃x�N�g��
      __out tgir::Vector2
          *const pvFilmPosition  ///< [out] �t�B������̍��W(X��[0,1],Y��[0,1])
      ) const;

  //==========================================================================
  //
  // Thin-lens camera model.
  //
  //==========================================================================

  void GetPrimaryRayOriginAndCameraBasis(
      __in tgir::Real const &s,  ///< [in]
      ///�����Y��̍��W���v�Z���邽�߂̃p�����[�^(���a�p)�D�l���[0,1]
      __in tgir::Real const &t,  ///< [in]
      ///�����Y��̍��W���v�Z���邽�߂̃p�����[�^(�p�x�p)�D�l���[0,1]
      __out tgir::Vector2
          *const vLensCoordinate,  ///< [out] �����Y��̍��W(X,Y)
      __out tgir::Vector3 *const
          vRayOrigin,  ///< [out] �����̎n�_(�O������ԓ��̃����Y��̈ʒu)
      __out tgir::Vector3
          *const vRightwardVector,  ///< [out] �J�����̉E�����x�N�g��
      __out tgir::Vector3
          *const vUpwardVector,  ///< [out] �J�����̏�����x�N�g��
      __out tgir::Vector3
          *const vBackwardVector  ///< [out] �J�����̌�����x�N�g��
      ) const;

  void GetPrimaryRayDirection(
      __in tgir::Real const &s,  ///< [in] �t�B�����ʏ�̈ʒu(X)�D�l���[0,1]
      __in tgir::Real const &t,  ///< [in] �t�B�����ʏ�̈ʒu(Y)�D�l���[0,1]
      __in tgir::Vector2 const &vLensCoordinate,  ///< [in] �����Y��̍��W(X,Y)
      __out tgir::Vector3 *const pvRayDirection   ///< [out] �ꎟ�����̕���
      ) const;

  bool GetFilmPosition(
      __in tgir::Vector2 const &vLensCoordinate,  ///< [in] �����Y��̍��W(X,Y)
      __in tgir::Vector3 const &vDirection,       ///< [in]
      ///�����Y��̍��W����O������ԓ��̕��̕\�ʏ�̓_�ւ̃x�N�g��
      __out tgir::Vector2
          *const pvFilmPosition  ///< [out] �t�B������̍��W(X��[0,1],Y��[0,1])
      ) const;

 private:
  tgir::Vector3 vPrincipalPoint_;  ///< �����Y�̎�_�̍��W
  tgir::Vector3 vGazingPoint_;     ///< �����_

  tgir::Real fFNumber_;      ///< �����Y����: F�l
  tgir::Real fFocalLength_;  ///< �����Y����: �œ_����(f)

  tgir::Real fFilmSize_;     ///< �t�B��������: �t�B�����T�C�Y(����(H))
  tgir::Real fAspectRatio_;  ///< �t�B��������:
  ///�t�B�����̃A�X�y�N�g��(�t�B�����̕�(H)���t�B�����̍���(H)�~�A�X�y�N�g��(A))

  tgir::Real fFocusDistance_;  ///<
  ///�����Y���u���Ă��镽�ʂƃs���g���������ʂ̊Ԃ̋���(L)
  tgir::Real fLensFilmDistance_;  ///< �����Y�E�t�B�����Ԃ̋���(V)
  tgir::Real fLensRadius_;        ///< �����Y�̗L�����a(R)

  tgir::Vector3 vRightwardVector_;  ///< �J�����̉E�����x�N�g��
  tgir::Vector3 vUpwardVector_;     ///< �J�����̏�����x�N�g��
  tgir::Vector3 vBackwardVector_;   ///< �J�����̌�����x�N�g��
};

}  // end of namespace tgir

#endif __TGIR_CORE_CAMERA_HPP__
