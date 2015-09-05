#ifndef XYZ_THINLENSCAMERA_HPP_
#define XYZ_THINLENSCAMERA_HPP_

namespace xyz {
class PrimarySample;

/**
 * The thin-lens ThinLensCamera
 * ���w�藝: 1/L + 1/V = 1/f
 *   L: �����Y����s���g�������ʂ܂ł̋���
 *   V: �t�B�����ʂ��烌���Y�܂ł̋���
 *   f: �œ_����(�����Y�̒��S�ʒu���畽�s�������œ_�����Ԉʒu�܂ł̋���)
 * F�l: F = f/2R
 *   R: �����Y�̔��a(=f/2F)
 */
class ThinLensCamera {
 public:
  ThinLensCamera();

  //
  // GUI
  //
 public:
  void Render() const;

  inline float3_t Eye() const { return vPrincipalPoint_; }
  inline float3_t At() const { return vGazingPoint_; }

  void MoveTo(float3_t const &v);
  void MoveBy(float3_t const &v);
  void LookAt(float3_t const &v);
  void LookBy(float3_t const &v);

  inline void SetAspectRatio(float_t const &aspect_ratio) {
    fAspectRatio_ = aspect_ratio;
  }

  inline float_t GetLensArea() const {
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
  inline float_t GetConstFactor() const {
    return hi::square_of(fFocalLength_ / fFilmSize_) / fAspectRatio_;
  }

  inline float_t FluxToRadianceCoefficient() const { return GetConstFactor(); }

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
      __out float3_t
          *const pvRayOrigin,  ///< �����̎n�_(�O������ԓ��̃����Y��̈ʒu)
      __out float3_t *const pvRightwardVector,  ///< �J�����̉E�����x�N�g��
      __out float3_t *const pvUpwardVector,     ///< �J�����̏�����x�N�g��
      __out float3_t *const pvBackwardVector    ///< �J�����̌�����x�N�g��
      ) const;

  void GetPrimaryRayDirection(
      __in float_t const &s,  ///< �t�B�����ʏ�̈ʒu(X)�D�l���[0,1]
      __in float_t const &t,  ///< �t�B�����ʏ�̈ʒu(Y)�D�l���[0,1]
      __out float3_t *const pvRayDirection  ///< �ꎟ�����̕���
      ) const;

  bool GetFilmPosition(
      __in float3_t const &vDirection,  ///< [in]
      ///�����Y��̍��W����O������ԓ��̕��̕\�ʏ�̓_�ւ̃x�N�g��
      __out float2_t
          *const pvFilmPosition  ///< [out] �t�B������̍��W(X��[0,1],Y��[0,1])
      ) const;

  //==========================================================================
  //
  // Thin-lens camera model.
  //
  //==========================================================================

  void GetPrimaryRayOriginAndCameraBasis(
      __in float_t const &s,  ///< [in]
      ///�����Y��̍��W���v�Z���邽�߂̃p�����[�^(���a�p)�D�l���[0,1]
      __in float_t const &t,  ///< [in]
      ///�����Y��̍��W���v�Z���邽�߂̃p�����[�^(�p�x�p)�D�l���[0,1]
      __out float2_t *const vLensCoordinate,  ///< [out] �����Y��̍��W(X,Y)
      __out float3_t *const
          vRayOrigin,  ///< [out] �����̎n�_(�O������ԓ��̃����Y��̈ʒu)
      __out float3_t *const vRightwardVector,  ///< [out] �J�����̉E�����x�N�g��
      __out float3_t *const vUpwardVector,     ///< [out] �J�����̏�����x�N�g��
      __out float3_t *const vBackwardVector    ///< [out] �J�����̌�����x�N�g��
      ) const;

  void GetPrimaryRayDirection(
      __in float_t const &s,  ///< [in] �t�B�����ʏ�̈ʒu(X)�D�l���[0,1]
      __in float_t const &t,  ///< [in] �t�B�����ʏ�̈ʒu(Y)�D�l���[0,1]
      __in float2_t const &vLensCoordinate,  ///< [in] �����Y��̍��W(X,Y)
      __out float3_t *const pvRayDirection   ///< [out] �ꎟ�����̕���
      ) const;

  bool GetFilmPosition(
      __in float2_t const &vLensCoordinate,  ///< [in] �����Y��̍��W(X,Y)
      __in float3_t const &vDirection,       ///< [in]
      ///�����Y��̍��W����O������ԓ��̕��̕\�ʏ�̓_�ւ̃x�N�g��
      __out float2_t
          *const pvFilmPosition  ///< [out] �t�B������̍��W(X��[0,1],Y��[0,1])
      ) const;

 private:
  float3_t vPrincipalPoint_;  ///< �����Y�̎�_�̍��W
  float3_t vGazingPoint_;     ///< �����_

  float_t fFNumber_;      ///< �����Y����: F�l
  float_t fFocalLength_;  ///< �����Y����: �œ_����(f)

  float_t fFilmSize_;     ///< �t�B��������: �t�B�����T�C�Y(����(H))
  float_t fAspectRatio_;  ///< �t�B��������:
  ///�t�B�����̃A�X�y�N�g��(�t�B�����̕�(H)���t�B�����̍���(H)�~�A�X�y�N�g��(A))

  float_t fFocusDistance_;  ///<
  ///�����Y���u���Ă��镽�ʂƃs���g���������ʂ̊Ԃ̋���(L)
  float_t fLensFilmDistance_;  ///< �����Y�E�t�B�����Ԃ̋���(V)
  float_t fLensRadius_;        ///< �����Y�̗L�����a(R)

  float3_t vRightwardVector_;  ///< �J�����̉E�����x�N�g��
  float3_t vUpwardVector_;     ///< �J�����̏�����x�N�g��
  float3_t vBackwardVector_;   ///< �J�����̌�����x�N�g��
};

}  // end of namespace xyz

#endif XYZ_CAMERA_HPP_
