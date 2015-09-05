#ifndef __TGIR_CORE_CONFIG_HPP__
#define __TGIR_CORE_CONFIG_HPP__

#define RENDERER_NAME _TEXT("Tiny Global Illumination Renderer")
#define RENDERER_BUILD _TEXT(__DATE__)

// __FILE__ �t�@�C����
// __TIMESTAMP__ �t�@�C���̍ŏI�X�V����
// __TIME__ �R���p�C������
// __DATE__ �R���p�C����
// __LINE__ �s�ԍ�
// __FUNCTION__ �֐���(C99��__func__)
// __cplusplus C++�ŃR���p�C������Ă���

#if 0
#define MQODIR_PATH _TEXT("C:/home/res/mqo/")
#define RENDIR_PATH _TEXT("C:/home/res/ren/")
#define OUTDIR_PATH _TEXT("C:/home/res/out/")
#else
#define MQODIR_PATH _TEXT("./")
#define RENDIR_PATH _TEXT("./")
#define OUTDIR_PATH _TEXT("./")
#endif

// �S�ẴA���S���Y���ɋ��ʂ���p�����[�^
#define TGIR_CONFIG_kMaxPathLength 8  ///< �o�����o�H�ǐՖ@�Ő�������ő�o�H��
#define TGIR_CONFIG_kSampleCount 40000  ///< �P���[�v�̃T���v����
#define TGIR_CONFIG_kSystemCount 16     ///< �ʂ̃T���v���[�̐�

#define TGIR_CONFIG_kLargeStepProb \
  0.25  ///< Simplified MLT �Ń��[�W�X�e�b�v�ψق��s���m��

// �o�H�ǐՖ@�Ɋւ���萔
#define TGIR_CONFIG_kMaxRandomWalkDepth 8  ///< �����_���E�H�[�N�ł̍ő�ǐՐ[�x
#define TGIR_CONFIG_kMaxBranchCount 10     ///< �����_���E�H�[�N�ł̍ő啪��
#define TGIR_CONFIG_kStratifiedSize 100    ///< Population Annealing �ł̑w����
#define CONFIG_SEEDS             \
  (TGIR_CONFIG_kStratifiedSize * \
   TGIR_CONFIG_kStratifiedSize)  // �ϕ��l�̐���̂��߂̃T���v����

#define CONFIG_BLOCK_SIZE 4  // ���q�t�B���^�̃t�B�����w���T�C�Y
#define CONFIG_LIGHT_SIZE 1  // ���q�t�B���^�̌����w���T�C�Y
#define CONFIG_COLOR_SIZE 1  // ���q�t�B���^�̃X�y�N�g���w���T�C�Y
#define CONFIG_FILTER_SIZE                                     \
  (CONFIG_BLOCK_SIZE * CONFIG_BLOCK_SIZE * CONFIG_LIGHT_SIZE * \
   CONFIG_LIGHT_SIZE * CONFIG_COLOR_SIZE)
#define CONFIG_SS_TYPE 3

// �t�H�g���}�b�v�p�̒萔
#define CONFIG_IMPORTON_STRATIFY 1024  // �C���|�[�g���̊e�����̑w�ʉ���
#define CONFIG_IMPORTON_SIZE  \
  (CONFIG_IMPORTON_STRATIFY * \
   CONFIG_IMPORTON_STRATIFY)   // �ˏo����C���|�[�g����
#define CONFIG_U_SIZE 4        //
#define CONFIG_V_SIZE 16       //
#define CONFIG_BUNDLE_SIZE 25  // �����ɒǐՂ���t�H�g����
#define CONFIG_N_SIZE (CONFIG_U_SIZE * CONFIG_V_SIZE)  //
#define CONFIG_M_SIZE (CONFIG_N_SIZE * CONFIG_BUNDLE_SIZE)
#define CONFIG_LIGHT_POWER 100000
#define CONFIG_SEARCH_RADIUS 0.2  // �t�H�g����C���|�[�g���̒T�����a

#define CONFIG_ALBEDO 0.65        // �g�U���˗�
#define CONFIG_EMITED_SIZE 10000  // �t�H�g���̎ˏo��
#define CONFIG_TERMINATE_THRESHOLD \
  0.2  // 0.2/0.5 // ���V�A�����[���b�g��臒l(sigma)

//#define CONFIG_PINHOLE_CAMERA_MODEL // �s���z�[���J�������ǂ����̃t���O
//#define CONFIG_DOUBLE_LIGHT // ���������ʗL�����ǂ����̃t���O
//#define CONFIG_BACKFACE_CULLING // �w�ʃJ�����O���s�����̃t���O
//#define CONFIG_RUSSIAN_RULLET //
//�o�H�ǐՂ̃��V�A�����[���b�g���s�����ǂ����̃t���O
#define CONFIG_MULTIPLE_INTEGRAION
#define CONFIG_CONCURRENT_LOCK
#define CONFIG_SPACE_SUBDIVISION

//#define CONFIG_HDR_STRICT

// ����̌o�H�������_�����O���邩�ǂ����̃t���O
#define CONFIG_EXPLICIT_GLOSSY  // ���򔽎˂Ŗ����I�Ȓ��ڏƖ��v�Z���s�����ǂ���
#define CONFIG_RENDER_CAUSTICS  // �o�H�ǐՖ@�ɂ�����W���͗l�o�H

//#define LIGHT_PROBE 900

namespace tgir {
inline tgir::SpectrumData const &GetLightPower() {
  return tgir::SpectrumData::LightMercuryArcLamp();
}

inline CIE_XYZ_Color LightMercuryArcLamp() {
  CIE_XYZ_Color value;
  hi::spd2xyz(GetLightPower(), value);
  return value;
}

inline CIE_XYZ_Color GetLightPowerXYZ() {
  static CIE_XYZ_Color const value(LightMercuryArcLamp());
  return value;
}

// �[���t���X�l����: fake fresnel term
inline tgir::Real FakeFresnelTerm(__in tgir::Real const hk,
                                  __in tgir::Real const fresnel0 = 0.9) {
  return fresnel0 + (1 - fresnel0) * hi::fifth_power_of(1 - hk);
}

// get random vector with importance sampleing (p(theta,phi) = cos(theta)/pi)
inline tgir::Real ComputeDiffusedVector(
    __in tgir::Real const fMuTheta, __in tgir::Real const fMuPhi,
    __in tgir::Vector3 const &vTangent, __in tgir::Vector3 const &vNormal,
    __in tgir::Vector3 const &vBinormal,
    __out tgir::Vector3 *const pvDiffusedVector) {
  tgir::Real const fPhi = tgir::Real(2 * M_PI) * fMuPhi;
  tgir::Real const fCosTheta = std::sqrt(fMuTheta);
  tgir::Real const fSinTheta = std::sqrt(1 - fMuTheta);
  *pvDiffusedVector = vTangent * (fSinTheta * std::cos(fPhi)) +
                      vBinormal * (fSinTheta * std::sin(fPhi)) +
                      vNormal * (fCosTheta);
  return fCosTheta;
}

inline tgir::Real ComputeDiffusedVector(
    __in tgir::Vector3 const &vMu, __in tgir::Vector3 const &vTangent,
    __in tgir::Vector3 const &vNormal, __in tgir::Vector3 const &vBinormal,
    __out tgir::Vector3 *const pvDiffusedVector) {
  return ComputeDiffusedVector(vMu[0], vMu[1], vTangent, vNormal, vBinormal,
                               pvDiffusedVector);
}

/// <summary>
/// �@���x�N�g���E�@���Ɠ��˕����x�N�g���̓��ρE���˕����x�N�g�����狾�ʔ��˃x�N�g�����v�Z����D
/// </summary>
/// <param name="vNormal">�@���x�N�g��</param>
/// <param name="rCosTheta">�@���Ɠ��˕����x�N�g���̓���</param>
/// <param
/// name="vReflectedVector">���͂Ƃ��Ă͓��˕����̃x�N�g��/�o�͂Ƃ��Ă͋��ʔ��˃x�N�g��</param>
/// <returns></returns>
/// <remarks>
/// �����Ɋ�Â��ċ��ʔ��˃x�N�g�����v�Z����D
/// <newpara>
/// vReflectedVector += 2 * rCosTheta * vNormalVector
/// </newpara>
/// </remarks>
inline void ComputeReflectedVector(
    __in tgir::Vector3 const &vNormalVector, __in tgir::Real const rCosTheta,
    __inout tgir::Vector3 *const vReflectedVector) {
  (*vReflectedVector) += (tgir::Real(2) * rCosTheta) * vNormalVector;
}

/// <summary>
/// ���΋��ܗ��Ɋ�Â������܃x�N�g�����v�Z����D
/// </summary>
/// <param name="vNormalVector">�@���x�N�g��</param>
/// <param name="gcm">�ڍׂ͌�q</param>
/// <param
/// name="rIndexOfRefraction">���΋��ܗ�(���ˌ�̔}���̋��ܗ�/���ˑO�̔}���̋��ܗ�)</param>
/// <param
/// name="vRefractedVector">���͂Ƃ��Ă͓��˕����̃x�N�g��/�o�͂Ƃ��Ă͑��΋��ܗ��Ɋ�Â������܃x�N�g��</param>
/// <returns></returns>
/// <remarks>
/// �����Ɋ�Â��ċ��܃x�N�g�����v�Z����D
/// <newpara>
/// vRefractedVector = (vRefractedVector - gmc * vNormalVector) /
/// rIndexOfRefraction
/// </newpara>
/// �����ŁCgmc�̒�`�͉����̒ʂ�ł���D
/// <newpara>
/// gmc = c - sqrt(rIndexOfRefraction^{2} + c^{2} - 1)
/// </newpara>
/// �������C
/// <newpara>
/// c = dot(vRefractedVector, vNormal)
/// </newpara>
/// �ł���D
/// </remarks>
inline void ComputeRefractedVector(__in tgir::Vector3 const &vNormalVector,
                                   __in tgir::Real const gmc,
                                   __in tgir::Real const rIndexOfRefraction,
                                   __inout tgir::Vector3 &vRefractedVector) {
  vRefractedVector =
      (vRefractedVector - gmc * vNormalVector) / rIndexOfRefraction;
}

/// <summary>
/// �m�����zph(vHalfwayVector)�Ɋ�Â��Ē���(halfway)�x�N�g�����v�Z����D
/// </summary>
/// <param
/// name="rMuTheta">���ʂ̃x�N�g���̈ܓx���������肷�邽�߂̗����l</param>
/// <param name="rMuPhi">���ʂ̃x�N�g���̌o�x���������肷�邽�߂̗����l</param>
/// <param
/// name="rShininess">���ʔ��˃��[�u�̌`������肷��ϐ�(�傫���قǉs�����˂ɂȂ�)</param>
/// <param name="vTangent">�ڃx�N�g��(���x�N�g����X���ɑ�������)</param>
/// <param name="vNormal">�@���x�N�g��(���x�N�g����Y���ɑ�������)</param>
/// <param name="vBinormal">�]�@���x�N�g��(���x�N�g����Z���ɑ�������)</param>
/// <param
/// name="vHalfwayVector">���͂Ƃ��Ă͓��˕����̃x�N�g��/�o�͂Ƃ��Ă͊m�����zph(vHalfwayVector)�Ɋ�Â������ԃx�N�g��</param>
/// <returns></returns>
/// <remarks>
/// �m�����zph(vHalfwayVector)�Ɋ�Â��Ē��ԃx�N�g�����v�Z����D
/// <newpara>
/// ph(vHalfwayVector) = (rShininess + 1) / 2pi * dot(vNormal,
/// vHalfwayVector)^{rShininess}
/// </newpara>
/// ���̒��ԃx�N�g�����g���Đ�������锽�˃x�N�g���́C�m�����zp(vGlossyVector)�ɏ]���D
/// <newpara>
/// p(vGlossyVector) = ph(vHalfwayVector) / (4 * dot(vIncomingVector,
/// vHalfwayVector))
/// </newpara>
/// </remarks>
inline tgir::Real ComputeGlossyHalfwayVector(
    __in tgir::Real const rMuTheta, __in tgir::Real const rMuPhi,
    __in tgir::Real const rShininess, __in tgir::Vector3 const &vTangent,
    __in tgir::Vector3 const &vNormal, __in tgir::Vector3 const &vBinormal,
    __out tgir::Vector3 *const pvHalfwayVector) {
  tgir::Real const rPhi = tgir::Real(2 * M_PI) * rMuPhi;
  tgir::Real const rCosTheta = std::pow(rMuTheta, hi::rcp(rShininess));
  tgir::Real const rSinTheta =
      std::sqrt(tgir::Real(1) - hi::square_of(rCosTheta));

  (*pvHalfwayVector) = vTangent * (rSinTheta * std::cos(rPhi)) +
                       vBinormal * (rSinTheta * std::sin(rPhi)) +
                       vNormal * (rCosTheta);

  return rCosTheta;
}

template <typename T>
inline T SafeDivide(T const &x, T const &y) {
  return (y < x) ? (y / x) : ((x > T(0)) ? T(1) : T(0));
}

}  // end of namespace tgir

#endif __TGIR_CORE_CONFIG_HPP__
