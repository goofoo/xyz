#ifndef XYZ_CONFIG_HPP_
#define XYZ_CONFIG_HPP_

#include "../stdafx.hpp"

#define RENDERER_NAME _TEXT("XYZ")
#define RENDERER_BUILD _TEXT(__DATE__)

// __FILE__ �t�@�C����
// __TIMESTAMP__ �t�@�C���̍ŏI�X�V����
// __TIME__ �R���p�C������
// __DATE__ �R���p�C����
// __LINE__ �s�ԍ�
// __FUNCTION__ �֐���(C99��__func__)
// __cplusplus C++�ŃR���p�C������Ă���

#define MQODIR_PATH _TEXT("./")
#define RENDIR_PATH _TEXT("./")
#define OUTDIR_PATH _TEXT("./")

// �S�ẴA���S���Y���ɋ��ʂ���p�����[�^
#define XYZ_CONFIG_kMaxPathLength 8    ///< �o�����o�H�ǐՖ@�Ő�������ő�o�H��
#define XYZ_CONFIG_kSampleCount 40000  ///< �P���[�v�̃T���v����
#define XYZ_CONFIG_kSystemCount 16     ///< �ʂ̃T���v���[�̐�

#define XYZ_CONFIG_kLargeStepProb \
  float_t(0.25)  ///< Simplified MLT �Ń��[�W�X�e�b�v�ψق��s���m��

// �o�H�ǐՖ@�Ɋւ���萔
#define XYZ_CONFIG_kMaxRandomWalkDepth 16  ///< �����_���E�H�[�N�ł̍ő�ǐՐ[�x
#define XYZ_CONFIG_kMaxBranchCount 10      ///< �����_���E�H�[�N�ł̍ő啪��
#define XYZ_CONFIG_kStratifiedSize 100     ///< Population Annealing �ł̑w����
#define CONFIG_SEEDS            \
  (XYZ_CONFIG_kStratifiedSize * \
   XYZ_CONFIG_kStratifiedSize)  // �ϕ��l�̐���̂��߂̃T���v����

#define CONFIG_BLOCK_SIZE 10  // ���q�t�B���^�̃t�B�����w���T�C�Y
#define CONFIG_LIGHT_SIZE 1   // ���q�t�B���^�̌����w���T�C�Y
#define CONFIG_COLOR_SIZE 1   // ���q�t�B���^�̃X�y�N�g���w���T�C�Y
#define CONFIG_FILTER_SIZE                                     \
  (CONFIG_BLOCK_SIZE * CONFIG_BLOCK_SIZE * CONFIG_LIGHT_SIZE * \
   CONFIG_LIGHT_SIZE * CONFIG_COLOR_SIZE)
#define CONFIG_SS_TYPE 3

// �t�H�g���}�b�v�p�̒萔
#define CONFIG_IMPORTON_STRATIFY 100  // �C���|�[�g���̊e�����̑w�ʉ���
#define CONFIG_IMPORTON_SIZE  \
  (CONFIG_IMPORTON_STRATIFY * \
   CONFIG_IMPORTON_STRATIFY)  // �ˏo����C���|�[�g����

#define CONFIG_BUNDLE_SIZE 25  // �����ɒǐՂ���t�H�g����
#define CONFIG_N_SIZE (CONFIG_U_SIZE * CONFIG_V_SIZE)  //
#define CONFIG_M_SIZE (CONFIG_N_SIZE * CONFIG_BUNDLE_SIZE)
#define CONFIG_LIGHT_POWER 100000

#define CONFIG_SEARCH_RADIUS float_t(0.1)  // �t�H�g����C���|�[�g���̒T�����a
#define CONFIG_EMITED_SIZE 100000          // �t�H�g���̎ˏo��

#define CONFIG_ALBEDO float_t(0.65)  // �g�U���˗�
#define CONFIG_TERMINATE_THRESHOLD \
  float_t(0.2)  // 0.2/0.5 // ���V�A�����[���b�g��臒l(sigma)

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

namespace xyz {
inline hi::basic_spectrum const &GetLightPower() {
  return hi::basic_spectrum::LightMercuryArcLamp();
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
}
#endif XYZ_CONFIG_HPP_
