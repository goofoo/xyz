// stdafx.hpp : �W���̃V�X�e�� �C���N���[�h �t�@�C���̃C���N���[�h
// �t�@�C���A�܂���
// �Q�Ɖ񐔂������A�����܂�ύX����Ȃ��A�v���W�F�N�g��p�̃C���N���[�h
// �t�@�C��
// ���L�q���܂��B
//

#pragma once

#pragma warning(disable : 4505)

#include <hi/lang.hpp>
#include <hi/thread.hpp>
#include <hi/math.hpp>
#include <hi/image.hpp>
#include <hi/file.hpp>
#include <hi/util.hpp>
#include <hi/tgl.hpp>
#include <hi/sds.hpp>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

namespace xyz {
typedef float float_t;

typedef hi::basic_vector2<float_t> float2_t;
typedef hi::basic_vector3<float_t> float3_t;
typedef hi::basic_vector4<float_t> float4_t;
typedef hi::basic_vector3<float_t> CIE_XYZ_Color;

float_t const EPSILON = float_t(1e-6);
// float_t const EPSILON = float_t(1e-15);

//__declspec(thread) int value; // ����ŃX���b�h���Ƀf�[�^�����������

typedef std::pair<std::size_t, CIE_XYZ_Color> pixel_descriptor_t;

// �[���t���X�l����: fake fresnel term
inline float_t FakeFresnelTerm(__in float_t const hk,
                               __in float_t const fresnel0 = 0.9) {
  return fresnel0 + (1 - fresnel0) * hi::fifth_power_of(1 - hk);
}

// get random vector with importance sampleing (p(theta,phi) = cos(theta)/pi)
inline float_t ComputeDiffusedVector(__in float_t const fMuTheta,
                                     __in float_t const fMuPhi,
                                     __in float3_t const &vTangent,
                                     __in float3_t const &vNormal,
                                     __in float3_t const &vBinormal,
                                     __out float3_t *const pvDiffusedVector) {
  float_t const fPhi = float_t(2 * M_PI) * fMuPhi;
  float_t const fCosTheta = std::sqrt(fMuTheta);
  float_t const fSinTheta = std::sqrt(1 - fMuTheta);
  *pvDiffusedVector = vTangent * (fSinTheta * std::cos(fPhi)) +
                      vBinormal * (fSinTheta * std::sin(fPhi)) +
                      vNormal * (fCosTheta);
  return fCosTheta;
}

inline float_t ComputeDiffusedVector(__in float3_t const &vMu,
                                     __in float3_t const &vTangent,
                                     __in float3_t const &vNormal,
                                     __in float3_t const &vBinormal,
                                     __out float3_t *const pvDiffusedVector) {
  return ComputeDiffusedVector(vMu[0], vMu[1], vTangent, vNormal, vBinormal,
                               pvDiffusedVector);
}

/// <summary>
/// �@���x�N�g���E�@���Ɠ��˕����x�N�g���̓��ρE���˕����x�N�g�����狾�ʔ��˃x�N�g�����v�Z����D
/// </summary>
/// <param name="vNormal">�@���x�N�g��</param>
/// <param name="fCosTheta">�@���Ɠ��˕����x�N�g���̓���</param>
/// <param
/// name="vReflectedVector">���͂Ƃ��Ă͓��˕����̃x�N�g��/�o�͂Ƃ��Ă͋��ʔ��˃x�N�g��</param>
/// <returns></returns>
/// <remarks>
/// �����Ɋ�Â��ċ��ʔ��˃x�N�g�����v�Z����D
/// <newpara>
/// vReflectedVector += 2 * fCosTheta * vNormalVector
/// </newpara>
/// </remarks>
inline void ComputeReflectedVector(__in float3_t const &vNormalVector,
                                   __in float_t const fCosTheta,
                                   __inout float3_t *const vReflectedVector) {
  (*vReflectedVector) += (float_t(2) * fCosTheta) * vNormalVector;
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
inline void ComputeRefractedVector(__in float3_t const &vNormalVector,
                                   __in float_t const gmc,
                                   __in float_t const rIndexOfRefraction,
                                   __inout float3_t *const pvRefractedVector) {
  *pvRefractedVector =
      (*pvRefractedVector - gmc * vNormalVector) / rIndexOfRefraction;
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
inline float_t ComputeGlossyHalfwayVector(
    __in float_t const rMuTheta, __in float_t const rMuPhi,
    __in float_t const rShininess, __in float3_t const &vTangent,
    __in float3_t const &vNormal, __in float3_t const &vBinormal,
    __out float3_t *const pvHalfwayVector) {
  float_t const rPhi = float_t(2 * M_PI) * rMuPhi;
  float_t const fCosTheta = std::pow(rMuTheta, hi::rcp(rShininess));
  float_t const rSinTheta = std::sqrt(float_t(1) - hi::square_of(fCosTheta));

  (*pvHalfwayVector) = vTangent * (rSinTheta * std::cos(rPhi)) +
                       vBinormal * (rSinTheta * std::sin(rPhi)) +
                       vNormal * (fCosTheta);

  return fCosTheta;
}

template <typename T>
inline T SafeDivide(T const &x, T const &y) {
  return (y < x) ? (y / x) : ((x > T(0)) ? T(1) : T(0));
}

inline bool Contains(__in float3_t const &pnt, __in float3_t const &min,
                     __in float3_t const &max) {
  return (min[0] <= pnt[0]) && (pnt[0] <= max[0]) && (min[1] <= pnt[1]) &&
         (pnt[1] <= max[1]) && (min[2] <= pnt[2]) && (pnt[2] <= max[2]);
}

inline bool Overlap(__in float3_t const &a_min, __in float3_t const &a_max,
                    __in float3_t const &b_min, __in float3_t const &b_max) {
  return (a_min[0] <= b_max[0]) && (a_max[0] > b_min[0]) &&
         (a_min[1] <= b_max[1]) && (a_max[1] > b_min[1]) &&
         (a_min[2] <= b_max[2]) && (a_max[2] > b_min[2]);
}

inline float_t SurfaceArea(__in float3_t const &box_size) {
  return box_size[0] * box_size[1] + box_size[1] * box_size[2] +
         box_size[2] * box_size[0];
}

}  // end of namespace xyz

// �������[���[�N�̃`�F�b�N���s��
// main�̍ŏ���
// ::_CrtSetDbgFlag(_CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF);
// ���ĂԂ���
#ifndef NDEBUG
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#define new new (_NORMAL_BLOCK, __FILE__, __LINE__)
#endif

/*
MFC�𗘗p����ꍇ

#if _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#define new new(__FILE__, __LINE__) //< (��)
#endif
*/
