// stdafx.hpp : 標準のシステム インクルード ファイルのインクルード
// ファイル、または
// 参照回数が多く、かつあまり変更されない、プロジェクト専用のインクルード
// ファイル
// を記述します。
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

//__declspec(thread) int value; // これでスレッド毎にデータが生成される

typedef std::pair<std::size_t, CIE_XYZ_Color> pixel_descriptor_t;

// 擬似フレスネル項: fake fresnel term
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
/// 法線ベクトル・法線と入射方向ベクトルの内積・入射方向ベクトルから鏡面反射ベクトルを計算する．
/// </summary>
/// <param name="vNormal">法線ベクトル</param>
/// <param name="fCosTheta">法線と入射方向ベクトルの内積</param>
/// <param
/// name="vReflectedVector">入力としては入射方向のベクトル/出力としては鏡面反射ベクトル</param>
/// <returns></returns>
/// <remarks>
/// 下式に基づいて鏡面反射ベクトルを計算する．
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
/// 相対屈折率に基づいた屈折ベクトルを計算する．
/// </summary>
/// <param name="vNormalVector">法線ベクトル</param>
/// <param name="gcm">詳細は後述</param>
/// <param
/// name="rIndexOfRefraction">相対屈折率(入射後の媒質の屈折率/入射前の媒質の屈折率)</param>
/// <param
/// name="vRefractedVector">入力としては入射方向のベクトル/出力としては相対屈折率に基づいた屈折ベクトル</param>
/// <returns></returns>
/// <remarks>
/// 下式に基づいて屈折ベクトルを計算する．
/// <newpara>
/// vRefractedVector = (vRefractedVector - gmc * vNormalVector) /
/// rIndexOfRefraction
/// </newpara>
/// ここで，gmcの定義は下式の通りである．
/// <newpara>
/// gmc = c - sqrt(rIndexOfRefraction^{2} + c^{2} - 1)
/// </newpara>
/// ただし，
/// <newpara>
/// c = dot(vRefractedVector, vNormal)
/// </newpara>
/// である．
/// </remarks>
inline void ComputeRefractedVector(__in float3_t const &vNormalVector,
                                   __in float_t const gmc,
                                   __in float_t const rIndexOfRefraction,
                                   __inout float3_t *const pvRefractedVector) {
  *pvRefractedVector =
      (*pvRefractedVector - gmc * vNormalVector) / rIndexOfRefraction;
}

/// <summary>
/// 確率分布ph(vHalfwayVector)に基づいて中間(halfway)ベクトルを計算する．
/// </summary>
/// <param
/// name="rMuTheta">結果のベクトルの緯度方向を決定するための乱数値</param>
/// <param name="rMuPhi">結果のベクトルの経度方向を決定するための乱数値</param>
/// <param
/// name="rShininess">鏡面反射ローブの形状を決定する変数(大きいほど鋭い反射になる)</param>
/// <param name="vTangent">接ベクトル(基底ベクトルのX軸に相当する)</param>
/// <param name="vNormal">法線ベクトル(基底ベクトルのY軸に相当する)</param>
/// <param name="vBinormal">従法線ベクトル(基底ベクトルのZ軸に相当する)</param>
/// <param
/// name="vHalfwayVector">入力としては入射方向のベクトル/出力としては確率分布ph(vHalfwayVector)に基づいた中間ベクトル</param>
/// <returns></returns>
/// <remarks>
/// 確率分布ph(vHalfwayVector)に基づいて中間ベクトルを計算する．
/// <newpara>
/// ph(vHalfwayVector) = (rShininess + 1) / 2pi * dot(vNormal,
/// vHalfwayVector)^{rShininess}
/// </newpara>
/// この中間ベクトルを使って生成される反射ベクトルは，確率分布p(vGlossyVector)に従う．
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

// メモリーリークのチェックを行う
// mainの最初で
// ::_CrtSetDbgFlag(_CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF);
// を呼ぶこと
#ifndef NDEBUG
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#define new new (_NORMAL_BLOCK, __FILE__, __LINE__)
#endif

/*
MFCを利用する場合

#if _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#define new new(__FILE__, __LINE__) //< (※)
#endif
*/
