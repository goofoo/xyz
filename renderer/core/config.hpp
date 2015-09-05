#ifndef __TGIR_CORE_CONFIG_HPP__
#define __TGIR_CORE_CONFIG_HPP__

#define RENDERER_NAME _TEXT("Tiny Global Illumination Renderer")
#define RENDERER_BUILD _TEXT(__DATE__)

// __FILE__ ファイル名
// __TIMESTAMP__ ファイルの最終更新日時
// __TIME__ コンパイル時間
// __DATE__ コンパイル日
// __LINE__ 行番号
// __FUNCTION__ 関数名(C99の__func__)
// __cplusplus C++でコンパイルされている

#if 0
#define MQODIR_PATH _TEXT("C:/home/res/mqo/")
#define RENDIR_PATH _TEXT("C:/home/res/ren/")
#define OUTDIR_PATH _TEXT("C:/home/res/out/")
#else
#define MQODIR_PATH _TEXT("./")
#define RENDIR_PATH _TEXT("./")
#define OUTDIR_PATH _TEXT("./")
#endif

// 全てのアルゴリズムに共通するパラメータ
#define TGIR_CONFIG_kMaxPathLength 8  ///< 双方向経路追跡法で生成する最大経路長
#define TGIR_CONFIG_kSampleCount 40000  ///< １ループのサンプル数
#define TGIR_CONFIG_kSystemCount 16     ///< 個別のサンプラーの数

#define TGIR_CONFIG_kLargeStepProb \
  0.25  ///< Simplified MLT でラージステップ変異を行う確率

// 経路追跡法に関する定数
#define TGIR_CONFIG_kMaxRandomWalkDepth 8  ///< ランダムウォークでの最大追跡深度
#define TGIR_CONFIG_kMaxBranchCount 10     ///< ランダムウォークでの最大分岐数
#define TGIR_CONFIG_kStratifiedSize 100    ///< Population Annealing での層化数
#define CONFIG_SEEDS             \
  (TGIR_CONFIG_kStratifiedSize * \
   TGIR_CONFIG_kStratifiedSize)  // 積分値の推定のためのサンプル数

#define CONFIG_BLOCK_SIZE 4  // 粒子フィルタのフィルム層化サイズ
#define CONFIG_LIGHT_SIZE 1  // 粒子フィルタの光源層化サイズ
#define CONFIG_COLOR_SIZE 1  // 粒子フィルタのスペクトル層化サイズ
#define CONFIG_FILTER_SIZE                                     \
  (CONFIG_BLOCK_SIZE * CONFIG_BLOCK_SIZE * CONFIG_LIGHT_SIZE * \
   CONFIG_LIGHT_SIZE * CONFIG_COLOR_SIZE)
#define CONFIG_SS_TYPE 3

// フォトンマップ用の定数
#define CONFIG_IMPORTON_STRATIFY 1024  // インポートンの各次元の層別化数
#define CONFIG_IMPORTON_SIZE  \
  (CONFIG_IMPORTON_STRATIFY * \
   CONFIG_IMPORTON_STRATIFY)   // 射出するインポートン数
#define CONFIG_U_SIZE 4        //
#define CONFIG_V_SIZE 16       //
#define CONFIG_BUNDLE_SIZE 25  // 同時に追跡するフォトン数
#define CONFIG_N_SIZE (CONFIG_U_SIZE * CONFIG_V_SIZE)  //
#define CONFIG_M_SIZE (CONFIG_N_SIZE * CONFIG_BUNDLE_SIZE)
#define CONFIG_LIGHT_POWER 100000
#define CONFIG_SEARCH_RADIUS 0.2  // フォトンやインポートンの探索半径

#define CONFIG_ALBEDO 0.65        // 拡散反射率
#define CONFIG_EMITED_SIZE 10000  // フォトンの射出数
#define CONFIG_TERMINATE_THRESHOLD \
  0.2  // 0.2/0.5 // ロシアンルーレットの閾値(sigma)

//#define CONFIG_PINHOLE_CAMERA_MODEL // ピンホールカメラかどうかのフラグ
//#define CONFIG_DOUBLE_LIGHT // 光源が両面有効かどうかのフラグ
//#define CONFIG_BACKFACE_CULLING // 背面カリングを行うかのフラグ
//#define CONFIG_RUSSIAN_RULLET //
//経路追跡のロシアンルーレットを行うかどうかのフラグ
#define CONFIG_MULTIPLE_INTEGRAION
#define CONFIG_CONCURRENT_LOCK
#define CONFIG_SPACE_SUBDIVISION

//#define CONFIG_HDR_STRICT

// 特定の経路をレンダリングするかどうかのフラグ
#define CONFIG_EXPLICIT_GLOSSY  // 光沢反射で明示的な直接照明計算を行うかどうか
#define CONFIG_RENDER_CAUSTICS  // 経路追跡法における集光模様経路

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

// 擬似フレスネル項: fake fresnel term
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
/// 法線ベクトル・法線と入射方向ベクトルの内積・入射方向ベクトルから鏡面反射ベクトルを計算する．
/// </summary>
/// <param name="vNormal">法線ベクトル</param>
/// <param name="rCosTheta">法線と入射方向ベクトルの内積</param>
/// <param
/// name="vReflectedVector">入力としては入射方向のベクトル/出力としては鏡面反射ベクトル</param>
/// <returns></returns>
/// <remarks>
/// 下式に基づいて鏡面反射ベクトルを計算する．
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
inline void ComputeRefractedVector(__in tgir::Vector3 const &vNormalVector,
                                   __in tgir::Real const gmc,
                                   __in tgir::Real const rIndexOfRefraction,
                                   __inout tgir::Vector3 &vRefractedVector) {
  vRefractedVector =
      (vRefractedVector - gmc * vNormalVector) / rIndexOfRefraction;
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
