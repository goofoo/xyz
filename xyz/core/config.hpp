#ifndef XYZ_CONFIG_HPP_
#define XYZ_CONFIG_HPP_

#include "../stdafx.hpp"

#define RENDERER_NAME _TEXT("XYZ")
#define RENDERER_BUILD _TEXT(__DATE__)

// __FILE__ ファイル名
// __TIMESTAMP__ ファイルの最終更新日時
// __TIME__ コンパイル時間
// __DATE__ コンパイル日
// __LINE__ 行番号
// __FUNCTION__ 関数名(C99の__func__)
// __cplusplus C++でコンパイルされている

#define MQODIR_PATH _TEXT("./")
#define RENDIR_PATH _TEXT("./")
#define OUTDIR_PATH _TEXT("./")

// 全てのアルゴリズムに共通するパラメータ
#define XYZ_CONFIG_kMaxPathLength 8    ///< 双方向経路追跡法で生成する最大経路長
#define XYZ_CONFIG_kSampleCount 40000  ///< １ループのサンプル数
#define XYZ_CONFIG_kSystemCount 16     ///< 個別のサンプラーの数

#define XYZ_CONFIG_kLargeStepProb \
  float_t(0.25)  ///< Simplified MLT でラージステップ変異を行う確率

// 経路追跡法に関する定数
#define XYZ_CONFIG_kMaxRandomWalkDepth 16  ///< ランダムウォークでの最大追跡深度
#define XYZ_CONFIG_kMaxBranchCount 10      ///< ランダムウォークでの最大分岐数
#define XYZ_CONFIG_kStratifiedSize 100     ///< Population Annealing での層化数
#define CONFIG_SEEDS            \
  (XYZ_CONFIG_kStratifiedSize * \
   XYZ_CONFIG_kStratifiedSize)  // 積分値の推定のためのサンプル数

#define CONFIG_BLOCK_SIZE 10  // 粒子フィルタのフィルム層化サイズ
#define CONFIG_LIGHT_SIZE 1   // 粒子フィルタの光源層化サイズ
#define CONFIG_COLOR_SIZE 1   // 粒子フィルタのスペクトル層化サイズ
#define CONFIG_FILTER_SIZE                                     \
  (CONFIG_BLOCK_SIZE * CONFIG_BLOCK_SIZE * CONFIG_LIGHT_SIZE * \
   CONFIG_LIGHT_SIZE * CONFIG_COLOR_SIZE)
#define CONFIG_SS_TYPE 3

// フォトンマップ用の定数
#define CONFIG_IMPORTON_STRATIFY 100  // インポートンの各次元の層別化数
#define CONFIG_IMPORTON_SIZE  \
  (CONFIG_IMPORTON_STRATIFY * \
   CONFIG_IMPORTON_STRATIFY)  // 射出するインポートン数

#define CONFIG_BUNDLE_SIZE 25  // 同時に追跡するフォトン数
#define CONFIG_N_SIZE (CONFIG_U_SIZE * CONFIG_V_SIZE)  //
#define CONFIG_M_SIZE (CONFIG_N_SIZE * CONFIG_BUNDLE_SIZE)
#define CONFIG_LIGHT_POWER 100000

#define CONFIG_SEARCH_RADIUS float_t(0.1)  // フォトンやインポートンの探索半径
#define CONFIG_EMITED_SIZE 100000          // フォトンの射出数

#define CONFIG_ALBEDO float_t(0.65)  // 拡散反射率
#define CONFIG_TERMINATE_THRESHOLD \
  float_t(0.2)  // 0.2/0.5 // ロシアンルーレットの閾値(sigma)

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
