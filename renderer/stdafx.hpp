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

namespace tgir {
typedef double Real;

typedef hi::basic_vector2<tgir::Real> Vector2;
typedef hi::basic_vector3<tgir::Real> Vector3;
typedef hi::basic_vector4<tgir::Real> Vector4;

typedef hi::basic_spectrum SpectrumData;
typedef SpectrumData::value_type Spectrum;
typedef hi::basic_vector3<Spectrum> SpectrumVector;
typedef hi::basic_vector3<tgir::Real> CIE_XYZ_Color;

tgir::Real const EPSILON = 1e-9;

//__declspec(thread) int value; // これでスレッド毎にデータが生成される

}  // end of namespace tgir

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
