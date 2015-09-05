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

//__declspec(thread) int value; // ����ŃX���b�h���Ƀf�[�^�����������

}  // end of namespace tgir

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
