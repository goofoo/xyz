#pragma once

#define NOMINMAX  ///< windows.h���� min(a,b) �� max(a,b) �}�N�������O���܂�.
#define VC_EXTRALEAN         ///< �������� NOservice ��`���`���܂�.
#define WIN32_LEAN_AND_MEAN  ///< windows.h����g�p����Ă��Ȃ����������O���܂�.
#define _USE_MATH_DEFINES    ///< cmath�̐��w�萔���g�p���܂�.

// C �����^�C�� �w�b�_�[ �t�@�C��
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cerrno>

#include <array>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include <vector>
#include <deque>
#include <queue>
#include <stack>
#include <map>
#include <set>
#include <algorithm>
#include <functional>
#include <limits>

// Windows �w�b�_�[ �t�@�C��:
#include <tchar.h>
#include <windows.h>
#include <process.h>
#include <shellapi.h>

#ifndef INOUT
#define INOUT
#endif

#define HI_TO_STRING(x) #x
#define HI_ARRAYS_LENGTH(ary) (sizeof(ary) / sizeof(ary[0]))

#if defined(NDEBUG)
#define HI_RANGE_CHECK(a, i)
#else
#define HI_RANGE_CHECK(a, i) \
  assert(((i) < (a).size()) && _TEXT(HI_TO_STRING(__LINE__)))
#endif

#define HI_DISALLOW_COPY_AND_ASSIGN(TypeName) \
  \
private:                                      \
  TypeName(TypeName const &);                 \
  TypeName &operator=(TypeName const &)

/// <summary>
/// ��{�^
/// </summary>
namespace hi {

typedef unsigned __int16 word;
typedef unsigned __int32 dword;
typedef unsigned __int64 qword;

typedef unsigned __int8 ubyte;
typedef unsigned __int16 ushort;
typedef unsigned __int32 uint;
typedef unsigned __int64 ulong;

typedef signed __int8 sbyte;
typedef signed __int16 sshort;
typedef signed __int32 sint;
typedef signed __int64 slong;

}  // end of namespace hi

#if 0
namespace hi
{
  /// <summary>
  /// nullptr
  /// </summary>
  class sealed
  {
  public:
    template <typename T> // ������^�̔񃁃��o�ւ̃|�C���^�ɕϊ��ł���
      operator T * () const
    { return 0; }

    template <class C, typename T> // ������^�̃����o�ւ̃|�C���^�ɕϊ��ł���
      operator T C::*() const 
    { return 0; }

  private:
    void operator & () const; // �A�h���X�͎擾�ł��Ȃ�
  }
  const nullptr;

} // end of namespace hi
#endif

namespace hi {
/// <summary>
/// ���Ԑ��_�E���L���X�g(Release���[�h�̂Ƃ��́C������static_cast�ɒu����������)
/// </summary>
template <typename T, typename S>
inline T polymorphic_downcast(S p) {
#if defined(NDEBUG)
  return static_cast<T>(p);
#else
  return dynamic_cast<T>(p);
#endif
}
}  // end of namespace hi

/// <summary>
/// �ȈՔŃ����_��
/// </summary>
namespace hi {
class __lambda {
 public:
  inline __lambda() {}

 private:
  __lambda(__lambda const &);
  __lambda &operator=(__lambda const &);
  void operator&() const;  // �A�h���X�͎擾�ł��Ȃ�
} const lambda;

#define HI_IMPLEMENT_LAMBDA(Op, Fn)                                         \
  template <typename T>                                                     \
  inline std::binder1st<Fn> operator Op(__lambda const &, T const &value) { \
    return std::binder1st(Fn(), value);                                     \
  }                                                                         \
  template <typename T>                                                     \
  inline std::binder2nd<Fn> operator Op(T const &value, __lambda const &) { \
    return std::binder2nd(Fn(), value);                                     \
  }
HI_IMPLEMENT_LAMBDA(==, std::equal_to<T>)
HI_IMPLEMENT_LAMBDA(!=, std::not_equal_to<T>)
HI_IMPLEMENT_LAMBDA(>, std::greater<T>)
HI_IMPLEMENT_LAMBDA(>=, std::greater_equal<T>)
HI_IMPLEMENT_LAMBDA(<, std::less<T>)
HI_IMPLEMENT_LAMBDA(<=, std::less_equal<T>)
HI_IMPLEMENT_LAMBDA(&&, std::logical_and<T>)
HI_IMPLEMENT_LAMBDA(||, std::logical_or<T>)
#undef HI_IMPLEMENT_LAMBDA
}

/// <summary>
/// �����񏈗��Q�̃��l�C��
/// </summary>
namespace std {
typedef _TCHAR tchar_t;

#if defined(_UNICODE)
typedef wstring tstring;
typedef wios tios;
typedef wstreambuf tstreambuf;
typedef wistream tistream;
typedef wostream tostream;
typedef wiostream tiostream;
typedef wstringbuf tstringbuf;
typedef wistringstream tistringstream;
typedef wostringstream tostringstream;
typedef wstringstream tstringstream;
typedef wfilebuf tfilebuf;
typedef wifstream tifstream;
typedef wofstream tofstream;
typedef wfstream tfstream;
#else
typedef string tstring;
typedef ios tios;
typedef streambuf tstreambuf;
typedef istream tistream;
typedef ostream tostream;
typedef iostream tiostream;
typedef stringbuf tstringbuf;
typedef istringstream tistringstream;
typedef ostringstream tostringstream;
typedef stringstream tstringstream;
typedef filebuf tfilebuf;
typedef ifstream tifstream;
typedef ofstream tofstream;
typedef fstream tfstream;
#endif
}

#if defined(__HI_LANG_MATERIALIZE_GLOBAL_OBJECT__)
#define __HI_LANG_GLOBAL_EXTERN__
#define __HI_LANG_GLOBAL_DEFAULT__(x) (x)
#else
#define __HI_LANG_GLOBAL_EXTERN__ extern
#define __HI_LANG_GLOBAL_DEFAULT__(x)
#endif

namespace std {
#if defined(_UNICODE)
__HI_LANG_GLOBAL_EXTERN__ tistream &tcin __HI_LANG_GLOBAL_DEFAULT__(wcin);
__HI_LANG_GLOBAL_EXTERN__ tostream &tcout __HI_LANG_GLOBAL_DEFAULT__(wcout);
__HI_LANG_GLOBAL_EXTERN__ tostream &tcerr __HI_LANG_GLOBAL_DEFAULT__(wcerr);
__HI_LANG_GLOBAL_EXTERN__ tostream &tclog __HI_LANG_GLOBAL_DEFAULT__(wclog);
#else
__HI_LANG_GLOBAL_EXTERN tistream &tcin __HI_LANG_GLOBAL_DEFAULT__(cin);
__HI_LANG_GLOBAL_EXTERN tostream &tcout __HI_LANG_GLOBAL_DEFAULT__(cout);
__HI_LANG_GLOBAL_EXTERN tostream &tcerr __HI_LANG_GLOBAL_DEFAULT__(cerr);
__HI_LANG_GLOBAL_EXTERN tostream &tclog __HI_LANG_GLOBAL_DEFAULT__(clog);
#endif
}  // end of namespace std
