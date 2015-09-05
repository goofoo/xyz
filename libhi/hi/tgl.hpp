#pragma once

#include <hi/lang.hpp>
#include <GL/glut.h>

namespace hi {
/// <summary>
/// glColor
/// </summary>
template <typename T>
inline void glColor(T const, T const, T const);

template <typename T>
inline void glColor(GLbyte const red, GLbyte const green, GLbyte const blue) {
  ::glColor3b(red, green, blue);
}
template <typename T>
inline void glColor(GLdouble const red, GLdouble const green,
                    GLdouble const blue) {
  ::glColor3d(red, green, blue);
}
template <typename T>
inline void glColor(GLfloat const red, GLfloat const green,
                    GLfloat const blue) {
  ::glColor3f(red, green, blue);
}
template <typename T>
inline void glColor(GLint const red, GLint const green, GLint const blue) {
  ::glColor3i(red, green, blue);
}
template <typename T>
inline void glColor(GLshort const red, GLshort const green,
                    GLshort const blue) {
  ::glColor3s(red, green, blue);
}
template <typename T>
inline void glColor(GLubyte const red, GLubyte const green,
                    GLubyte const blue) {
  ::glColor3ub(red, green, blue);
}
template <typename T>
inline void glColor(GLuint const red, GLuint const green, GLuint const blue) {
  ::glColor3ui(red, green, blue);
}
template <typename T>
inline void glColor(GLushort const red, GLushort const green,
                    GLushort const blue) {
  ::glColor3us(red, green, blue);
}

template <typename T>
inline void glColor(T const, T const, T const, T const);

template <typename T>
inline void glColor(GLbyte const red, GLbyte const green, GLbyte const blue,
                    GLbyte const alpha) {
  ::glColor4b(red, green, blue, alpha);
}
template <typename T>
inline void glColor(GLdouble const red, GLdouble const green,
                    GLdouble const blue, GLdouble const alpha) {
  ::glColor4d(red, green, blue, alpha);
}
template <typename T>
inline void glColor(GLfloat const red, GLfloat const green, GLfloat const blue,
                    GLfloat const alpha) {
  ::glColor4f(red, green, blue, alpha);
}
template <typename T>
inline void glColor(GLint const red, GLint const green, GLint const blue,
                    GLint const alpha) {
  ::glColor4i(red, green, blue, alpha);
}
template <typename T>
inline void glColor(GLshort const red, GLshort const green, GLshort const blue,
                    GLshort const alpha) {
  ::glColor4s(red, green, blue, alpha);
}
template <typename T>
inline void glColor(GLubyte const red, GLubyte const green, GLubyte const blue,
                    GLubyte const alpha) {
  ::glColor4ub(red, green, blue, alpha);
}
template <typename T>
inline void glColor(GLuint const red, GLuint const green, GLuint const blue,
                    GLuint const alpha) {
  ::glColor4ui(red, green, blue, alpha);
}
template <typename T>
inline void glColor(GLushort const red, GLushort const green,
                    GLushort const blue, GLushort const alpha) {
  ::glColor4us(red, green, blue, alpha);
}

/// <summary>
/// glColorN
/// </summary>
template <std::size_t Sz, typename T>
inline void glColor(T const *);

template <>
inline void glColor<3, GLbyte>(GLbyte const *const v) {
  ::glColor3bv(v);
}
template <>
inline void glColor<3, GLdouble>(GLdouble const *const v) {
  ::glColor3dv(v);
}
template <>
inline void glColor<3, GLfloat>(GLfloat const *const v) {
  ::glColor3fv(v);
}
template <>
inline void glColor<3, GLint>(GLint const *const v) {
  ::glColor3iv(v);
}
template <>
inline void glColor<3, GLshort>(GLshort const *const v) {
  ::glColor3sv(v);
}
template <>
inline void glColor<3, GLubyte>(GLubyte const *const v) {
  ::glColor3ubv(v);
}
template <>
inline void glColor<3, GLuint>(GLuint const *const v) {
  ::glColor3uiv(v);
}
template <>
inline void glColor<3, GLushort>(GLushort const *const v) {
  ::glColor3usv(v);
}

template <>
inline void glColor<4, GLbyte>(GLbyte const *const v) {
  ::glColor4bv(v);
}
template <>
inline void glColor<4, GLdouble>(GLdouble const *const v) {
  ::glColor4dv(v);
}
template <>
inline void glColor<4, GLfloat>(GLfloat const *const v) {
  ::glColor4fv(v);
}
template <>
inline void glColor<4, GLint>(GLint const *const v) {
  ::glColor4iv(v);
}
template <>
inline void glColor<4, GLshort>(GLshort const *const v) {
  ::glColor4sv(v);
}
template <>
inline void glColor<4, GLubyte>(GLubyte const *const v) {
  ::glColor4ubv(v);
}
template <>
inline void glColor<4, GLuint>(GLuint const *const v) {
  ::glColor4uiv(v);
}
template <>
inline void glColor<4, GLushort>(GLushort const *const v) {
  ::glColor4usv(v);
}

/// <summary>
/// glLight
/// </summary>
template <typename T>
inline void glLight(GLenum const, GLenum const, T const);

template <>
inline void glLight(GLenum const light, GLenum const pname,
                    GLfloat const param) {
  ::glLightf(light, pname, param);
}
template <>
inline void glLight(GLenum const light, GLenum const pname,
                    GLfloat *const params) {
  ::glLightfv(light, pname, params);
}
template <>
inline void glLight(GLenum const light, GLenum const pname,
                    GLfloat const *const params) {
  ::glLightfv(light, pname, params);
}
template <>
inline void glLight(GLenum const light, GLenum const pname, GLint const param) {
  ::glLighti(light, pname, param);
}
template <>
inline void glLight(GLenum const light, GLenum const pname,
                    GLint *const params) {
  ::glLightiv(light, pname, params);
}
template <>
inline void glLight(GLenum const light, GLenum const pname,
                    GLint const *const params) {
  ::glLightiv(light, pname, params);
}

/// <summary>
/// glMaterial
/// </summary>
template <typename T>
inline void glMaterial(GLenum const face, GLenum const pname, T const);

template <>
inline void glMaterial(GLenum face, GLenum pname, GLfloat const param) {
  ::glMaterialf(face, pname, param);
}
template <>
inline void glMaterial(GLenum face, GLenum pname, GLfloat *const params) {
  ::glMaterialfv(face, pname, params);
}
template <>
inline void glMaterial(GLenum face, GLenum pname, GLfloat const *const params) {
  ::glMaterialfv(face, pname, params);
}
template <>
inline void glMaterial(GLenum face, GLenum pname, GLint const param) {
  ::glMateriali(face, pname, param);
}
template <>
inline void glMaterial(GLenum face, GLenum pname, GLint *const params) {
  ::glMaterialiv(face, pname, params);
}
template <>
inline void glMaterial(GLenum face, GLenum pname, GLint const *const params) {
  ::glMaterialiv(face, pname, params);
}

/// <summary>
/// glNormal
/// </summary>
template <typename T>
inline void glNormal(T const, T const, T const);

template <>
inline void glNormal(GLbyte const nx, GLbyte const ny, GLbyte const nz) {
  ::glNormal3b(nx, ny, nz);
}
template <>
inline void glNormal(GLdouble const nx, GLdouble const ny, GLdouble const nz) {
  ::glNormal3d(nx, ny, nz);
}
template <>
inline void glNormal(GLfloat const nx, GLfloat const ny, GLfloat const nz) {
  ::glNormal3f(nx, ny, nz);
}
template <>
inline void glNormal(GLint const nx, GLint const ny, GLint const nz) {
  ::glNormal3i(nx, ny, nz);
}
template <>
inline void glNormal(GLshort const nx, GLshort const ny, GLshort const nz) {
  ::glNormal3s(nx, ny, nz);
}

template <typename T>
inline void glNormal(T const *);

template <>
inline void glNormal(GLbyte const *const v) {
  ::glNormal3bv(v);
}
template <>
inline void glNormal(GLdouble const *const v) {
  ::glNormal3dv(v);
}
template <>
inline void glNormal(GLfloat const *const v) {
  ::glNormal3fv(v);
}
template <>
inline void glNormal(GLint const *const v) {
  ::glNormal3iv(v);
}
template <>
inline void glNormal(GLshort const *const v) {
  ::glNormal3sv(v);
}

/// <summary>
/// glTranslate
/// </summary>
template <typename T>
inline void glTranslate(T const, T const, T const);

template <>
inline void glTranslate(GLdouble const x, GLdouble const y, GLdouble const z) {
  ::glTranslated(x, y, z);
}
template <>
inline void glTranslate(GLfloat const x, GLfloat const y, GLfloat const z) {
  ::glTranslatef(x, y, z);
}
template <>
inline void glTranslate(GLint const x, GLint const y, GLint const z) {
  ::glTranslated(x, y, z);
}

/// <summary>
/// glVertex2d
/// </summary>
template <typename T>
inline void glVertex(T const, T const);

template <>
inline void glVertex(GLdouble const x, GLdouble const y) {
  ::glVertex2d(x, y);
}
template <>
inline void glVertex(GLfloat const x, GLfloat const y) {
  ::glVertex2f(x, y);
}
template <>
inline void glVertex(GLint const x, GLint const y) {
  ::glVertex2i(x, y);
}
template <>
inline void glVertex(GLshort const x, GLshort const y) {
  ::glVertex2s(x, y);
}

/// <summary>
/// glVertex3d
/// </summary>
template <typename T>
inline void glVertex(T const, T const, T const);

template <>
inline void glVertex(GLdouble const x, GLdouble const y, GLdouble const z) {
  ::glVertex3d(x, y, z);
}
template <>
inline void glVertex(GLfloat const x, GLfloat const y, GLfloat const z) {
  ::glVertex3f(x, y, z);
}
template <>
inline void glVertex(GLint const x, GLint const y, GLint const z) {
  ::glVertex3i(x, y, z);
}
template <>
inline void glVertex(GLshort const x, GLshort const y, GLshort const z) {
  ::glVertex3s(x, y, z);
}

/// <summary>
/// glVertex4d
/// </summary>
template <typename T>
inline void glVertex(T const, T const, T const, T const);

template <>
inline void glVertex(GLdouble const x, GLdouble const y, GLdouble const z,
                     GLdouble const w) {
  ::glVertex4d(x, y, z, w);
}
template <>
inline void glVertex(GLfloat const x, GLfloat const y, GLfloat const z,
                     GLfloat const w) {
  ::glVertex4f(x, y, z, w);
}
template <>
inline void glVertex(GLint const x, GLint const y, GLint const z,
                     GLint const w) {
  ::glVertex4i(x, y, z, w);
}
template <>
inline void glVertex(GLshort const x, GLshort const y, GLshort const z,
                     GLshort const w) {
  ::glVertex4s(x, y, z, w);
}

/// <summary>
/// glVertexNd
/// </summary>
template <std::size_t Sz, typename T>
inline void glVertex(T const *);

template <>
inline void glVertex<2, GLdouble>(GLdouble const *const v) {
  ::glVertex2dv(v);
}
template <>
inline void glVertex<2, GLfloat>(GLfloat const *const v) {
  ::glVertex2fv(v);
}
template <>
inline void glVertex<2, GLint>(GLint const *const v) {
  ::glVertex2iv(v);
}
template <>
inline void glVertex<2, GLshort>(GLshort const *const v) {
  ::glVertex2sv(v);
}

template <>
inline void glVertex<3, GLdouble>(GLdouble const *const v) {
  ::glVertex3dv(v);
}
template <>
inline void glVertex<3, GLfloat>(GLfloat const *const v) {
  ::glVertex3fv(v);
}
template <>
inline void glVertex<3, GLint>(GLint const *const v) {
  ::glVertex3iv(v);
}
template <>
inline void glVertex<3, GLshort>(GLshort const *const v) {
  ::glVertex3sv(v);
}

template <>
inline void glVertex<4, GLdouble>(GLdouble const *const v) {
  ::glVertex4dv(v);
}
template <>
inline void glVertex<4, GLfloat>(GLfloat const *const v) {
  ::glVertex4fv(v);
}
template <>
inline void glVertex<4, GLint>(GLint const *const v) {
  ::glVertex4iv(v);
}
template <>
inline void glVertex<4, GLshort>(GLshort const *const v) {
  ::glVertex4sv(v);
}

}  // end of namespace tgl
