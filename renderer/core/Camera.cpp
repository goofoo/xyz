#include "Camera.hpp"
#include "sampler/PrimarySample.hpp"

namespace {
inline void SampleDisk(
    __in tgir::Real const &s,                   ///< [in] パラメータ1
    __in tgir::Real const &t,                   ///< [in] パラメータ2
    __in tgir::Real const &fDiskRadius,         ///< [in] ディスクの半径
    __out tgir::Real *const pfLensCoordinateX,  ///< [out] 結果のx座標
    __out tgir::Real *const pfLensCoordinateY)  ///< [out] 結果のy座標
{
  tgir::Real const fRadius = fDiskRadius * std::sqrt(s);
  tgir::Real const fAngle = (2 * M_PI) * t;
  *pfLensCoordinateX = fRadius * std::cos(fAngle);
  *pfLensCoordinateY = fRadius * std::sin(fAngle);
}

}  // end of unnamed namespace

namespace tgir {
Camera::Camera()
    : vPrincipalPoint_(0.00, 1.00, 4.00),
      vGazingPoint_(0.00, 1.00, 0.00),
      fFilmSize_(0.025),
      fAspectRatio_(1.000),
      fFocalLength_(0.035),
      fFNumber_(1.4) {
  Reset();
}

void Camera::Render() const {
  tgir::Real const fFieldOfViewY =
      hi::to_degrees(2 * std::atan2(fFilmSize_ / 2, fLensFilmDistance_));

  ::glMatrixMode(GL_PROJECTION);
  ::glLoadIdentity();
  ::gluPerspective(fFieldOfViewY, fAspectRatio_, 0.1, 100.0);

  ::glMatrixMode(GL_MODELVIEW);
  ::glLoadIdentity();
  ::gluLookAt(vPrincipalPoint_[0], vPrincipalPoint_[1], vPrincipalPoint_[2],
              vGazingPoint_[0], vGazingPoint_[1], vGazingPoint_[2], 0, 1, 0);

  ::glPushMatrix();
  {
    hi::basic_vector4<GLfloat> ambient_and_diffuse(0, 0, 1, 1);
    hi::glMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,
                   ambient_and_diffuse.data());
    hi::glTranslate(vGazingPoint_[0], vGazingPoint_[1], vGazingPoint_[2]);
    ::glutSolidSphere(fLensRadius_, 20, 20);
  }
  ::glPopMatrix();
}

void Camera::MoveTo(tgir::Vector3 const &v) {
  vPrincipalPoint_ = v;
  Reset();
}
void Camera::MoveBy(tgir::Vector3 const &v) {
  vPrincipalPoint_ += v;
  Reset();
}
void Camera::LookAt(tgir::Vector3 const &v) {
  vGazingPoint_ = v;
  Reset();
}
void Camera::LookBy(tgir::Vector3 const &v) {
  vGazingPoint_ += v;
  Reset();
}

void Camera::Reset() {
  tgir::Vector3 const vViewDirection = vPrincipalPoint_ - vGazingPoint_;

  fFocusDistance_ = hi::length(vViewDirection);  // L
  fLensFilmDistance_ = fFocusDistance_ * fFocalLength_ /
                       (fFocusDistance_ - fFocalLength_);  // V=Lf(L-f)
  fLensRadius_ = fFocalLength_ / (2 * fFNumber_);          // R=f/(2F)

  ::_ftprintf_s(stderr, _TEXT("  ピントが合う面までの距離: %gm\n"),
                fFocusDistance_);
  ::_ftprintf_s(stderr, _TEXT("  レンズ／フィルム間の距離: %gmm\n"),
                fLensFilmDistance_ * 1000);
  ::_ftprintf_s(stderr, _TEXT("  レンズの有効半径(mm): %g\n"),
                fLensRadius_ * 1000);

  vBackwardVector_ = hi::normalize(vViewDirection);
  vRightwardVector_ =
      hi::normalize(hi::cross(tgir::Vector3(0, 1, 0), vBackwardVector_));
  vUpwardVector_ = hi::cross(vBackwardVector_, vRightwardVector_);
}

//============================================================================
//
// Pinhole camera model.
//
//============================================================================

void Camera::GetPrimaryRayOriginAndCameraBasis(
    __out tgir::Vector3
        *const pvRayOrigin,  ///< [out] 光線の始点(三次元空間内のレンズ上の位置)
    __out tgir::Vector3
        *const pvRightwardVector,  ///< [out] カメラの右方向ベクトル
    __out tgir::Vector3
        *const pvUpwardVector,  ///< [out] カメラの上方向ベクトル
    __out tgir::Vector3
        *const pvBackwardVector  ///< [out] カメラの後方向ベクトル
    ) const {
  *pvRayOrigin = vPrincipalPoint_;
  *pvRightwardVector = vRightwardVector_;
  *pvUpwardVector = vUpwardVector_;
  *pvBackwardVector = vBackwardVector_;
}

void Camera::GetPrimaryRayDirection(
    __in tgir::Real const &s,  ///< [in] フィルム面上の位置(X)．値域は[0,1]
    __in tgir::Real const &t,  ///< [in] フィルム面上の位置(Y)．値域は[0,1]
    __out tgir::Vector3 *const pvRayDirection  ///< [out] 一次光線の方向
    ) const {
  // フィルム位置から方向を決める
  tgir::Vector3 const vLocalDirection((s - 0.5) * fFilmSize_ * fAspectRatio_,
                                      (t - 0.5) * fFilmSize_,
                                      -fLensFilmDistance_);

  // 世界座標系に変換した後，正規化する
  *pvRayDirection = hi::normalize(vLocalDirection[0] * vRightwardVector_ +
                                  vLocalDirection[1] * vUpwardVector_ +
                                  vLocalDirection[2] * vBackwardVector_);
}

bool Camera::GetFilmPosition(
    __in tgir::Vector3 const &vDirection,  ///< [in]
                                           ///レンズ上の座標から三次元空間内の物体表面上の点へのベクトル
    __out tgir::Vector2
        *const pvFilmPosition  ///< [out] フィルム上の座標(X∈[0,1],Y∈[0,1])
    ) const {
  // 点の位置をローカル座標系に変換する
  tgir::Vector3 const vLocalDirection(hi::dot(vDirection, vRightwardVector_),
                                      hi::dot(vDirection, vUpwardVector_),
                                      hi::dot(vDirection, vBackwardVector_));

  assert(vLocalDirection[2] < 0);

  // フィルム面上に射影する(合同変換・フィルムサイズで値域を[0,1]に変換)
  tgir::Real const fConstantFactor =
      (fLensFilmDistance_ / -vLocalDirection[2]) / fFilmSize_;
  (*pvFilmPosition)[0] =
      tgir::Real(0.5) + vLocalDirection[0] * fConstantFactor / fAspectRatio_;
  (*pvFilmPosition)[1] = tgir::Real(0.5) + vLocalDirection[1] * fConstantFactor;

  // フィルム領域内かどうか判定して結果を返す
  return (0 <= (*pvFilmPosition)[0]) && ((*pvFilmPosition)[0] < 1) &&
         (0 <= (*pvFilmPosition)[1]) && ((*pvFilmPosition)[1] < 1);
}

//============================================================================
//
// Thin-lens camera model.
//
//============================================================================

void Camera::GetPrimaryRayOriginAndCameraBasis(
    __in tgir::Real const &s,                     ///< [in]
                                                  ///レンズ上の座標を計算するためのパラメータ(動径用)．値域は[0,1]
    __in tgir::Real const &t,                     ///< [in]
                                                  ///レンズ上の座標を計算するためのパラメータ(角度用)．値域は[0,1]
    __out tgir::Vector2 *const pvLensCoordinate,  ///< [out] レンズ上の座標(X,Y)
    __out tgir::Vector3
        *const pvRayOrigin,  ///< [out] 光線の始点(三次元空間内のレンズ上の位置)
    __out tgir::Vector3
        *const pvRightwardVector,  ///< [out] カメラの右方向ベクトル
    __out tgir::Vector3
        *const pvUpwardVector,  ///< [out] カメラの上方向ベクトル
    __out tgir::Vector3
        *const pvBackwardVector  ///< [out] カメラの後方向ベクトル
    ) const {
  // レンズ位置をサンプリング
  ::SampleDisk(s, t, fLensRadius_, &(*pvLensCoordinate)[0],
               &(*pvLensCoordinate)[1]);

  // 三次元空間内の座標に変換
  *pvRayOrigin = vPrincipalPoint_ + (*pvLensCoordinate)[0] * vRightwardVector_ +
                 (*pvLensCoordinate)[1] * vUpwardVector_;

  // 基底をセット
  *pvRightwardVector = vRightwardVector_;
  *pvUpwardVector = vUpwardVector_;
  *pvBackwardVector = vBackwardVector_;
}

void Camera::GetPrimaryRayDirection(
    __in tgir::Real const &s,  ///< [in] フィルム面上の位置(X)．値域は[0,1]
    __in tgir::Real const &t,  ///< [in] フィルム面上の位置(Y)．値域は[0,1]
    __in tgir::Vector2 const &vLensCoordinate,  ///< [in] レンズ上の座標(X,Y)
    __out tgir::Vector3 *const pvRayDirection  ///< [out] 一次光線の方向
    ) const {
  // フィルム位置から，ピントが合う面上での位置を求めて，レンズ上の位置で補正して方向を決める
  tgir::Real const fFocasPlaneSize =
      fFilmSize_ * fFocusDistance_ / fLensFilmDistance_;
  tgir::Vector3 const vLocalDirection(
      (s - 0.5) * fFocasPlaneSize * fAspectRatio_ - vLensCoordinate[0],
      (t - 0.5) * fFocasPlaneSize - vLensCoordinate[1], -fFocusDistance_);

  // 世界座標系に変換した後，正規化する
  (*pvRayDirection) = hi::normalize(vLocalDirection[0] * vRightwardVector_ +
                                    vLocalDirection[1] * vUpwardVector_ +
                                    vLocalDirection[2] * vBackwardVector_);
}

/**
 *
 * @retval true フィルム範囲内
 * @retval false フィルム範囲外
 */
bool Camera::GetFilmPosition(
    __in tgir::Vector2 const &vLensCoordinate,  ///< [in] レンズ上の座標(X,Y)
    __in tgir::Vector3 const &vDirection,       ///< [in]
                                                ///レンズ上の座標から三次元空間内の物体表面上の点へのベクトル
    __out tgir::Vector2
        *const pvFilmPosition  ///< [out] フィルム上の座標(X∈[0,1],Y∈[0,1])
    ) const {
  // 点の位置をローカル座標系に変換する
  tgir::Vector3 const vLocalDirection(hi::dot(vDirection, vRightwardVector_),
                                      hi::dot(vDirection, vUpwardVector_),
                                      hi::dot(vDirection, vBackwardVector_));

  assert(vLocalDirection[2] < 0);

  // ピントが合う面上での座標を求める
  tgir::Real const fToFocusPlane = fFocusDistance_ / -vLocalDirection[2];
  tgir::Real const fFocusPlaneX =
      vLocalDirection[0] * fToFocusPlane + vLensCoordinate[0];
  tgir::Real const fFocusPlaneY =
      vLocalDirection[1] * fToFocusPlane + vLensCoordinate[1];

  // フィルム面上に射影する(合同変換・フィルムサイズで値域を[0,1]に変換)
  tgir::Real const fToFilmPlane =
      (fLensFilmDistance_ / fFocusDistance_) / fFilmSize_;
  (*pvFilmPosition)[0] =
      tgir::Real(0.5) + fFocusPlaneX * fToFilmPlane / fAspectRatio_;
  (*pvFilmPosition)[1] = tgir::Real(0.5) + fFocusPlaneY * fToFilmPlane;

  // フィルム領域内かどうか判定して結果を返す
  return (0 <= (*pvFilmPosition)[0]) && ((*pvFilmPosition)[0] < 1) &&
         (0 <= (*pvFilmPosition)[1]) && ((*pvFilmPosition)[1] < 1);
}

}  // end of namespace tgir

/*
namespace
{
  void concentric_map(tgir::Real & u, tgir::Real & v)
  {
    u = 2 * u - 1;
    v = 2 * v - 1;

    // Avoid Zero Divide
    if ((0 == u) && (0 == v))
    {
      return;
    }

    tgir::Real radius;
    tgir::Real theta;

    if (u > -v)
    {
      if (u > v) { radius =  u; theta = (u == 0) ? 0 : 0 + v / u; }
      else         { radius =  v; theta = (v == 0) ? 0 : 2 - u / v; }
    }
    else
    {
      if (u < v) { radius = -u; theta = (u == 0) ? 0 : 4 + v / u; }
      else         { radius = -v; theta = (v == 0) ? 0 : 6 - u / v; }
    }

    theta *= (M_PI/4);

    u = radius * std::cos(theta);
    v = radius * std::sin(theta);
  }

} // end of unnamed namespace
*/
