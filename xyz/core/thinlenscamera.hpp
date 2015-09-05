#ifndef XYZ_THINLENSCAMERA_HPP_
#define XYZ_THINLENSCAMERA_HPP_

namespace xyz {
class PrimarySample;

/**
 * The thin-lens ThinLensCamera
 * 光学定理: 1/L + 1/V = 1/f
 *   L: レンズからピントが合う面までの距離
 *   V: フィルム面からレンズまでの距離
 *   f: 焦点距離(レンズの中心位置から平行光線が焦点を結ぶ位置までの距離)
 * F値: F = f/2R
 *   R: レンズの半径(=f/2F)
 */
class ThinLensCamera {
 public:
  ThinLensCamera();

  //
  // GUI
  //
 public:
  void Render() const;

  inline float3_t Eye() const { return vPrincipalPoint_; }
  inline float3_t At() const { return vGazingPoint_; }

  void MoveTo(float3_t const &v);
  void MoveBy(float3_t const &v);
  void LookAt(float3_t const &v);
  void LookBy(float3_t const &v);

  inline void SetAspectRatio(float_t const &aspect_ratio) {
    fAspectRatio_ = aspect_ratio;
  }

  inline float_t GetLensArea() const {
    return hi::square_of(fLensRadius_) * M_PI;
  }

  /**
   * ピンホールカメラモデルにおける初期視線ベクトルの
   * サンプリング用の確率密度分布(測度ds(w)) p(w) について．
   *
   * フィルム平面上のサンプリング用の確率密度分布関数(測度dA(x,y))は
   *   p(x,y) = 1/(W*H) ... (フィルム平面上の一様分布の場合)
   * である．ここで，w は世界座標系におけるフィルムの幅(m単位)で，
   * 同様に h はフィルムの高さ(同じくm単位)である．
   * つまり，フィルム面上のある微小領域 dA(x,y) がサンプリングされる確率は
   *   p(x,y) dA(x,y) = 1/(W*H) dA(x,y)
   * となる．
   *
   * 球面における方向のサンプリングの測度(ds(w))から
   * 平面における点のサンプリングの測度(dA(x))への変換式は，
   * 立体角の定義から ds(w) = cos(theta)/r^2 dA(x) となる．
   * ここで，theta は，方向 w と平面の法線との角度である．
   * また，D は，方向のサンプリングの始点から
   * 平面上のサンプリング点までの距離である．
   *
   * つまり，フィルム平面上で一様分布は，方向に関する確率密度関数に変換すると
   *   p(x,y) dA(x,y)
   *     = 1/(W*H) dA(x,y)                 ...
   *サンプリング用の確率密度関数は平面上の一様分布
   *     = 1/(W*H) D^2/cos(theta) ds(w)    ... 測度を面積測度から方向測度に変換
   *     = V^2/(W*H) 1/cos^3(theta) ds(w)  ... V = D * cos(theta)
   *を代入して，変数部と定数部をまとめた
   * ここで，V はレンズ・フィルム間の距離である．
   *
   * p(x,y) dA(x,y) = p(w) ds(w) より，
   * p(w) = V^2/(W*H) 1/cos^3(theta) であり，
   * 射影確率密度関数は p⊥(w) = d^2/(W*H) 1/cos^4(theta)
   * となる．
   *
   * この関数は，この定数部 V^2/(W*H) を返す．
   *
   * 当然，フォーカスが合う面でサンプリングした場合も同じ値になる．
   * W' = (L/V)*W  ... フォーカスがあう面の幅
   * H' = (L/V)*H  ... フォーカスがあう面の高さ
   *
   * p(w) = L^2/(W'*H') 1/cos^3(theta)
   *      = L^2/((L/V)^2*W*H) 1/cos^3(theta)
   *      = V^2/(W*H) 1/cos^3(theta)
   */
  inline float_t GetConstFactor() const {
    return hi::square_of(fFocalLength_ / fFilmSize_) / fAspectRatio_;
  }

  inline float_t FluxToRadianceCoefficient() const { return GetConstFactor(); }

 private:
  void Reset();

  //
  // renderer
  //
 public:
  //==========================================================================
  //
  // Pinhole camera model.
  //
  //==========================================================================

  void GetPrimaryRayOriginAndCameraBasis(
      __out float3_t
          *const pvRayOrigin,  ///< 光線の始点(三次元空間内のレンズ上の位置)
      __out float3_t *const pvRightwardVector,  ///< カメラの右方向ベクトル
      __out float3_t *const pvUpwardVector,     ///< カメラの上方向ベクトル
      __out float3_t *const pvBackwardVector    ///< カメラの後方向ベクトル
      ) const;

  void GetPrimaryRayDirection(
      __in float_t const &s,  ///< フィルム面上の位置(X)．値域は[0,1]
      __in float_t const &t,  ///< フィルム面上の位置(Y)．値域は[0,1]
      __out float3_t *const pvRayDirection  ///< 一次光線の方向
      ) const;

  bool GetFilmPosition(
      __in float3_t const &vDirection,  ///< [in]
      ///レンズ上の座標から三次元空間内の物体表面上の点へのベクトル
      __out float2_t
          *const pvFilmPosition  ///< [out] フィルム上の座標(X∈[0,1],Y∈[0,1])
      ) const;

  //==========================================================================
  //
  // Thin-lens camera model.
  //
  //==========================================================================

  void GetPrimaryRayOriginAndCameraBasis(
      __in float_t const &s,  ///< [in]
      ///レンズ上の座標を計算するためのパラメータ(動径用)．値域は[0,1]
      __in float_t const &t,  ///< [in]
      ///レンズ上の座標を計算するためのパラメータ(角度用)．値域は[0,1]
      __out float2_t *const vLensCoordinate,  ///< [out] レンズ上の座標(X,Y)
      __out float3_t *const
          vRayOrigin,  ///< [out] 光線の始点(三次元空間内のレンズ上の位置)
      __out float3_t *const vRightwardVector,  ///< [out] カメラの右方向ベクトル
      __out float3_t *const vUpwardVector,     ///< [out] カメラの上方向ベクトル
      __out float3_t *const vBackwardVector    ///< [out] カメラの後方向ベクトル
      ) const;

  void GetPrimaryRayDirection(
      __in float_t const &s,  ///< [in] フィルム面上の位置(X)．値域は[0,1]
      __in float_t const &t,  ///< [in] フィルム面上の位置(Y)．値域は[0,1]
      __in float2_t const &vLensCoordinate,  ///< [in] レンズ上の座標(X,Y)
      __out float3_t *const pvRayDirection   ///< [out] 一次光線の方向
      ) const;

  bool GetFilmPosition(
      __in float2_t const &vLensCoordinate,  ///< [in] レンズ上の座標(X,Y)
      __in float3_t const &vDirection,       ///< [in]
      ///レンズ上の座標から三次元空間内の物体表面上の点へのベクトル
      __out float2_t
          *const pvFilmPosition  ///< [out] フィルム上の座標(X∈[0,1],Y∈[0,1])
      ) const;

 private:
  float3_t vPrincipalPoint_;  ///< レンズの主点の座標
  float3_t vGazingPoint_;     ///< 注視点

  float_t fFNumber_;      ///< レンズ特性: F値
  float_t fFocalLength_;  ///< レンズ特性: 焦点距離(f)

  float_t fFilmSize_;     ///< フィルム特性: フィルムサイズ(高さ(H))
  float_t fAspectRatio_;  ///< フィルム特性:
  ///フィルムのアスペクト比(フィルムの幅(H)＝フィルムの高さ(H)×アスペクト比(A))

  float_t fFocusDistance_;  ///<
  ///レンズが置いてある平面とピントが合う平面の間の距離(L)
  float_t fLensFilmDistance_;  ///< レンズ・フィルム間の距離(V)
  float_t fLensRadius_;        ///< レンズの有効半径(R)

  float3_t vRightwardVector_;  ///< カメラの右方向ベクトル
  float3_t vUpwardVector_;     ///< カメラの上方向ベクトル
  float3_t vBackwardVector_;   ///< カメラの後方向ベクトル
};

}  // end of namespace xyz

#endif XYZ_CAMERA_HPP_
