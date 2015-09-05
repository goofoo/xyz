#include "glass.hpp"
#include "../core/scene.hpp"
#include "../geom/pathvertex.hpp"
#include "../geom/intersection.hpp"

namespace xyz {
void Glass::Render() const {
  hi::basic_vector4<GLfloat> const color(0.3465f, 0.222f, 0.4402f,
                                         0.5f);  // Patch 10, Purple
  ::glEnable(GL_LIGHTING);
  ::glEnable(GL_BLEND);
  ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  hi::glMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color.data());
}

bool Glass::BoundsImportance(PathVertex *const pLightPathVertex,
                             std::mt19937_64 &random) const {
  // vOutgoingDirection には 入射してきた方向が入っている
  float_t const fCosTheta_s = -hi::dot(pLightPathVertex->vOutgoingDirection,
                                       pLightPathVertex->vShadingNormal);
  if (fCosTheta_s <= 0) {
    return false;
  }
  float_t const fCosTheta_g = -hi::dot(pLightPathVertex->vOutgoingDirection,
                                       pLightPathVertex->vGeometricNormal);

  float_t const ior = 1.5;
  float_t const ior_sq =
      hi::square_of(pLightPathVertex->bBackSide ? hi::rcp(ior) : ior);

  float_t const nk = fCosTheta_s;
  float_t const b = ior_sq + nk * nk - 1;

  pLightPathVertex->bSpecular = true;

  if (b > 0) {
    float_t const g = std::sqrt(b);
    float_t const gpc = g + nk;
    float_t const gmc = g - nk;
    float_t const alpha = gmc / gpc;
    float_t const beta = (nk * gpc - 1) / (nk * gmc + 1);
    float_t const fr = float_t(0.5) * alpha * alpha * (beta * beta + 1);

    if (random.next<float_t>() < fr) {
      ComputeReflectedVector(pLightPathVertex->vShadingNormal, nk,
                             &pLightPathVertex->vOutgoingDirection);  // 反射

      float_t const fCosTheta_o = hi::dot(pLightPathVertex->vOutgoingDirection,
                                          pLightPathVertex->vGeometricNormal);
      if (fCosTheta_o <= 0) {
        return false;
      }

      // NOTE: ここでエネルギーが増えている可能性がある
      float_t const fShadingFactor = fCosTheta_o / fCosTheta_g;
      if (random.next<float_t>() >= fShadingFactor) {
        return false;
      }
    } else {
      ComputeRefractedVector(pLightPathVertex->vShadingNormal, gmc,
                             std::sqrt(ior_sq),
                             &pLightPathVertex->vOutgoingDirection);  // 屈折

      float_t const fCosTheta_o = -hi::dot(pLightPathVertex->vOutgoingDirection,
                                           pLightPathVertex->vGeometricNormal);
      if (fCosTheta_o <= 0) {
        return false;
      }

      // NOTE: ここでエネルギーが増えている可能性がある
      float_t const fShadingFactor = fCosTheta_o / fCosTheta_g;
      if (random.next<float_t>() >= fShadingFactor) {
        return false;
      }
    }
  } else {
    ComputeReflectedVector(pLightPathVertex->vShadingNormal, nk,
                           &pLightPathVertex->vOutgoingDirection);  // 完全反射

    float_t const fCosTheta_o = hi::dot(pLightPathVertex->vOutgoingDirection,
                                        pLightPathVertex->vGeometricNormal);
    if (fCosTheta_o <= 0) {
      return false;
    }

    // NOTE: ここでエネルギーが増えている可能性がある
    float_t const fShadingFactor = fCosTheta_o / fCosTheta_g;
    if (random.next<float_t>() >= fShadingFactor) {
      return false;
    }
  }

  return true;
}

float3_t Glass::CalculateWeightedDirectLighting(
    __inout PathVertex *const pEyePathVertex,
    __in PathVertex const &lightPathVertex, Scene const &) const {
  UNREFERENCED_PARAMETER(pEyePathVertex);
  UNREFERENCED_PARAMETER(lightPathVertex);
  return float3_t(0);
}

bool Glass::NextScatteringDirection(__inout PathVertex *const pEyePathVertex,
                                    __inout IPrimarySample &sample) const {
  /*
      static float_t const B1 = 1.738484030;
      static float_t const B2 = 3.11168974e-1;
      static float_t const B3 = 1.174908710;
      static float_t const C1 = 1.36068604e-2;
      static float_t const C2 = 6.15960463e-2;
      static float_t const C3 = 1.21922711e+2;
      float_t const L_sq = hi::square_of(0.4 + (0.7 - 0.4) * wavelength); // μm
      float_t const mu_sq = (B1*L_sq/(L_sq-C1)) + (B2*L_sq/(L_sq-C2)) +
     (B3*L_sq/(L_sq-C3)) + 1;
      float_t const ior_sq = pathVertex.back ? hi::rcp(mu_sq) : mu_sq;
  */
  // float_t const ior = 1.400 + 50.0 / ((400-230) + (700-400) * wavelength);
  float_t const ior = 1.5;
  float_t const ior_sq =
      hi::square_of(pEyePathVertex->bBackSide ? hi::rcp(ior) : ior);

  float_t const nk = -hi::dot(pEyePathVertex->vOutgoingDirection,
                              pEyePathVertex->vShadingNormal);
  float_t const b = ior_sq + nk * nk - 1;

  if (b > 0) {
    float_t const g = std::sqrt(b);
    float_t const gpc = g + nk;
    float_t const gmc = g - nk;
    float_t const alpha = gmc / gpc;
    float_t const beta = (nk * gpc - 1) / (nk * gmc + 1);
    float_t const fr = float_t(0.5) * alpha * alpha * (beta * beta + 1);

    if (sample.next() < fr) {
      ComputeReflectedVector(pEyePathVertex->vShadingNormal, nk,
                             &pEyePathVertex->vOutgoingDirection);  // 反射
      if (hi::dot(pEyePathVertex->vOutgoingDirection,
                  pEyePathVertex->vGeometricNormal) <= 0) {
        return false;  // 打ち止め
      }
    } else {
      ComputeRefractedVector(pEyePathVertex->vShadingNormal, gmc,
                             std::sqrt(ior_sq),
                             &pEyePathVertex->vOutgoingDirection);  // 屈折
      if (hi::dot(pEyePathVertex->vOutgoingDirection,
                  pEyePathVertex->vGeometricNormal) >= 0) {
        return false;  // 打ち止め
      }
    }
  } else {
    ComputeReflectedVector(pEyePathVertex->vShadingNormal, nk,
                           &pEyePathVertex->vOutgoingDirection);  // 完全反射
    if (hi::dot(pEyePathVertex->vOutgoingDirection,
                pEyePathVertex->vGeometricNormal) <= 0) {
      return false;  // 打ち止め
    }
  }

  pEyePathVertex->bSpecular = true;
  pEyePathVertex->fSamplingPrev = 1;
  pEyePathVertex->fSamplingNext = 1;

  return true;
}

float_t Glass::GetDensityVariance() const { return 2; }

bool Glass::NextScatteringDirection(PathVertex &pathVertex, float3_t const &mu,
                                    bool const &) const {
  float_t const ior = 1.4;
  float_t const ior_sq =
      hi::square_of(pathVertex.bBackSide ? hi::rcp(ior) : ior);
  float_t const nk =
      -hi::dot(pathVertex.vOutgoingDirection, pathVertex.vShadingNormal);
  float_t const b = ior_sq + nk * nk - float_t(1);

  if (b > 0) {
    float_t const g = std::sqrt(b);
    float_t const gpc = g + nk;
    float_t const gmc = g - nk;
    float_t const alpha = gmc / gpc;
    float_t const beta = (nk * gpc - 1) / (nk * gmc + 1);
    float_t const fr = float_t(0.5) * alpha * alpha * (beta * beta + 1);

    if (mu[2] < fr) {
      ComputeReflectedVector(pathVertex.vShadingNormal, nk,
                             &pathVertex.vOutgoingDirection);  // 反射
      if (hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal) <=
          0) {
        return false;  // 打ち止め
      }

      pathVertex.fSamplingPrev = fr;
      pathVertex.fSamplingNext = fr;
    } else {
      ComputeRefractedVector(pathVertex.vShadingNormal, gmc, std::sqrt(ior_sq),
                             &pathVertex.vOutgoingDirection);  // 屈折
      if (hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal) >=
          0) {
        return false;  // 打ち止め
      }

      pathVertex.fSamplingPrev = 1 - fr;
      pathVertex.fSamplingNext = 1 - fr;
    }
  } else {
    ComputeReflectedVector(pathVertex.vShadingNormal, nk,
                           &pathVertex.vOutgoingDirection);  // 完全反射
    if (hi::dot(pathVertex.vOutgoingDirection, pathVertex.vGeometricNormal) <=
        0) {
      return false;  // 打ち止め
    }

    pathVertex.fSamplingPrev = 1;
    pathVertex.fSamplingNext = 1;
  }

  pathVertex.bSpecular = true;

  return true;
}

}  // end of namespace xyz

namespace xyz {
// pPathVertex には，交点とその基底ベクトル，vIncomingDirection，
// および前の頂点と同じ power の値が設定されているものとする．
//   pPathVertex->vIncomingDirection = -pPreviousPathVertex->vOutgoingDirection;
//   pPathVertex->power              =  pPreviousPathVertex->power;
// また，条件として[hi::dot(pPathVertex->vIncomingDirection,
// pPathVertex->vShadingNormal) > 0]を課す．
//
// pPathVertex->fIncomingCosThetaShading は設定済み
// pPathVertex->fOutgoingCosThetaGeometric を設定する．
//
// NOTE: Incoming や Outgoing は光エネルギーの流れる方向を示しており
//       実装上の光線追跡の方向を示すものではない
bool Glass::LightsToCameraScatteringDirection(
    __inout PathVertex *const pPathVertex,
    __inout IPrimarySample &sample) const {
  float_t const ior = 1.5;
  float_t const ior_sq =
      hi::square_of(pPathVertex->bBackSide ? hi::rcp(ior) : ior);

  float_t const nk =
      hi::dot(pPathVertex->vIncomingDirection, pPathVertex->vShadingNormal);
  float_t const b = ior_sq + nk * nk - 1;

  pPathVertex->vOutgoingDirection = -pPathVertex->vIncomingDirection;
  pPathVertex->bSpecular = true;

  if (b > 0) {
    float_t const g = std::sqrt(b);
    float_t const gpc = g + nk;
    float_t const gmc = g - nk;
    float_t const alpha = gmc / gpc;
    float_t const beta = (nk * gpc - 1) / (nk * gmc + 1);
    float_t const fr = float_t(0.5) * alpha * alpha * (beta * beta + 1);

    if (sample.next() >= fr) {
      pPathVertex->fSamplingPrev = pPathVertex->fSamplingNext = 1 - fr;

      // 透過屈折
      {
        ComputeRefractedVector(pPathVertex->vShadingNormal, gmc,
                               std::sqrt(ior_sq),
                               &(pPathVertex->vOutgoingDirection));

        // 幾何学的法線と散乱方向ベクトルの内積を求める．
        pPathVertex->fOutgoingCosThetaGeometric =
            -hi::dot(pPathVertex->vOutgoingDirection,
                     pPathVertex->vGeometricNormal);
        if (pPathVertex->fOutgoingCosThetaGeometric <= 0) {
          return false;
        }

        // pPathVertex->fSamplingPrev *=
        // hi::rcp(pPathVertex->fIncomingCosThetaShading );
        // pPathVertex->fSamplingNext *=
        // hi::rcp(pPathVertex->fOutgoingCosThetaGeometric);

        // 幾何学的法線と入射方向ベクトルの内積を求める．
        float_t const fIncomingCosThetaGeometric = hi::dot(
            pPathVertex->vIncomingDirection, pPathVertex->vGeometricNormal);

        pPathVertex->fBSDFxIPDF = pPathVertex->fOutgoingCosThetaGeometric /
                                  fIncomingCosThetaGeometric;
        // pPathVertex->power *= albedo_;
      }
      return true;
    }

    pPathVertex->fSamplingPrev = pPathVertex->fSamplingNext = fr;
  } else {
    pPathVertex->fSamplingPrev = pPathVertex->fSamplingNext = 1;
  }

  // 鏡面反射
  {
    ComputeReflectedVector(pPathVertex->vShadingNormal, nk,
                           &(pPathVertex->vOutgoingDirection));

    // 幾何学的法線と散乱方向ベクトルの内積を求める．
    pPathVertex->fOutgoingCosThetaGeometric =
        hi::dot(pPathVertex->vOutgoingDirection, pPathVertex->vGeometricNormal);
    if (pPathVertex->fOutgoingCosThetaGeometric <= 0) {
      return false;
    }

    // pPathVertex->fSamplingPrev *=
    // hi::rcp(pPathVertex->fIncomingCosThetaShading );
    // pPathVertex->fSamplingNext *=
    // hi::rcp(pPathVertex->fOutgoingCosThetaGeometric);

    // 幾何学的法線と入射方向ベクトルの内積を求める．
    float_t const fIncomingCosThetaGeometric =
        hi::dot(pPathVertex->vIncomingDirection, pPathVertex->vGeometricNormal);

    pPathVertex->fBSDFxIPDF =
        pPathVertex->fOutgoingCosThetaGeometric / fIncomingCosThetaGeometric;
    // pPathVertex->power *= albedo_;
  }
  return true;
}

// pPathVertex には，交点とその基底ベクトル，vOutgoingDirection，
// および前の頂点と同じ power の値が設定されているものとする．
//   pPathVertex->vOutgoingDirection = -pPreviousPathVertex->vIncomingDirection;
//   pPathVertex->power              =  pPreviousPathVertex->power;
// また，条件として[hi::dot(pPathVertex->vOutgingDirection,
// pPathVertex->vShadingNormal) > 0]を課す．
//
// pPathVertex->fOutgoingCosThetaGeometric は設定済み
// pPathVertex->fIncomingCosThetaShading を設定する．
//
// NOTE: Incoming や Outgoing は光エネルギーの流れる方向を示しており
//       実装上の光線追跡の方向を示すものではない
bool Glass::CameraToLightsScatteringDirection(
    __inout PathVertex *const pPathVertex,
    __inout IPrimarySample &sample) const {
  float_t const ior = 1.5;
  float_t const ior_sq =
      hi::square_of(pPathVertex->bBackSide ? hi::rcp(ior) : ior);

  float_t const nk =
      hi::dot(pPathVertex->vOutgoingDirection, pPathVertex->vShadingNormal);
  float_t const b = ior_sq + nk * nk - 1;

  pPathVertex->vIncomingDirection = -pPathVertex->vOutgoingDirection;
  pPathVertex->bSpecular = true;

  if (b > 0) {
    float_t const g = std::sqrt(b);
    float_t const gpc = g + nk;
    float_t const gmc = g - nk;
    float_t const alpha = gmc / gpc;
    float_t const beta = (nk * gpc - 1) / (nk * gmc + 1);
    float_t const fr = float_t(0.5) * alpha * alpha * (beta * beta + 1);

    if (sample.next() >= fr) {
      pPathVertex->fSamplingPrev = pPathVertex->fSamplingNext = 1 - fr;

      // 透過屈折
      {
        ComputeRefractedVector(pPathVertex->vShadingNormal, gmc,
                               std::sqrt(ior_sq),
                               &(pPathVertex->vIncomingDirection));
        pPathVertex->fIncomingCosThetaShading =
            -hi::dot(pPathVertex->vIncomingDirection,
                     pPathVertex->vGeometricNormal);
        if (pPathVertex->fIncomingCosThetaShading <= 0) {
          return false;
        }
        // pPathVertex->fSamplingPrev *=
        // hi::rcp(pPathVertex->fIncomingCosThetaShading );
        // pPathVertex->fSamplingNext *=
        // hi::rcp(pPathVertex->fOutgoingCosThetaGeometric);
        pPathVertex->fBSDFxIPDF = 1;
      }
      return true;
    }

    pPathVertex->fSamplingPrev = pPathVertex->fSamplingNext = fr;
  } else {
    pPathVertex->fSamplingPrev = pPathVertex->fSamplingNext = 1;
  }

  // 鏡面反射
  {
    ComputeReflectedVector(pPathVertex->vShadingNormal, nk,
                           &(pPathVertex->vIncomingDirection));
    pPathVertex->fIncomingCosThetaShading =
        hi::dot(pPathVertex->vIncomingDirection, pPathVertex->vGeometricNormal);
    if (pPathVertex->fIncomingCosThetaShading <= 0) {
      return false;
    }
    // pPathVertex->fSamplingPrev *=
    // hi::rcp(pPathVertex->fIncomingCosThetaShading );
    // pPathVertex->fSamplingNext *=
    // hi::rcp(pPathVertex->fOutgoingCosThetaGeometric);
    pPathVertex->fBSDFxIPDF = 1;
    return true;
  }
}
}  // end of namespace xyz
