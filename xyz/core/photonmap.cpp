#include "PhotonMap.hpp"
#include "../geom/pathvertex.hpp"

#include "../bsdf/bsdf.hpp"

namespace {
xyz::float_t const EPSILON_IMPORTON = 1;

}  // end of unnamed namespace

namespace xyz {
void ImportonQuery::begin_query() {
  for (std::size_t i = 1; i <= CONFIG_N_SIZE; ++i) {
    cdfDiffuse_[i] = EPSILON_IMPORTON;
  }
}

void ImportonQuery::push_back(Importance const &data,
                              Importance::vector_type const &vDistance) {
  float_t const fCosTheta = -hi::dot(data.dir(), pathVertex_.vGeometricNormal);
  if (fCosTheta <= 0) {
    return;
  }

  float3_t const vPoint =
      vDistance + data.dir() * hi::dot(vDistance, pathVertex_.vGeometricNormal);
  if (hi::length_squared(vPoint) >= ImportonQuery::max_search_range()) {
    return;
  }

  // (u, v) = (cos^2\theta, \phi/2pi)
  int u = static_cast<int>(hi::square_of(fCosTheta) * CONFIG_U_SIZE);
  if (u >= CONFIG_U_SIZE) {
    u = CONFIG_U_SIZE - 1;
  }

  float_t phi =
      float_t(0.5) * std::atan2(-hi::dot(data.dir(), pathVertex_.vBinormal),
                                -hi::dot(data.dir(), pathVertex_.vTangent));
  if (phi < 0) {
    phi += M_PI;
  }

  int v = static_cast<int>(phi * (CONFIG_V_SIZE * M_1_PI));
  if (v >= CONFIG_V_SIZE) {
    v = CONFIG_V_SIZE - 1;
  }

  ++cdfDiffuse_[1 + u * CONFIG_V_SIZE + v];
}

void ImportonQuery::end_query() {
  for (std::size_t i = 2; i <= CONFIG_N_SIZE; ++i) {
    cdfDiffuse_[i] += cdfDiffuse_[i - 1];
  }
}

void ImportonQuery::SetDiffusedVector(std::mt19937_64 &random) {
  float_t const value = cdfDiffuse_[CONFIG_N_SIZE] * random.next<float_t>();
  std::vector<float_t>::const_iterator const it =
      std::upper_bound(cdfDiffuse_.begin(), cdfDiffuse_.end(), value);

  std::size_t const n = it - cdfDiffuse_.begin() - 1;
  std::size_t const u = n / CONFIG_V_SIZE;
  std::size_t const v = n % CONFIG_V_SIZE;

  float_t const rMuU = (u + random.next<float_t>()) / CONFIG_U_SIZE;
  float_t const rMuV = (v + random.next<float_t>()) / CONFIG_V_SIZE;

  float_t const fCosTheta = std::sqrt(rMuU);
  float_t const rSinTheta = std::sqrt(1 - rMuU);
  float_t const rPhi = float_t(2 * M_PI) * rMuV;

  pathVertex_.vOutgoingDirection =
      pathVertex_.vTangent * (rSinTheta * std::cos(rPhi)) +
      pathVertex_.vBinormal * (rSinTheta * std::sin(rPhi)) +
      pathVertex_.vGeometricNormal * (fCosTheta);

  pathVertex_.power *= cdfDiffuse_[CONFIG_N_SIZE] /
                       (CONFIG_N_SIZE * (cdfDiffuse_[n + 1] - cdfDiffuse_[n]));
}

}  // end of namespace xyz

namespace xyz {
void ParticleQuery::begin_query() {
  for (std::size_t i = 0; i < CONFIG_N_SIZE; ++i) {
    cdfDiffuse_[uVertexIndexOffset_ + i] = EPSILON_IMPORTON;
  }
}

void ParticleQuery::push_back(Importance const &data,
                              Importance::vector_type const &vDistance) {
  float_t const fCosTheta =
      -hi::dot(data.dir(), pathVertices_[uVertexIndex_].vGeometricNormal);
  if (fCosTheta <= 0) {
    return;
  }

  float3_t const vPoint =
      vDistance +
      data.dir() *
          hi::dot(vDistance, pathVertices_[uVertexIndex_].vGeometricNormal);
  if (hi::length_squared(vPoint) >= ImportonQuery::max_search_range()) {
    return;
  }

  // (u, v) = (cos^2\theta, \phi/2pi)
  int u = static_cast<int>(hi::square_of(fCosTheta) * CONFIG_U_SIZE);
  if (u >= CONFIG_U_SIZE) {
    u = CONFIG_U_SIZE - 1;
  }

  float_t phi =
      float_t(0.5) *
      std::atan2(-hi::dot(data.dir(), pathVertices_[uVertexIndex_].vBinormal),
                 -hi::dot(data.dir(), pathVertices_[uVertexIndex_].vTangent));
  if (phi < 0) {
    phi += M_PI;
  }

  int v = static_cast<int>(phi * (CONFIG_V_SIZE * M_1_PI));
  if (phi >= CONFIG_V_SIZE) {
    phi = CONFIG_V_SIZE - 1;
  }

  ++cdfDiffuse_[uVertexIndexOffset_ + u * CONFIG_V_SIZE + v];
}

void ParticleQuery::end_query() {
  for (std::size_t i = 0; i < CONFIG_N_SIZE; ++i) {
    cdfDiffuse_[uVertexIndexOffset_ + i] +=
        cdfDiffuse_[uVertexIndexOffset_ + i - 1];
  }
}

void ParticleQuery::SetDiffusedVector(std::size_t nRayCount,
                                      std::mt19937_64 &random) {
  std::size_t const uLastIndex = nRayCount * CONFIG_N_SIZE;
  std::vector<float_t>::const_iterator const begin = cdfDiffuse_.begin();
  std::vector<float_t>::const_iterator const end = begin + uLastIndex + 1;

  for (std::size_t j = 0; j < nRayCount; ++j) {
    float_t const value =
        cdfDiffuse_[uLastIndex] * (j + random.next<float_t>()) / nRayCount;
    std::vector<float_t>::const_iterator const it =
        std::upper_bound(begin, end, value);

    std::size_t const n = it - begin - 1;
    std::size_t const i = (n / CONFIG_N_SIZE);
    std::size_t const u = (n % CONFIG_N_SIZE) / CONFIG_V_SIZE;
    std::size_t const v = (n % CONFIG_N_SIZE) % CONFIG_V_SIZE;

    pathVertices_[j] = tempVertices_[i];

    float_t const rMuU = (u + random.next<float_t>()) / CONFIG_U_SIZE;
    float_t const rMuV = (v + random.next<float_t>()) / CONFIG_V_SIZE;

    float_t const fCosTheta = std::sqrt(rMuU);
    float_t const rSinTheta = std::sqrt(1 - rMuU);
    float_t const rPhi = float_t(2 * M_PI) * rMuV;

    pathVertices_[j].vOutgoingDirection =
        pathVertices_[j].vTangent * (rSinTheta * std::cos(rPhi)) +
        pathVertices_[j].vBinormal * (rSinTheta * std::sin(rPhi)) +
        pathVertices_[j].vGeometricNormal * (fCosTheta);

    pathVertices_[j].power *=
        cdfDiffuse_[uLastIndex] /
        (uLastIndex * (cdfDiffuse_[n + 1] - cdfDiffuse_[n]));
  }
}

}  // end of namespace xyz

namespace xyz {
void PhotonQuery::begin_query() { fPower_ = 0; }

void PhotonQuery::push_back(Photon const &data,
                            Photon::vector_type const &vDistance) {
  float_t const fCosTheta = -hi::dot(data.dir(), pathVertex_.vGeometricNormal);
  if (fCosTheta <= 0) {
    return;
  }

  float3_t const vPoint =
      vDistance + data.dir() * hi::dot(vDistance, pathVertex_.vGeometricNormal);
  float_t const rDistanceSquared = hi::length_squared(vPoint);
  if (rDistanceSquared >= PhotonQuery::max_search_range()) {
    return;
  }
#if 0  // フィルタをかけるか
#if 1  // コーンフィルタか
    static float_t const k = 1.1; // 円錐フィルタのパラメータ
    float_t const w = 1 - std::sqrt(rDistanceSquared) / (k * PhotonQuery::max_distance());
#else
    static float_t const alpha = 0.918; // ガウスフィルタのパラメータ
    static float_t const beta = 1.953;
    float_t const u = 1 - std::exp(-beta * rDistanceSquared / (2 * PhotonQuery::max_search_range()));
    float_t const s = 1 - std::exp(-beta);
    float_t const w = alpha * (1 - u / s);
#endif
    fPower_ += data.pow() * w;
#else
  fPower_ += data.pow();
#endif
}

void PhotonQuery::push_back(Importance const &data,
                            Importance::vector_type const &vDistance) {
  float_t const fCosTheta = -hi::dot(data.dir(), pathVertex_.vGeometricNormal);
  if (fCosTheta <= 0) {
    return;
  }

  float3_t const vPoint =
      vDistance + data.dir() * hi::dot(vDistance, pathVertex_.vGeometricNormal);
  float_t const rDistanceSquared = hi::length_squared(vPoint);
  if (rDistanceSquared >= PhotonQuery::max_search_range()) {
    return;
  }
#if 0  // フィルタをかけるか
#if 1  // コーンフィルタか
    static float_t const k = 1.1; // 円錐フィルタのパラメータ
    float_t const w = 1 - std::sqrt(rDistanceSquared) / (k * PhotonQuery::max_distance());
#else
    static float_t const alpha = 0.918; // ガウスフィルタのパラメータ
    static float_t const beta = 1.953;
    float_t const u = 1 - std::exp(-beta * rDistanceSquared / (2 * PhotonQuery::max_search_range()));
    float_t const s = 1 - std::exp(-beta);
    float_t const w = alpha * (1 - u / s);
#endif
    fPower_ += (data.dir() * 0.5 + 0.5) * w;
#else
  fPower_ += (data.dir() * 0.5 + 0.5);
#endif
}

void PhotonQuery::end_query() {}

float3_t PhotonQuery::GetEstimatedRadiance() {
  float3_t rgb;
#if 0  // フィルタをかけるか
#if 1  // コーンフィルタか
    static float_t const k = 1.1; // 円錐フィルタのパラメータ
    static float_t const n = 1.0 - 2.0 / (3.0 * k); // 円錐フィルタの正規化定数
    rgb = fPower_ / (n * M_PI * PhotonQuery::max_search_range());
#else
    rgb = fPower_;
#endif
#else
  rgb = fPower_ * M_1_PI / PhotonQuery::max_search_range();
#endif
  float3_t xyz;
  hi::CCIR601_1_RGB2XYZ(rgb, xyz);
  return xyz;
}

}  // end of namespace xyz

namespace xyz {
PhotonMapProperty::PhotonMapProperty()
    : nPhotons_(0),
      rPowerAvg_(0),
      rPowerVar_(0),
      rPowerMin_(std::numeric_limits<float_t>::infinity()),
      rPowerMax_(-std::numeric_limits<float_t>::infinity()) {}

void PhotonMapProperty::operator()(Photon const &photon) {
  ++nPhotons_;

  float_t const rPower = photon.pow()[1];

  rPowerAvg_ += rPower;
  rPowerVar_ += hi::square_of(rPower);

  if (rPowerMin_ > rPower) {
    rPowerMin_ = rPower;
  }

  if (rPowerMax_ < rPower) {
    rPowerMax_ = rPower;
  }
}

void PhotonMapProperty::operator()(Importance const &impoton) {
  UNREFERENCED_PARAMETER(impoton);

  ++nPhotons_;

  float_t const fPower = 1;

  rPowerAvg_ += fPower;
  rPowerVar_ += hi::square_of(fPower);

  if (rPowerMin_ > fPower) {
    rPowerMin_ = fPower;
  }

  if (rPowerMax_ < fPower) {
    rPowerMax_ = fPower;
  }
}

void PhotonMapProperty::print() {
  rPowerAvg_ /= nPhotons_;
  rPowerVar_ /= nPhotons_;
  rPowerVar_ -= hi::square_of(rPowerAvg_);
  rPowerVar_ = std::abs(rPowerVar_);

  ::_ftprintf_s(stdout, _TEXT("  >> NUM: %u\n"), nPhotons_);
  ::_ftprintf_s(stdout, _TEXT("  >> VALUE: [%g, %g]\n"), rPowerAvg_,
                rPowerVar_);
  ::_ftprintf_s(stdout, _TEXT("  >> RANGE: [%g, %g]\n"), rPowerMin_,
                rPowerMax_);
}

}  // end of namespace xyz

namespace xyz {
void ImportanceQuery::begin_query() {
  for (std::size_t i = 1; i <= CONFIG_N_SIZE; ++i) {
#if 0
      cdfDiffuse_[i] = EPSILON_IMPORTON;
#else
    cdfDiffuse_[i] = CONFIG_N_SIZE;
#endif
  }
}

void ImportanceQuery::push_back(Importance const &data,
                                Importance::vector_type const &vDistance) {
  float_t const fCosTheta = -hi::dot(data.dir(), pathVertex_.vShadingNormal);
  if (fCosTheta <= 0) {
    return;
  }

  float3_t const vPoint =
      vDistance + data.dir() * hi::dot(vDistance, pathVertex_.vShadingNormal);
  if (hi::length_squared(vPoint) >= ImportonQuery::max_search_range()) {
    return;
  }

  // (u, v) = (cos^2\theta, \phi/2pi)
  int u = static_cast<int>(hi::square_of(fCosTheta) * CONFIG_U_SIZE);
  if (u >= CONFIG_U_SIZE) {
    u = CONFIG_U_SIZE - 1;
  }

  float_t phi =
      float_t(0.5) * std::atan2(-hi::dot(data.dir(), pathVertex_.vBinormal),
                                -hi::dot(data.dir(), pathVertex_.vTangent));
  if (phi < 0) {
    phi += M_PI;
  }

  int v = static_cast<int>(phi * (CONFIG_V_SIZE * M_1_PI));
  if (v >= CONFIG_V_SIZE) {
    v = CONFIG_V_SIZE - 1;
  }

  ++cdfDiffuse_[1 + u * CONFIG_V_SIZE + v];
}

void ImportanceQuery::end_query() {
#if 1
  /*
      cdfU_[0] = 0;
      for (std::size_t u = 0; u < CONFIG_U_SIZE; ++u)
      {
        // P(v|u) を計算する
        cdfV_[u*(CONFIG_V_SIZE+1)] = 0;
        for (std::size_t v = 1; v <= CONFIG_V_SIZE; ++v)
        {
          cdfV_[u*(CONFIG_V_SIZE+1)+v] = cdfDiffuse_[u*CONFIG_V_SIZE+v] +
     cdfV_[u*(CONFIG_V_SIZE+1)+v-1];
        }

        // P(u) を計算する
        cdfU_[u+1] = cdfV_[u*(CONFIG_V_SIZE+1)+CONFIG_V_SIZE] + cdfU_[u];
      }
  */
  for (std::size_t i = 2; i <= CONFIG_N_SIZE; ++i) {
    cdfDiffuse_[i] += cdfDiffuse_[i - 1];
  }
/*
    // ある一定数以上集まってないと推定に使わない
    if (cdfDiffuse_[CONFIG_N_SIZE] < 2 * CONFIG_N_SIZE)
    {
      for (std::size_t i = 1; i <= CONFIG_N_SIZE; ++i)
      {
        cdfDiffuse_[i] = i;
      }
    }
*/
#elif 1
  // 逆に考えるんだ，明るいところは逆にサンプルする必要がないと考えるんだ
  float_t fSumOfPower = float_t(0);
  for (std::size_t i = 1; i <= CONFIG_N_SIZE; ++i) {
    fSumOfPower += cdfDiffuse_[i];
  }
  for (std::size_t i = 1; i <= CONFIG_N_SIZE; ++i) {
    cdfDiffuse_[i] = fSumOfPower - cdfDiffuse_[i];
    // cdfDiffuse_[i]  = fSumOfPower / cdfDiffuse_[i];
    cdfDiffuse_[i] += cdfDiffuse_[i - 1];
  }
#else
  float_t fSumOfPhoton = float_t(0);
  for (std::size_t i = 1; i <= CONFIG_N_SIZE; ++i) {
    fSumOfPhoton += cdfDiffuse_[i];
  }
  if (fSumOfPhoton > 0) {
    fSumOfPhoton /= CONFIG_N_SIZE;

    float_t const fAlpha = std::exp(-fSumOfPhoton);
    assert(fAlpha > 0);
    assert(fAlpha <= 1);
    for (std::size_t i = 1; i <= CONFIG_N_SIZE; ++i) {
      cdfDiffuse_[i] += (fSumOfPhoton - cdfDiffuse_[i]) * fAlpha;
      cdfDiffuse_[i] += cdfDiffuse_[i - 1];
    }
  } else {
    for (std::size_t i = 1; i <= CONFIG_N_SIZE; ++i) {
      cdfDiffuse_[i] = i;
    }
  }
#endif
}

float_t ImportanceQuery::GetDensity(float3_t const &dir) const {
  float_t const fCosTheta = -hi::dot(dir, pathVertex_.vShadingNormal);
  if (fCosTheta <= 0) {
    return 0.0;
  }

  // (u, v) = (cos^2\theta, \phi/2pi)
  int u = static_cast<int>(hi::square_of(fCosTheta) * CONFIG_U_SIZE);
  if (u >= CONFIG_U_SIZE) {
    u = CONFIG_U_SIZE - 1;
  }

  float_t phi = float_t(0.5) * std::atan2(-hi::dot(dir, pathVertex_.vBinormal),
                                          -hi::dot(dir, pathVertex_.vTangent));
  if (phi < 0) {
    phi += M_PI;
  }

  int v = static_cast<int>(phi * (CONFIG_V_SIZE * M_1_PI));
  if (v >= CONFIG_V_SIZE) {
    v = CONFIG_V_SIZE - 1;
  }

  std::size_t const n = u * CONFIG_V_SIZE + v;

  return CONFIG_N_SIZE * (cdfDiffuse_[n + 1] - cdfDiffuse_[n]) /
         cdfDiffuse_[CONFIG_N_SIZE];
}

void ImportanceQuery::SetDiffusedVector(__inout class IPrimarySample &sample,
                                        __out float3_t *const pVector,
                                        __out float_t *const pCosTheta,
                                        __out float_t *const pPDF,
                                        __out float_t *const pAlpha) const {
#if 1
  float_t const value = cdfDiffuse_[CONFIG_N_SIZE] * sample.next();
  std::vector<float_t>::const_iterator const it =
      std::upper_bound(cdfDiffuse_.begin(), cdfDiffuse_.end(), value);

  assert(it != cdfDiffuse_.begin());
  assert(it != cdfDiffuse_.end());

  std::size_t const n = it - cdfDiffuse_.begin() - 1;
  std::size_t const u = n / CONFIG_V_SIZE;
  std::size_t const v = n % CONFIG_V_SIZE;

  float_t const fMuU = (u + sample.next()) / CONFIG_U_SIZE;
  float_t const fMuV = (v + sample.next()) / CONFIG_V_SIZE;

  float_t const fCosTheta = std::sqrt(fMuU);
  float_t const fSinTheta = std::sqrt(1 - fMuU);
  float_t const fPhi = float_t(2 * M_PI) * fMuV;

  // 反射方向を設定．
  (*pVector) = pathVertex_.vTangent * (fSinTheta * std::cos(fPhi)) +
               pathVertex_.vBinormal * (fSinTheta * std::sin(fPhi)) +
               pathVertex_.vShadingNormal * (fCosTheta);

  float_t const fDensity = CONFIG_N_SIZE *
                           (cdfDiffuse_[n + 1] - cdfDiffuse_[n]) /
                           cdfDiffuse_[CONFIG_N_SIZE];

  (*pCosTheta) = fCosTheta;
  (*pPDF) = fDensity * fCosTheta * M_1_PI;
  (*pAlpha) = hi::rcp(fDensity);
#else
  // p(u) からサンプリング
  float_t const du = cdfU_[CONFIG_U_SIZE] * sample.next();
  std::vector<float_t>::const_iterator const itU =
      std::upper_bound(cdfU_.begin(), cdfU_.end(), du);
  assert(itU != cdfU_.begin());
  assert(itU != cdfU_.end());
  std::size_t const u = itU - cdfU_.begin() - 1;

  // p(v|u) からのサンプリング
  std::size_t const bv = u * (CONFIG_V_SIZE + 1);
  std::size_t const ev = bv + CONFIG_V_SIZE;
  float_t const dv = cdfV_[ev] * sample.next();
  std::vector<float_t>::const_iterator const itV =
      std::upper_bound(cdfV_.begin() + bv, cdfV_.begin() + ev + 1, dv);
  assert(itV != cdfV_.begin() + bv);
  assert(itV != cdfV_.begin() + ev + 1);
  std::size_t const v = itV - (cdfV_.begin() + bv) - 1;

  float_t const fMuU =
      (u + (du - cdfU_[u]) / (cdfU_[u + 1] - cdfU_[u])) / CONFIG_U_SIZE;
  float_t const fMuV =
      (v + (dv - cdfV_[v]) / (cdfV_[v + 1] - cdfV_[v])) / CONFIG_V_SIZE;

  float_t const fCosTheta = std::sqrt(fMuU);
  float_t const fSinTheta = std::sqrt(1 - fMuU);
  float_t const fPhi = float_t(2 * M_PI) * fMuV;

  pathVertex_.vOutgoingDirection =
      pathVertex_.vTangent * (fSinTheta * std::cos(fPhi)) +
      pathVertex_.vBinormal * (fSinTheta * std::sin(fPhi)) +
      pathVertex_.vShadingNormal * (fCosTheta);

  std::size_t const n = u * CONFIG_V_SIZE + v;
  float_t const fDensity = CONFIG_N_SIZE *
                           (cdfDiffuse_[n + 1] - cdfDiffuse_[n]) /
                           cdfDiffuse_[CONFIG_N_SIZE];

  pathVertex_.power /= fDensity;

  return fDensity * fCosTheta * M_1_PI;
#endif
}

}  // end of namespace xyz
