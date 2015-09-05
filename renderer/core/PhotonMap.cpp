#include "PhotonMap.hpp"
#include "geom/PathVertex.hpp"

namespace {
tgir::Real const EPSILON_IMPORTON = 1;

}  // end of unnamed namespace

namespace tgir {
void ImportonQuery::begin_query() {
  for (std::size_t i = 1; i <= CONFIG_N_SIZE; ++i) {
    cdfDiffuse_[i] = EPSILON_IMPORTON;
  }
}

void ImportonQuery::push_back(Importance const &data,
                              Importance::vector_type const &vDistance) {
  tgir::Real const rCosTheta =
      -hi::dot(data.dir(), pathVertex_.vGeometricNormal);
  if (rCosTheta <= 0) {
    return;
  }

  tgir::Vector3 const vPoint =
      vDistance + data.dir() * hi::dot(vDistance, pathVertex_.vGeometricNormal);
  if (hi::length_squared(vPoint) >= ImportonQuery::max_search_range()) {
    return;
  }

  // (u, v) = (cos^2\theta, \phi/2pi)
  int u = static_cast<int>(hi::square_of(rCosTheta) * CONFIG_U_SIZE);
  if (u >= CONFIG_U_SIZE) {
    u = CONFIG_U_SIZE - 1;
  }

  tgir::Real phi =
      tgir::Real(0.5) * std::atan2(-hi::dot(data.dir(), pathVertex_.vBinormal),
                                   -hi::dot(data.dir(), pathVertex_.vTangent));
  if (phi < 0) {
    phi += M_PI;
  }

  int v = static_cast<int>(phi * (CONFIG_V_SIZE * M_1_PI));
  if (phi >= CONFIG_V_SIZE) {
    phi = CONFIG_V_SIZE - 1;
  }

  ++cdfDiffuse_[1 + u * CONFIG_V_SIZE + v];
}

void ImportonQuery::end_query() {
  for (std::size_t i = 2; i <= CONFIG_N_SIZE; ++i) {
    cdfDiffuse_[i] += cdfDiffuse_[i - 1];
  }
}

void ImportonQuery::SetDiffusedVector(std::mt19937_64 &random) {
  tgir::Real const value = cdfDiffuse_[CONFIG_N_SIZE] *
                           std::uniform_real_distribution<tgir::Real>()(random);
  std::vector<tgir::Real>::const_iterator const it =
      std::upper_bound(cdfDiffuse_.begin(), cdfDiffuse_.end(), value);

  std::size_t const n = it - cdfDiffuse_.begin() - 1;
  std::size_t const u = n / CONFIG_V_SIZE;
  std::size_t const v = n % CONFIG_V_SIZE;

  tgir::Real const rMuU =
      (u + std::uniform_real_distribution<tgir::Real>()(random)) /
      CONFIG_U_SIZE;
  tgir::Real const rMuV =
      (v + std::uniform_real_distribution<tgir::Real>()(random)) /
      CONFIG_V_SIZE;

  tgir::Real const rCosTheta = std::sqrt(rMuU);
  tgir::Real const rSinTheta = std::sqrt(1 - rMuU);
  tgir::Real const rPhi = tgir::Real(2 * M_PI) * rMuV;

  pathVertex_.vOutgoingDirection =
      pathVertex_.vTangent * (rSinTheta * std::cos(rPhi)) +
      pathVertex_.vBinormal * (rSinTheta * std::sin(rPhi)) +
      pathVertex_.vGeometricNormal * (rCosTheta);

  pathVertex_.sQuantum *=
      cdfDiffuse_[CONFIG_N_SIZE] /
      (CONFIG_N_SIZE * (cdfDiffuse_[n + 1] - cdfDiffuse_[n]));
}

}  // end of namespace tgir

namespace tgir {
void ParticleQuery::begin_query() {
  for (std::size_t i = 0; i < CONFIG_N_SIZE; ++i) {
    cdfDiffuse_[uVertexIndexOffset_ + i] = EPSILON_IMPORTON;
  }
}

void ParticleQuery::push_back(Importance const &data,
                              Importance::vector_type const &vDistance) {
  tgir::Real const rCosTheta =
      -hi::dot(data.dir(), pathVertices_[uVertexIndex_].vGeometricNormal);
  if (rCosTheta <= 0) {
    return;
  }

  tgir::Vector3 const vPoint =
      vDistance +
      data.dir() *
          hi::dot(vDistance, pathVertices_[uVertexIndex_].vGeometricNormal);
  if (hi::length_squared(vPoint) >= ImportonQuery::max_search_range()) {
    return;
  }

  // (u, v) = (cos^2\theta, \phi/2pi)
  int u = static_cast<int>(hi::square_of(rCosTheta) * CONFIG_U_SIZE);
  if (u >= CONFIG_U_SIZE) {
    u = CONFIG_U_SIZE - 1;
  }

  tgir::Real phi =
      tgir::Real(0.5) *
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
  std::vector<tgir::Real>::const_iterator const begin = cdfDiffuse_.begin();
  std::vector<tgir::Real>::const_iterator const end = begin + uLastIndex + 1;

  for (std::size_t j = 0; j < nRayCount; ++j) {
    tgir::Real const value =
        cdfDiffuse_[uLastIndex] *
        (j + std::uniform_real_distribution<tgir::Real>()(random)) / nRayCount;
    std::vector<tgir::Real>::const_iterator const it =
        std::upper_bound(begin, end, value);

    std::size_t const n = it - begin - 1;
    std::size_t const i = (n / CONFIG_N_SIZE);
    std::size_t const u = (n % CONFIG_N_SIZE) / CONFIG_V_SIZE;
    std::size_t const v = (n % CONFIG_N_SIZE) % CONFIG_V_SIZE;

    pathVertices_[j] = tempVertices_[i];

    tgir::Real const rMuU =
        (u + std::uniform_real_distribution<tgir::Real>()(random)) /
        CONFIG_U_SIZE;
    tgir::Real const rMuV =
        (v + std::uniform_real_distribution<tgir::Real>()(random)) /
        CONFIG_V_SIZE;

    tgir::Real const rCosTheta = std::sqrt(rMuU);
    tgir::Real const rSinTheta = std::sqrt(1 - rMuU);
    tgir::Real const rPhi = tgir::Real(2 * M_PI) * rMuV;

    pathVertices_[j].vOutgoingDirection =
        pathVertices_[j].vTangent * (rSinTheta * std::cos(rPhi)) +
        pathVertices_[j].vBinormal * (rSinTheta * std::sin(rPhi)) +
        pathVertices_[j].vGeometricNormal * (rCosTheta);

    pathVertices_[j].sQuantum *=
        cdfDiffuse_[uLastIndex] /
        (uLastIndex * (cdfDiffuse_[n + 1] - cdfDiffuse_[n]));
  }
}

}  // end of namespace tgir

namespace tgir {
void PhotonQuery::begin_query() { rPower_ = 0; }

void PhotonQuery::push_back(Photon const &data,
                            Photon::vector_type const &vDistance) {
  tgir::Real const rCosTheta =
      -hi::dot(data.dir(), pathVertex_.vGeometricNormal);
  if (rCosTheta <= 0) {
    return;
  }

  tgir::Vector3 const vPoint =
      vDistance + data.dir() * hi::dot(vDistance, pathVertex_.vGeometricNormal);
  tgir::Real const rDistanceSquared = hi::length_squared(vPoint);
  if (rDistanceSquared >= PhotonQuery::max_search_range()) {
    return;
  }
#if 0  // フィルタをかけるか
#if 1  // コーンフィルタか
    static tgir::Real const k = 1.1; // 円錐フィルタのパラメータ
    tgir::Real const w = 1 - std::sqrt(rDistanceSquared) / (k * PhotonQuery::max_distance());
#else
    static tgir::Real const alpha = 0.918; // ガウスフィルタのパラメータ
    static tgir::Real const beta = 1.953;
    tgir::Real const u = 1 - std::exp(-beta * rDistanceSquared / (2 * PhotonQuery::max_search_range()));
    tgir::Real const s = 1 - std::exp(-beta);
    tgir::Real const w = alpha * (1 - u / s);
#endif
    rPower_ += data.pow() * w;
#else
  rPower_ += data.pow();
#endif
}

void PhotonQuery::push_back(Importance const &data,
                            Importance::vector_type const &vDistance) {
  tgir::Real const rCosTheta =
      -hi::dot(data.dir(), pathVertex_.vGeometricNormal);
  if (rCosTheta <= 0) {
    return;
  }

  tgir::Vector3 const vPoint =
      vDistance + data.dir() * hi::dot(vDistance, pathVertex_.vGeometricNormal);
  tgir::Real const rDistanceSquared = hi::length_squared(vPoint);
  if (rDistanceSquared >= PhotonQuery::max_search_range()) {
    return;
  }
#if 0  // フィルタをかけるか
#if 1  // コーンフィルタか
    static tgir::Real const k = 1.1; // 円錐フィルタのパラメータ
    tgir::Real const w = 1 - std::sqrt(rDistanceSquared) / (k * PhotonQuery::max_distance());
#else
    static tgir::Real const alpha = 0.918; // ガウスフィルタのパラメータ
    static tgir::Real const beta = 1.953;
    tgir::Real const u = 1 - std::exp(-beta * rDistanceSquared / (2 * PhotonQuery::max_search_range()));
    tgir::Real const s = 1 - std::exp(-beta);
    tgir::Real const w = alpha * (1 - u / s);
#endif
    rPower_ += w;
#else
  rPower_ += 1;
#endif
}

void PhotonQuery::end_query() {}

tgir::Real PhotonQuery::GetEstimatedRadiance() {
#if 0  // フィルタをかけるか
#if 1  // コーンフィルタか
    static tgir::Real const k = 1.1; // 円錐フィルタのパラメータ
    static tgir::Real const n = 1.0 - 2.0 / (3.0 * k); // 円錐フィルタの正規化定数
    return rPower_ / (n * M_PI * PhotonQuery::max_search_range());
#else
    return rPower_;
#endif
#else
  return rPower_ * M_1_PI / PhotonQuery::max_search_range();
#endif
}

}  // end of namespace tgir

namespace tgir {
PhotonMapProperty::PhotonMapProperty()
    : nPhotons_(0),
      rPowerAvg_(0),
      rPowerVar_(0),
      rPowerMin_(std::numeric_limits<tgir::Real>::infinity()),
      rPowerMax_(-std::numeric_limits<tgir::Real>::infinity()) {}

void PhotonMapProperty::operator()(Photon const &photon) {
  ++nPhotons_;

  tgir::Real const rPower = photon.pow();

  rPowerAvg_ += rPower;
  rPowerVar_ += hi::square_of(rPower);

  if (rPowerMin_ > rPower) {
    rPowerMin_ = rPower;
  }

  if (rPowerMax_ < rPower) {
    rPowerMax_ = rPower;
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

}  // end of namespace tgir

namespace tgir {
void ImportanceQuery::begin_query() {
  for (std::size_t i = 1; i <= CONFIG_N_SIZE; ++i) {
    cdfDiffuse_[i] = EPSILON_IMPORTON;
  }
}

void ImportanceQuery::push_back(Importance const &data,
                                Importance::vector_type const &vDistance) {
  tgir::Real const rCosTheta = -hi::dot(data.dir(), pathVertex_.vShadingNormal);
  if (rCosTheta <= 0) {
    return;
  }

  tgir::Vector3 const vPoint =
      vDistance + data.dir() * hi::dot(vDistance, pathVertex_.vShadingNormal);
  if (hi::length_squared(vPoint) >= ImportonQuery::max_search_range()) {
    return;
  }

  // (u, v) = (cos^2\theta, \phi/2pi)
  int u = static_cast<int>(hi::square_of(rCosTheta) * CONFIG_U_SIZE);
  if (u >= CONFIG_U_SIZE) {
    u = CONFIG_U_SIZE - 1;
  }

  tgir::Real phi =
      tgir::Real(0.5) * std::atan2(-hi::dot(data.dir(), pathVertex_.vBinormal),
                                   -hi::dot(data.dir(), pathVertex_.vTangent));
  if (phi < 0) {
    phi += M_PI;
  }

  int v = static_cast<int>(phi * (CONFIG_V_SIZE * M_1_PI));
  if (phi >= CONFIG_V_SIZE) {
    phi = CONFIG_V_SIZE - 1;
  }

  ++cdfDiffuse_[1 + u * CONFIG_V_SIZE + v];
}

void ImportanceQuery::end_query() {
  for (std::size_t i = 2; i <= CONFIG_N_SIZE; ++i) {
    cdfDiffuse_[i] += cdfDiffuse_[i - 1];
  }
}

void ImportanceQuery::SetDiffusedVector(std::mt19937_64 &random) {
  tgir::Real const value = cdfDiffuse_[CONFIG_N_SIZE] *
                           std::uniform_real_distribution<tgir::Real>()(random);
  std::vector<tgir::Real>::const_iterator const it =
      std::upper_bound(cdfDiffuse_.begin(), cdfDiffuse_.end(), value);

  std::size_t const n = it - cdfDiffuse_.begin() - 1;
  std::size_t const u = n / CONFIG_V_SIZE;
  std::size_t const v = n % CONFIG_V_SIZE;

  tgir::Real const rMuU =
      (u + std::uniform_real_distribution<tgir::Real>()(random)) /
      CONFIG_U_SIZE;
  tgir::Real const rMuV =
      (v + std::uniform_real_distribution<tgir::Real>()(random)) /
      CONFIG_V_SIZE;

  tgir::Real const rCosTheta = std::sqrt(rMuU);
  tgir::Real const rSinTheta = std::sqrt(1 - rMuU);
  tgir::Real const rPhi = tgir::Real(2 * M_PI) * rMuV;

  pathVertex_.vOutgoingDirection =
      pathVertex_.vTangent * (rSinTheta * std::cos(rPhi)) +
      pathVertex_.vBinormal * (rSinTheta * std::sin(rPhi)) +
      pathVertex_.vShadingNormal * (rCosTheta);

  pathVertex_.sQuantum *=
      cdfDiffuse_[CONFIG_N_SIZE] /
      (CONFIG_N_SIZE * (cdfDiffuse_[n + 1] - cdfDiffuse_[n]));
}

}  // end of namespace tgir
