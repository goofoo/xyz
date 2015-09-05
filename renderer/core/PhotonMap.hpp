#ifndef __TGIR_CORE_PHOTONMAP_HPP__
#define __TGIR_CORE_PHOTONMAP_HPP__

#include "core/config.hpp"
#include "geom/PathVertex.hpp"

namespace tgir {
/// <summary>
/// 視点からの重要度
/// </summary>
class Importance {
 public:
  // ポリシーの定義
  typedef tgir::Vector3 vector_type;
  typedef vector_type::value_type value_type;

  inline Importance() {}
  inline vector_type const &val() const { return pos_; }
  inline std::size_t &key() { return key_; }
  inline std::size_t const &key() const { return key_; }

 public:
  // その他の実装
  inline Importance(vector_type const &pos, vector_type const &dir)
      : pos_(pos), dir_(dir) {}

  inline vector_type const &dir() const { return dir_; }

 private:
  vector_type pos_;  // 交点位置
  std::size_t key_;  // 分割軸
  vector_type dir_;  // 入射方向
};

}  // end of namespace tgir

namespace tgir {
/// <summary>
/// フォトン
/// </summary>
class Photon {
 public:
  // ポリシーの定義
  typedef tgir::Vector3 vector_type;
  typedef vector_type::value_type value_type;

  inline Photon() {}
  inline vector_type const &val() const { return pos_; }
  inline std::size_t &key() { return key_; }
  inline std::size_t const &key() const { return key_; }

 public:
  // その他の実装
  inline Photon(vector_type const &pos, vector_type const &dir,
                tgir::Real const pow)
      : pos_(pos), dir_(dir), pow_(pow) {}

  inline vector_type const &dir() const { return dir_; }
  inline tgir::Real pow() const { return pow_; }

 private:
  vector_type pos_;  // 交点位置
  std::size_t key_;  // 分割軸
  vector_type dir_;  // 入射方向
  tgir::Real pow_;   // 放射束
};

}  // end of namespace tgir

namespace tgir {
typedef hi::array_kdtree<Importance> ImportonMap;
typedef hi::array_kdtree<Photon> PhotonMap;

}  // end of namespace tgir

namespace tgir {
/// <summary>
/// Importance Driven Photon Tracing 用のクエリ
/// </summary>
class ImportonQuery {
 public:
  inline static Importance::value_type max_search_range() {
    return hi::square_of(CONFIG_SEARCH_RADIUS);
  }

  inline ImportonQuery(tgir::PathVertex &pathVertex)
      : pathVertex_(pathVertex), cdfDiffuse_(CONFIG_N_SIZE + 1) {
    cdfDiffuse_[0] = 0;
  }

  void begin_query();
  void push_back(Importance const &data,
                 Importance::vector_type const &vDistance);
  void end_query();

  // pathVertex_ に拡散反射ベクトルを設定する
  void SetDiffusedVector(std::mt19937_64 &random);

 private:
  tgir::PathVertex &pathVertex_;
  std::vector<tgir::Real> cdfDiffuse_;

  ImportonQuery(ImportonQuery const &);
  ImportonQuery &operator=(ImportonQuery const &);
};

}  // end of namespace tgir

namespace tgir {
/// <summary>
/// Importance Driven Photon Tracing with Particle Filter 用のクエリ
/// </summary>
class ParticleQuery {
 public:
  inline static Importance::value_type max_search_range() {
    return hi::square_of(CONFIG_SEARCH_RADIUS);
  }

  inline ParticleQuery(std::vector<tgir::PathVertex> &pathVertices,
                       std::vector<tgir::PathVertex> &tempVertices)
      : pathVertices_(pathVertices),
        tempVertices_(tempVertices),
        cdfDiffuse_(CONFIG_M_SIZE + 1) {
    cdfDiffuse_[0] = 0;
  }

  void begin_query();
  void push_back(Importance const &data,
                 Importance::vector_type const &vDistance);
  void end_query();

  inline void SetVertexIndex(std::size_t uVertexIndex,
                             std::size_t uProcessIndex) {
    uVertexIndex_ = uVertexIndex;
    uVertexIndexOffset_ = uProcessIndex * CONFIG_N_SIZE + 1;
  }

  void SetDiffusedVector(std::size_t nRayCount, std::mt19937_64 &random);

 private:
  std::vector<tgir::PathVertex> &pathVertices_;
  std::vector<tgir::PathVertex> &tempVertices_;
  std::vector<tgir::Real> cdfDiffuse_;
  std::size_t uVertexIndex_;
  std::size_t uVertexIndexOffset_;

  ParticleQuery(ParticleQuery const &);
  ParticleQuery &operator=(ParticleQuery const &);
};

}  // end of namespace tgir

namespace tgir {
/// <summary>
/// Photon もしくは Importance の直接可視化用のクエリ
/// </summary>
class PhotonQuery {
 public:
  inline static Importance::value_type max_distance() {
    return CONFIG_SEARCH_RADIUS;
  }

  inline static Importance::value_type max_search_range() {
    return hi::square_of(max_distance());
  }

 public:
  inline PhotonQuery(tgir::PathVertex const &pathVertex)
      : pathVertex_(pathVertex) {}

  void begin_query();
  void push_back(Photon const &data, Photon::vector_type const &vDistance);
  void push_back(Importance const &data,
                 Importance::vector_type const &vDistance);
  void end_query();

  tgir::Real GetEstimatedRadiance();

 private:
  tgir::PathVertex const &pathVertex_;
  tgir::Real rPower_;

  PhotonQuery(PhotonQuery const &);
  PhotonQuery &operator=(PhotonQuery const &);
};

}  // end of namespace tgir

namespace tgir {
/// <summary>
/// フォトンマップの統計的性質を計算するオペレータ
/// </summary>
struct PhotonMapProperty {
  PhotonMapProperty();

  void operator()(Photon const &);
  void print();

 private:
  std::size_t nPhotons_;   // フォトン数
  long double rPowerAvg_;  // 平均
  long double rPowerVar_;  // 分散
  tgir::Real rPowerMin_;   // 最小値
  tgir::Real rPowerMax_;   // 最大値
};

}  // end of namespace tgir

namespace tgir {
/// <summary>
/// Importance Driven Path Tracing 用のクエリ
/// </summary>
class ImportanceQuery {
 public:
  inline static Importance::value_type max_search_range() {
    return hi::square_of(CONFIG_SEARCH_RADIUS);
  }

  inline ImportanceQuery(tgir::PathVertex &pathVertex)
      : pathVertex_(pathVertex), cdfDiffuse_(CONFIG_N_SIZE + 1) {
    cdfDiffuse_[0] = 0;
  }

  void begin_query();
  void push_back(Importance const &data,
                 Importance::vector_type const &vDistance);
  void end_query();

  // pathVertex_ に拡散反射ベクトルを設定する
  void SetDiffusedVector(std::mt19937_64 &random);

 private:
  tgir::PathVertex &pathVertex_;
  std::vector<tgir::Real> cdfDiffuse_;

  ImportanceQuery(ImportanceQuery const &);
  ImportanceQuery &operator=(ImportanceQuery const &);
};

}  // end of namespace tgir

#endif __TGIR_CORE_PHOTONMAP_HPP__
