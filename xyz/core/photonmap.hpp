#ifndef XYZ_PHOTONMAP_HPP_
#define XYZ_PHOTONMAP_HPP_

//#include "core/config.hpp"
#include "../geom/pathvertex.hpp"

//#define CONFIG_U_SIZE 2 //
//#define CONFIG_V_SIZE 4 //
#define CONFIG_U_SIZE 2  //
#define CONFIG_V_SIZE 8  //
//#define CONFIG_U_SIZE 3 //
//#define CONFIG_V_SIZE 8 //
//#define CONFIG_U_SIZE 4 //
//#define CONFIG_V_SIZE 16 //

namespace xyz {
/// <summary>
/// 視点/光源からの重要度
/// </summary>
class Importance {
 public:
  // ポリシーの定義
  typedef float3_t vector_type;
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

}  // end of namespace xyz

namespace xyz {
/// <summary>
/// フォトン
/// </summary>
class Photon {
 public:
  // ポリシーの定義
  typedef float3_t vector_type;
  typedef vector_type::value_type value_type;

  inline Photon() {}
  inline vector_type const &val() const { return pos_; }
  inline std::size_t &key() { return key_; }
  inline std::size_t const &key() const { return key_; }

 public:
  // その他の実装
  inline Photon(vector_type const &pos, vector_type const &dir,
                float3_t const pow)
      : pos_(pos), dir_(dir), pow_(pow) {}

  inline vector_type const &dir() const { return dir_; }
  inline float3_t pow() const { return pow_; }

 private:
  vector_type pos_;  // 交点位置
  std::size_t key_;  // 分割軸
  vector_type dir_;  // 入射方向
  float3_t pow_;     // 放射束
};

}  // end of namespace xyz

namespace xyz {
typedef hi::array_kdtree<Importance> ImportonMap;
typedef hi::array_kdtree<Photon> PhotonMap;

}  // end of namespace xyz

namespace xyz {
/// <summary>
/// Importance Driven Photon Tracing 用のクエリ
/// </summary>
class ImportonQuery {
 public:
  inline static Importance::value_type max_search_range() {
    return hi::square_of(CONFIG_SEARCH_RADIUS);
  }

  inline ImportonQuery(PathVertex &pathVertex)
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
  PathVertex &pathVertex_;
  std::vector<float_t> cdfDiffuse_;

  HI_DISALLOW_COPY_AND_ASSIGN(ImportonQuery);
};

}  // end of namespace xyz

namespace xyz {
/// <summary>
/// Importance Driven Photon Tracing with Particle Filter 用のクエリ
/// </summary>
class ParticleQuery {
 public:
  inline static Importance::value_type max_search_range() {
    return hi::square_of(CONFIG_SEARCH_RADIUS);
  }

  inline ParticleQuery(std::vector<PathVertex> &pathVertices,
                       std::vector<PathVertex> &tempVertices)
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
  std::vector<PathVertex> &pathVertices_;
  std::vector<PathVertex> &tempVertices_;
  std::vector<float_t> cdfDiffuse_;
  std::size_t uVertexIndex_;
  std::size_t uVertexIndexOffset_;

  ParticleQuery(ParticleQuery const &);
  ParticleQuery &operator=(ParticleQuery const &);
};

}  // end of namespace xyz

namespace xyz {
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
  inline PhotonQuery(PathVertex const &pathVertex) : pathVertex_(pathVertex) {}

  void begin_query();
  void push_back(Photon const &data, Photon::vector_type const &vDistance);
  void push_back(Importance const &data,
                 Importance::vector_type const &vDistance);
  void end_query();

  float3_t GetEstimatedRadiance();

 private:
  PathVertex const &pathVertex_;
  float3_t fPower_;

  PhotonQuery(PhotonQuery const &);
  PhotonQuery &operator=(PhotonQuery const &);
};

}  // end of namespace xyz

namespace xyz {
/// <summary>
/// フォトンマップの統計的性質を計算するオペレータ
/// </summary>
struct PhotonMapProperty {
  PhotonMapProperty();

  void operator()(Photon const &);
  void operator()(Importance const &);
  void print();

 private:
  std::size_t nPhotons_;   // フォトン数
  long double rPowerAvg_;  // 平均
  long double rPowerVar_;  // 分散
  float_t rPowerMin_;      // 最小値
  float_t rPowerMax_;      // 最大値
};

}  // end of namespace xyz

namespace xyz {
/// <summary>
/// Importance Driven Path Tracing 用のクエリ
/// </summary>
class ImportanceQuery {
 public:
  inline static Importance::value_type max_search_range() {
    return hi::square_of(CONFIG_SEARCH_RADIUS);
  }

  inline ImportanceQuery(PathVertex const &pathVertex)
      : pathVertex_(pathVertex),
        cdfDiffuse_(CONFIG_N_SIZE + 1)
  //, cdfU_(CONFIG_U_SIZE + 1)
  //, cdfV_(CONFIG_U_SIZE * (CONFIG_V_SIZE+1))
  {
    cdfDiffuse_[0] = 0;
  }

  void begin_query();
  void push_back(Importance const &data,
                 Importance::vector_type const &vDistance);
  void end_query();

 public:
  float_t GetDensity(float3_t const &dir) const;

  // pathVertex_ に拡散反射ベクトルを設定する
  // @return サンプル点の確率密度
  void SetDiffusedVector(__inout class IPrimarySample &sample,
                         __out float3_t *const pVector,
                         __out float_t *const pCosTheta,
                         __out float_t *const pDensity,
                         __out float_t *const pAlpha) const;

 private:
  PathVertex const &pathVertex_;
  std::vector<float_t> cdfDiffuse_;
  // std::vector<float_t> cdfU_; // p(u)
  // std::vector<float_t> cdfV_; // p(v|u)

  ImportanceQuery(ImportanceQuery const &);
  ImportanceQuery &operator=(ImportanceQuery const &);
};

}  // end of namespace xyz

#endif XYZ_PHOTONMAP_HPP_
