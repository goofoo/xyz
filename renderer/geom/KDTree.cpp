#include "KDTree.hpp"
#include "Triangle.hpp"
#include "Intersection.hpp"

namespace {
inline bool Contains(__in tgir::Vector3 const &pnt,
                     __in tgir::Vector3 const &min,
                     __in tgir::Vector3 const &max) {
  return (min[0] <= pnt[0]) && (pnt[0] <= max[0]) && (min[1] <= pnt[1]) &&
         (pnt[1] <= max[1]) && (min[2] <= pnt[2]) && (pnt[2] <= max[2]);
}

inline bool Overlap(__in tgir::Vector3 const &a_min,
                    __in tgir::Vector3 const &a_max,
                    __in tgir::Vector3 const &b_min,
                    __in tgir::Vector3 const &b_max) {
  return (a_min[0] <= b_max[0]) && (a_max[0] > b_min[0]) &&
         (a_min[1] <= b_max[1]) && (a_max[1] > b_min[1]) &&
         (a_min[2] <= b_max[2]) && (a_max[2] > b_min[2]);
}

inline tgir::Real SurfaceArea(__in tgir::Vector3 const &box_size) {
  return box_size[0] * box_size[1] + box_size[1] * box_size[2] +
         box_size[2] * box_size[0];
}

}  // end of unnamed namespace

namespace tgir {
class AbstractKDNode {
 public:
  virtual ~AbstractKDNode() {}

  //
  // traversal
  //
  virtual tgir::Triangle const *FindIntersection(
      __in tgir::Vector3 const &vOrigin, __in tgir::Vector3 const &vDirection,
      __out tgir::Intersection *const pParam) const = 0;
};

}  // end of namespace tgir

namespace {
//
// kD-node
//
class KDNode : public tgir::AbstractKDNode {
  HI_DISALLOW_COPY_AND_ASSIGN(KDNode);

 public:
  KDNode(__in std::size_t const axis, __in tgir::Real const value);
  virtual ~KDNode();

  //
  // traversal
  //
  inline std::size_t Axis() const { return axis_; }

  tgir::Real Value() const { return value_; }

  void SetChild(__in std::size_t i, __out tgir::AbstractKDNode const *node);

  virtual tgir::Triangle const *FindIntersection(
      __in tgir::Vector3 const &vOrigin, __in tgir::Vector3 const &vDirection,
      __out tgir::Intersection *const pParam) const;

 private:
  tgir::Real DistanceToDivisionPlane(
      __in tgir::Vector3 const &vOrigin,
      __in tgir::Vector3 const &vDirection) const;

  void EnumChild(__in tgir::Vector3 const &org, __in tgir::Vector3 const &dir,
                 __out tgir::AbstractKDNode const *nearfar[]) const;

 private:
  // 分割次元
  std::size_t const axis_;

  // 分割次元の座標
  tgir::Real const value_;

  // 子ノード
  tgir::AbstractKDNode const *children_[2];
};

KDNode::KDNode(__in std::size_t const axis, __in tgir::Real const value)
    : axis_(axis), value_(value) {
  for (std::size_t i = 0; i < 2; ++i) {
    children_[i] = nullptr;
  }
}

KDNode::~KDNode() {
  for (std::size_t i = 0; i < 2; ++i) {
    delete children_[i];
    children_[i] = nullptr;
  }
}

void KDNode::SetChild(__in std::size_t i,
                      __out tgir::AbstractKDNode const *child) {
  delete children_[i];
  children_[i] = child;
}

tgir::Triangle const *KDNode::FindIntersection(
    tgir::Vector3 const &vOrigin, tgir::Vector3 const &vDirection,
    tgir::Intersection *const pParam) const {
  tgir::Real e = DistanceToDivisionPlane(vOrigin, vDirection);

  tgir::AbstractKDNode const *nearfar[2];
  EnumChild(vOrigin, vDirection, nearfar);

  if ((e > pParam->t_max()) || (e < 0))  // 範囲外 or 逆方向
  {
    return nearfar[0]->FindIntersection(vOrigin, vDirection, pParam);
  } else if (e < pParam->t_min())  // 超えてる
  {
    return nearfar[1]->FindIntersection(vOrigin, vDirection, pParam);
  } else {
    tgir::Real tmax = pParam->t_max();
    pParam->set_t_max(e);
    tgir::Triangle const *s =
        nearfar[0]->FindIntersection(vOrigin, vDirection, pParam);
    if (s) {
      return s;
    }
    pParam->set_t_min(e);
    pParam->set_t_max(tmax);
    return nearfar[1]->FindIntersection(vOrigin, vDirection, pParam);
  }
}

tgir::Real KDNode::DistanceToDivisionPlane(
    __in tgir::Vector3 const &vOrigin,
    __in tgir::Vector3 const &vDirection) const {
  return (value_ - vOrigin[axis_]) / vDirection[axis_];
}

void KDNode::EnumChild(tgir::Vector3 const &org, tgir::Vector3 const &,
                       tgir::AbstractKDNode const *nearfar[]) const {
  // [min,value), [value,max)
  int i = org[axis_] < value_;
  nearfar[0] = children_[1 - i];
  nearfar[1] = children_[i];
}

//
// kD-leaf
//
class KDLeaf : public tgir::AbstractKDNode {
 public:
  inline KDLeaf(std::vector<tgir::Triangle const *> const &set) : set_(set) {}
  virtual ~KDLeaf() {}

  virtual tgir::Triangle const *FindIntersection(
      __in tgir::Vector3 const &vOrigin, __in tgir::Vector3 const &vDirection,
      __out tgir::Intersection *const pParam) const;

 private:
  std::vector<tgir::Triangle const *> set_;
};

tgir::Triangle const *KDLeaf::FindIntersection(
    __in tgir::Vector3 const &vOrigin, __in tgir::Vector3 const &vDirection,
    __out tgir::Intersection *const pParam) const {
  tgir::Triangle const *s = nullptr;
  for (std::size_t i = 0, size = set_.size(); i < size; ++i) {
    if (set_[i]->FindIntersection(vOrigin, vDirection, pParam)) {
      s = set_[i];
    }
  }
  return s;
}

//
// kd-nullptr
//
class KDNull : public tgir::AbstractKDNode {
 public:
  inline KDNull() {}
  virtual ~KDNull() {}

  virtual tgir::Triangle *FindIntersectionCW(
      __in tgir::Vector3 const &vOrigin, __in tgir::Vector3 const &vDirection,
      __out tgir::Intersection *const pParam) const {
    UNREFERENCED_PARAMETER(vOrigin);
    UNREFERENCED_PARAMETER(vDirection);
    UNREFERENCED_PARAMETER(pParam);
    return nullptr;
  }

  virtual tgir::Triangle *FindIntersection(
      __in tgir::Vector3 const &vOrigin, __in tgir::Vector3 const &vDirection,
      __out tgir::Intersection *const pParam) const {
    UNREFERENCED_PARAMETER(vOrigin);
    UNREFERENCED_PARAMETER(vDirection);
    UNREFERENCED_PARAMETER(pParam);
    return nullptr;
  }
};

}  // end of unnamed namespace

namespace tgir {
//
// kD-tree
//
KDTree::KDTree(__in std::size_t max_buket_size,
               __in std::size_t max_hierarchy_depth)
    : max_bucket_size_(max_buket_size),
      max_hierarchy_depth_(max_hierarchy_depth),
      root_(nullptr),
      min_(std::numeric_limits<tgir::Real>::infinity()),
      max_(-std::numeric_limits<tgir::Real>::infinity()) {}

KDTree::~KDTree() { Clear(); }

void KDTree::SetMaxBucketSize(__in std::size_t max_bucket_size) {
  max_bucket_size_ = max_bucket_size;
}

void KDTree::SetMaxHierarchyDepth(__in std::size_t max_hierarchy_depth) {
  max_hierarchy_depth_ = max_hierarchy_depth;
}

//[TODO]ノード数の評価をする．
// std::size_t node_count;

// builder
void KDTree::Build(__in std::vector<tgir::Triangle const *> const &set) {
  // compute bounding box
  min_ = std::numeric_limits<tgir::Real>::infinity();
  max_ = -std::numeric_limits<tgir::Real>::infinity();
  for (std::size_t i = 0, size = set.size(); i < size; ++i) {
    set[i]->min(min_);
    set[i]->max(max_);
  }
  max_ += (max_ - min_) *
          tgir::Real(1e-5);  // 10mのスケールで1mmぐらいの幅だけ拡張する
  // つまりここでの max は正確には upper bound を表している

  Clear();

  // std::tcerr << "buket";
  // node_count = 0;
  root_ = Build(set, min_, max_, 0);
}

AbstractKDNode *KDTree::Build(
    __in std::vector<tgir::Triangle const *> const &set,
    __in tgir::Vector3 const &min, __in tgir::Vector3 const &max,
    __in std::size_t depth) {
  if ((set.size() <= max_bucket_size_) || (depth > max_hierarchy_depth_)) {
    // std::tcerr << set.size();

    if (set.empty()) {
      return new KDNull();
    } else {
      return new KDLeaf(set);
    }
  }
  depth++;

  tgir::Vector3 const range(max - min);

  std::size_t axis = ~0U;  // 分割次元
  tgir::Real value =
      std::numeric_limits<tgir::Real>::signaling_NaN();  // 分割位置

  // surface area heuristic
  tgir::Real min_cost = std::numeric_limits<tgir::Real>::max();

  // event
  struct open_close_event {
    bool operator<(open_close_event const &e) const {
      // value が等しい場合 open Event < close Event とする
      return (value < e.value) ? true : (value > e.value) ? false
                                                          : type < e.type;
    }

    bool type;  // event type: {false: open, true: close}
    tgir::Real value;
  };

#define EVENT_UPDATE(v)                                            \
  {                                                                \
    tgir::Real const r = v - min[k];                               \
    tgir::Real const sa_left = a + b * r;                          \
    tgir::Real const sa_right = a + b * (range[k] - r);            \
    tgir::Real const cost = n_left * sa_left + n_right * sa_right; \
    if (min_cost > cost) {                                         \
      min_cost = cost;                                             \
      axis = k;                                                    \
      value = v;                                                   \
    }                                                              \
  }

  std::vector<open_close_event> event_list(set.size() * 2);
  for (std::size_t k = 0; k < 3; ++k) {
    // init event list
    for (std::size_t i = 0, size = set.size(); i < size; ++i) {
      event_list[i].type = false;
      event_list[i].value = set[i]->min(k);  // 最小値
      event_list[i + size].type = true;
      event_list[i + size].value = set[i]->max(k) + EPSILON;  // 最大値
    }

    // sort event list
    std::sort(event_list.begin(), event_list.end());

    // compute min cost
    tgir::Real const a = range[(k + 1) % 3] * range[(k + 2) % 3];
    tgir::Real const b = range[(k + 1) % 3] + range[(k + 2) % 3];

    std::size_t n_left = 0;
    std::size_t n_right = set.size();
    {
      tgir::Real const v = event_list[0].value - EPSILON;
      EVENT_UPDATE(v);
    }
    for (std::size_t i = 0, size = event_list.size(); i < size; ++i) {
      if (event_list[i].type) {
        n_right--;  // close event
      } else {
        n_left++;  // open event
      }

      tgir::Real v = event_list[i].value;
      if ((i + 1) < size) {
        if (v == event_list[i + 1].value) {
          continue;  // do not evalute
        }
        // v += (event_list[i+1].value - v) * tgir::Real(0.5);
      } else {
        v += EPSILON;
      }

      EVENT_UPDATE(v);
    }
  }
#undef EVENT_UPDATE

  // 左右に分類
  std::vector<tgir::Triangle const *> set_left;
  std::vector<tgir::Triangle const *> set_right;

  for (std::size_t i = 0, size = set.size(); i < size; ++i) {
    switch (set[i]->Overlap(axis, value)) {
    case 1:
      set_left.push_back(set[i]);
      break;  // 左
    case 2:
      set_right.push_back(set[i]);
      break;  // 右
    case 3:
      set_left.push_back(set[i]);
      set_right.push_back(set[i]);
      break;  // 両方
    default:
      assert(!"bad Triangle::Overlap");
      break;
    }
  }

  // set.swap(std::vector<tgir::Triangle const *>()); // clear & trim

  // create nodes
  KDNode *node = new KDNode(axis, value);
  {
    tgir::Vector3 max_left(max);
    max_left[axis] = value;
    node->SetChild(0, Build(set_left, min, max_left, depth));
  }
  {
    tgir::Vector3 min_right(min);
    min_right[axis] = value;
    node->SetChild(1, Build(set_right, min_right, max, depth));
  }

  return node;
}

void KDTree::Clear() {
  if (root_) {
    delete root_;
    root_ = nullptr;
  }
}

tgir::Triangle const *KDTree::FindIntersection(
    __in tgir::Vector3 const &vOrigin, __in tgir::Vector3 const &vDirection,
    __out tgir::Intersection *const pParam) const {
  pParam->Init();
  return (root_) ? root_->FindIntersection(vOrigin, vDirection, pParam)
                 : nullptr;
}

}  // end of namespace tgir
